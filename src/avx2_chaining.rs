use std::arch::x86_64::*;
use crate::types::*;
use rust_lapper::Interval;
use fxhash::FxHashSet;

#[target_feature(enable = "avx2")]
pub unsafe fn dp_anchors_simd(
    matches: &[Anchor],
    reverse: bool,
    gap_cost: i32,
    match_score: i32,
    band: usize,
) -> Vec<(i32, Vec<Anchor>)> {
    if matches.is_empty() {
        return vec![(0, Vec::new())];
    }
    let mut dp = vec![0i32; matches.len()];
    let mut prev = vec![None; matches.len()];
    let mut max_score = 0i32;
    let mut max_index = 0;
    let max_skip = 10;
    
    // Create SIMD constants as i32
    let match_score_v = _mm256_set1_epi32(match_score);
    let gap_cost_v = _mm256_set1_epi32(gap_cost);
    
    for i in 0..matches.len() {
        let start1 = matches[i].pos1 as i32;
        let start2 = matches[i].pos2 as i32;
        dp[i] = match_score;
        let mut unimproved = 0;
        let back = if i > band { i - band } else { 0 };
        
        // Process 8 elements at a time using AVX2
        let chunks = ((i - back) / 8) * 8;
        if chunks > 0 {
            for j in ((i - chunks)..i).rev().step_by(8){
                // Load 8 previous positions
                let end1_v = _mm256_set_epi32(
                    matches[j-7].pos1 as i32,
                    matches[j-6].pos1 as i32,
                    matches[j-5].pos1 as i32,
                    matches[j-4].pos1 as i32,
                    matches[j-3].pos1 as i32,
                    matches[j-2].pos1 as i32,
                    matches[j-1].pos1 as i32,
                    matches[j].pos1 as i32,
                );
                let end2_v = _mm256_set_epi32(
                    matches[j-7].pos2 as i32,
                    matches[j-6].pos2 as i32,
                    matches[j-5].pos2 as i32,
                    matches[j-4].pos2 as i32,
                    matches[j-3].pos2 as i32,
                    matches[j-2].pos2 as i32,
                    matches[j-1].pos2 as i32,
                    matches[j].pos2 as i32,
                );
                
                let start1_v = _mm256_set1_epi32(start1);
                let start2_v = _mm256_set1_epi32(start2);
                
                // Check conditions using integer comparisons
                let mask = if reverse {
                    let cond1 = _mm256_cmpgt_epi32(end1_v, start1_v);  // end1 >= start1
                    let cond2 = _mm256_cmpgt_epi32(start2_v, end2_v);  // end2 <= start2
                    _mm256_or_si256(cond1, cond2)
                } else {
                    let cond1 = _mm256_cmpgt_epi32(end1_v, start1_v);  // end1 >= start1
                    let cond2 = _mm256_cmpgt_epi32(end2_v, start2_v);  // end2 >= start2
                    _mm256_or_si256(cond1, cond2)
                };
                
                // Calculate gap penalties
                let gap1 = _mm256_sub_epi32(start1_v, end1_v);
                let gap2 = _mm256_sub_epi32(start2_v, end2_v);
                
                let gap1_abs = _mm256_abs_epi32(gap1);
                let gap2_abs = _mm256_abs_epi32(gap2);
                
                let gap_diff = _mm256_sub_epi32(gap1_abs, gap2_abs);
                let gap_diff_abs = _mm256_abs_epi32(gap_diff);
                let gap_penalty = _mm256_mullo_epi32(gap_cost_v, gap_diff_abs);
                
                // Load previous dp values
                let dp_prev = _mm256_set_epi32(
                    dp[j-7],
                    dp[j-6],
                    dp[j-5],
                    dp[j-4],
                    dp[j-3],
                    dp[j-2],
                    dp[j-1],
                    dp[j],
                );
                
                // Calculate scores
                let scores = _mm256_sub_epi32(
                    _mm256_add_epi32(dp_prev, match_score_v),
                    gap_penalty
                );
                
                // Apply mask - set masked values to i32::MIN
                let min_score = _mm256_set1_epi32(i32::MIN);
                let masked_scores = _mm256_blendv_epi8(
                    scores,
                    min_score,
                    mask,
                );
                
                // Extract and compare scores
                let mut score_array = [0i32; 8];
                _mm256_storeu_si256(score_array.as_mut_ptr() as *mut __m256i, masked_scores);
                
                for (k, &score) in score_array.iter().enumerate() {
                    if score > dp[i]{
                        dp[i] = score;
                        prev[i] = Some(j - k);
                        if score > max_score {
                            max_score = score;
                            max_index = i;
                        }
                        unimproved = 0;
                    } else if score != i32::MIN {
                        unimproved += 1;
                    }
                }
                
                if unimproved > max_skip {
                    break;
                }
            }
        }

        if unimproved > max_skip {
            continue;
        }
        
        // Handle remaining elements scalar
        for j in (back..(i - chunks)).rev() {
            let end1 = matches[j].pos1 as i32;
            let end2 = matches[j].pos2 as i32;
            if reverse {
                if end1 >= start1 || end2 <= start2 {
                    continue;
                }
            } else {
                if end1 >= start1 || end2 >= start2 {
                    continue;
                }
            }
            let gap_penalty = (start1 - end1).abs() - (start2 - end2).abs();
            let score = dp[j] + match_score - gap_cost * gap_penalty.abs();
            if score > dp[i] {
                dp[i] = score;
                prev[i] = Some(j);
                if score > max_score {
                    max_score = score;
                    max_index = i;
                }
                unimproved = 0;
            } else {
                unimproved += 1;
            }
            if unimproved > max_skip {
                break;
            }
        }
    }
    
    // Reconstruct the chain
    let mut chains = Vec::new();
    let mut used_anchors = FxHashSet::default();
    let mut reference_ranges : Vec<Interval<u32,bool>> = Vec::new();
    let mut best_indices_ordered = (0..matches.len())
        .map(|i| (dp[i], i))
        .collect::<Vec<_>>();
    best_indices_ordered.sort_unstable_by_key(|&(score, _)| -score);
    assert!(dp[max_index] == best_indices_ordered[0].0);

    for (score, best_index) in best_indices_ordered{
        if used_anchors.contains(&best_index) {
            continue;
        }
        if score < max_score * 3 / 4 {
            break;
        }

        let mut chain = Vec::new();
        let mut i = Some(best_index);
        while let Some(idx) = i {
            used_anchors.insert(idx);
            chain.push(matches[idx].clone());
            i = prev[idx];
        }

        if chain.len() < 2 {
            break;
        }

        chain.reverse();

        let l = chain.first().unwrap().pos2;
        let r = chain.last().unwrap().pos2;
        let interval = Interval{start: l.min(r), stop: l.max(r), val: true};
        if reference_ranges.iter().any(|x| {
            let intersect = x.intersect(&interval);
            intersect as f64 / (interval.stop - interval.start) as f64 > 0.25
        }) {
            continue;
        }
        chains.push((score, chain));
        reference_ranges.push(interval);
    }
    
    return chains
}