use crate::{constants::*, graph::GraphEdge, types::*};
use bio_seq::prelude::*;
use std::panic;
use crate::unitig::*;

impl UnitigGraph{
    pub fn get_safe_edges_from_cov_threshold(&self, initial_edges: &[EdgeIndex], safety_cov_edge_ratio: Option<f64>) -> Vec<usize> {
        if initial_edges.len() == 0 {
            return vec![];
        }
        if let Some(safety_cov_edge_ratio_val) = safety_cov_edge_ratio {
            let mut cov_ratios = vec![];
            for edge_id in initial_edges{
                let edge = &self.edges[*edge_id].as_ref().unwrap();
                let uni1 = &self.nodes[&edge.from_unitig];
                let uni2 = &self.nodes[&edge.to_unitig];
                let pseudocov_ratio = pseudocount_cov_multi(uni1.min_read_depth_multi.unwrap(), uni2.min_read_depth_multi.unwrap());
                cov_ratios.push(pseudocov_ratio);
            }
            let best_ratio = *cov_ratios.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            let normalized_ratios = cov_ratios.into_iter().map(|x| x / best_ratio).collect::<Vec<_>>();
            return normalized_ratios.into_iter().enumerate().filter(|x| x.1 < safety_cov_edge_ratio_val).map(|x| initial_edges[x.0]).collect::<Vec<_>>();
        }
        else{
            return initial_edges.to_vec();
        }
    }

}

pub fn overhang_and_overlap_list(node: &UnitigNode) -> (Vec<usize>, Vec<usize>) {
    let mut overhangs = Vec::new();
    let mut overlaps = Vec::new();
    for i in 0..node.read_indices_ori.len() - 1 {
        let node1 = node.read_indices_ori[i].0;
        let node2 = node.read_indices_ori[i + 1].0;
        let internal_edge = &node.internal_overlaps[i];
        if node1 == internal_edge.node1 && node2 == internal_edge.node2 {
            overhangs.push(internal_edge.hang1);
            overhangs.push(internal_edge.hang2);
            overlaps.push(internal_edge.overlap1_len);
            overlaps.push(internal_edge.overlap2_len);
        } else if node1 == internal_edge.node2 && node2 == internal_edge.node1 {
            overhangs.push(internal_edge.hang2);
            overhangs.push(internal_edge.hang1);
            overlaps.push(internal_edge.overlap2_len);
            overlaps.push(internal_edge.overlap1_len);
        } else {
            dbg!(node1, node2, internal_edge.node1, internal_edge.node2, i);
            panic!("Internal overlap does not match read indices");
        }
    }
    return (overhangs, overlaps);
}

pub fn get_base_info_overlaps(
    left_cut: usize,
    right_cut: usize,
    node: &UnitigNode,
    reads: &[TwinRead],
    dna_seq_info: bool,
) -> BaseInfo {
    let mut length = 0;
    let mut base_seq = Seq::new();
    let mut ranges = vec![];
    let mut carryover = left_cut;
    let (overhangs, internal_ol_len) = overhang_and_overlap_list(node);
    let mut hang_ind = 0;
    for (i, indori) in node.read_indices_ori.iter().enumerate() {
        let ind = indori.0;
        let ori = indori.1;
        let mut range;
        if node.read_indices_ori.len() == 1 {
            let end;
            if right_cut > reads[ind].base_length {
                end = 0;
            } else {
                end = reads[ind].base_length - right_cut;
            }
            range = (carryover, end);
        } else if i == 0 {
            let mut end = reads[ind].base_length as i64 + 1
                - internal_ol_len[0] as i64
                - overhangs[0] as i64;
                //- overhangs[1] as i64; TODO check
            if end < 0 {
                end = 0;
            }
            range = (carryover, end as usize);
            hang_ind += 2;
        } else if i == node.read_indices_ori.len() - 1 {
            range = (
                carryover
                + overhangs[hang_ind - 1], // TODO    
                reads[ind].base_length - right_cut.min(reads[ind].base_length),
            );
            hang_ind += 2;
        } else {
            let mut end = reads[ind].base_length as i64 + 1
                - internal_ol_len[hang_ind] as i64
                - overhangs[hang_ind] as i64;
                //- overhangs[hang_ind + 1] as i64; TODO
            if end < 0 {
                end = 0;
            }
            range = (carryover + 
                overhangs[hang_ind - 1], // TODO
                end as usize);
                //, end as usize);
            hang_ind += 2;
        }
        if range.0 >= range.1 {
            if range.0 - range.1 > carryover{  
                carryover = 0;
            }
            else{
                carryover -= range.0 - range.1;
            }
            range = (0, 0);
        } else {
            carryover = 0;
        }
        length += range.1 - range.0;
        ranges.push(range);
        if dna_seq_info {
            if ori {
                base_seq.append(&reads[ind].dna_seq[range.0..range.1]);
            } else {
                base_seq.append(&reads[ind].dna_seq.to_revcomp()[range.0..range.1]);
            }
        }
    }

    let base_info = BaseInfo {
        base_seq: base_seq,
        read_positions_internal: ranges,
        length: length,
        left_cut: left_cut,
        right_cut: right_cut,
        present: true,
    };
    return base_info;
}


#[inline]
pub fn median(v: &mut [f64]) -> Option<f64> {
    if v.len() == 0 {
        return None;
    }
    let len = v.len();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if len % 2 == 0 {
        return Some((v[len / 2 - 1] + v[len / 2]) / 2.);
    }
    else{
        Some(v[len / 2])
    }

}

#[inline]
pub fn median_weight_multi(v: &[(MultiCov, usize)], quantile: f64) -> Option<MultiCov> {
    if v.len() == 0 {
        return None;
    }

    let mut multi_depths = [0.; ID_THRESHOLD_ITERS];

    for i in 0..ID_THRESHOLD_ITERS{
        let mut v = v.iter().map(|x| (x.0[i], x.1)).collect::<Vec<_>>();
        let len = v.len();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut median_ind = 0;
        let total_length = v.iter().map(|x| x.1).sum::<usize>();
        let mut curr_sum = 0;
        for i in 0..len {
            curr_sum += v[i].1;
            if curr_sum as f64 >= total_length as f64 * quantile {
                median_ind = i;
                break;
            }
        }
        multi_depths[i] = v[median_ind].0
    }

    return Some(multi_depths);
}

#[inline]
pub fn median_weight(v: &[(f64, usize)], quantile: f64) -> Option<f64> {
    let mut v = v.to_vec();
    if v.len() == 0 {
        return None;
    }

    let len = v.len();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut median_ind = 0;
    let total_length = v.iter().map(|x| x.1).sum::<usize>();
    let mut curr_sum = 0;
    for i in 0..len {
        curr_sum += v[i].1;
        if curr_sum as f64 >= total_length as f64 * quantile {
            median_ind = i;
            break;
        }
    }
    Some(v[median_ind].0)
}

pub fn quantile_dist(v1: &[([f64;ID_THRESHOLD_ITERS], usize)], v2: &[([f64;ID_THRESHOLD_ITERS], usize)]) -> Option<f64> {
    let quants = [0.25, 0.5, 0.75];
    let quantiles_v1 = quants.iter().map(|x| median_weight_multi(v1, *x)).collect::<Vec<_>>();
    let quantiles_v2 = quants.iter().map(|x| median_weight_multi(v2, *x)).collect::<Vec<_>>();
    let mut minimum_distance = vec![];
    for i in 0..quants.len() {
        for j in 0..quants.len() {
            if quantiles_v1[i].is_none() || quantiles_v2[j].is_none() {
                continue;
            }
            //minimum_distance.push((quantiles_v2[j].unwrap() - quantiles_v1[i].unwrap()).abs());
            let qvi = quantiles_v1[i].unwrap();
            let qvj = quantiles_v2[j].unwrap();
            minimum_distance.push(pseudocount_cov_multi(qvi, qvj).abs() - 1.);
        }
    }
    if minimum_distance.len() == 0 {
        return None;
    }
    Some(*minimum_distance.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap())
}

pub fn log_distribution_distance(v1: &[([f64;ID_THRESHOLD_ITERS], usize)], v2: &[([f64;ID_THRESHOLD_ITERS], usize)]) -> Option<f64> {
    let quants = [0.25, 0.5, 0.75];
    
    // Transform into logs and ratios
    let mut log_v1 = vec![];
    let mut log_v2 = vec![];

    for j in 0..2 {
        let iterable = if j == 0 { v1 } else { v2 };
        for tup in iterable.iter() {
            let mut new_tup = [0.; ID_THRESHOLD_ITERS];
            for i in 0..ID_THRESHOLD_ITERS {
                if i == 0 {
                    new_tup[i] = (tup.0[i] + PSEUDOCOUNT).ln();
                } else {
                    //new_tup[i] = ((tup.0[i] + PSEUDOCOUNT) / (tup.0[i - 1] + PSEUDOCOUNT)).ln();
                    new_tup[i] = (tup.0[i] + PSEUDOCOUNT).ln();
                }
            }
            if j == 0 {
                log_v1.push((new_tup, tup.1));
            } else {
                log_v2.push((new_tup, tup.1));
            }
        }
    }

    let larger;
    let smaller;

    if log_v1.len() >= log_v2.len(){
        larger = &log_v1;
        smaller = &log_v2;
    }
    //TODO doesn't use larger/smaller
    else{
        larger = &log_v2;
        smaller = &log_v1;
    }

    let quantiles_larger = quants.iter().map(|x| median_weight_multi(&larger, *x)).collect::<Vec<_>>();
    let quantiles_smaller = quants.iter().map(|x| median_weight_multi(&smaller, *x)).collect::<Vec<_>>();

    //median log ratios
    let median_larger_3 = &quantiles_larger[1];
    let median_smaller_3 = &quantiles_smaller[1];

    
    let mut distances = vec![];
    let mut weights = vec![];

    for i in 0..ID_THRESHOLD_ITERS{
        let d = (median_larger_3.unwrap()[i] - median_smaller_3.unwrap()[i]).abs();
        distances.push(d);

        
        let mut ratio_distribution = vec![];
        let mut count = 0;

        for (log_ratio, length) in larger.iter(){
            let index = count % smaller.len();
            let w = log_ratio[i] - smaller[index].0[i];
            ratio_distribution.push((w, *length));
            count += 1;
        }

        let iqr_weight  = median_weight(&ratio_distribution, 0.75).unwrap() - median_weight(&ratio_distribution, 0.25).unwrap();
        weights.push(iqr_weight + 1.);
    }

    //This is wrong TODO
    let weighted_distance = distances.iter().zip(weights.iter()).map(|(x, y)| x * y).sum::<f64>();

    return Some(weighted_distance);
}

pub fn log_shape_distance(
    v1: &[([f64;ID_THRESHOLD_ITERS], usize)], 
    v2: &[([f64;ID_THRESHOLD_ITERS], usize)],
) -> f64 {
    // Transform into logs and ratios
    let mut log_v1 = vec![];
    let mut log_v2 = vec![];

    for j in 0..2 {
        let iterable = if j == 0 { v1 } else { v2 };
        for tup in iterable.iter() {
            let mut new_tup = [0.; ID_THRESHOLD_ITERS];
            for i in 0..ID_THRESHOLD_ITERS {
                if i == 0 {
                    new_tup[i] = (tup.0[i] + PSEUDOCOUNT).ln();
                } else {
                    //new_tup[i] = ((tup.0[i] + PSEUDOCOUNT) / (tup.0[i - 1] + PSEUDOCOUNT)).ln();
                    new_tup[i] = (tup.0[i] + PSEUDOCOUNT).ln();
                }
            }
            if j == 0 {
                log_v1.push((new_tup, tup.1));
            } else {
                log_v2.push((new_tup, tup.1));
            }
        }
    }

    let larger;
    let smaller;

    if log_v1.len() >= log_v2.len(){
        larger = &log_v1;
        smaller = &log_v2;
    }
    //TODO doesn't use larger/smaller
    else{
        larger = &log_v2;
        smaller = &log_v1;
    }

    //median log ratios
     let mut shape_distribution_distances = vec![];
     let num_shapes = ID_THRESHOLD_ITERS * (ID_THRESHOLD_ITERS - 1) / 2;
     for i in 0..ID_THRESHOLD_ITERS{
         for j in i + 1..ID_THRESHOLD_ITERS{
             let mut shape_larger = vec![];
             let mut shape_smaller = vec![];
             for (log_ratio, _) in larger.iter(){
                 shape_larger.push(log_ratio[i] - log_ratio[j]);
             }
             for (log_ratio, _) in smaller.iter(){
                 shape_smaller.push(log_ratio[i] - log_ratio[j]);
             }

             let mut shape_distances_dist = vec![];

             for i in 0..shape_larger.len(){
                 shape_distances_dist.push((shape_larger[i] - shape_smaller[i % shape_smaller.len()]).abs());
             }

             shape_distances_dist.sort_by(|a, b| a.partial_cmp(b).unwrap());
             //dbg!(&shape_distances_dist);
             let lower_quartile = shape_distances_dist[shape_distances_dist.len() / 10];
             shape_distribution_distances.push(lower_quartile);
         }
     }

     let shape_distribution_distance;
     // distribution "shapes" have huge variance, limit this. 
    shape_distribution_distance = shape_distribution_distances.iter().sum::<f64>() / num_shapes as f64;

    return shape_distribution_distance;
}


pub fn log_distribution_distance_new(
    v1: &[([f64;ID_THRESHOLD_ITERS], usize)], 
    v2: &[([f64;ID_THRESHOLD_ITERS], usize)],
) -> Option<f64> {
    let quants = [0.25, 0.5, 0.75];
    
    // Transform into logs and ratios
    let mut log_v1 = vec![];
    let mut log_v2 = vec![];

    for j in 0..2 {
        let iterable = if j == 0 { v1 } else { v2 };
        for tup in iterable.iter() {
            let mut new_tup = [0.; ID_THRESHOLD_ITERS];
            for i in 0..ID_THRESHOLD_ITERS {
                if i == 0 {
                    new_tup[i] = (tup.0[i] + PSEUDOCOUNT).ln();
                } else {
                    //new_tup[i] = ((tup.0[i] + PSEUDOCOUNT) / (tup.0[i - 1] + PSEUDOCOUNT)).ln();
                    new_tup[i] = (tup.0[i] + PSEUDOCOUNT).ln();
                }
            }
            if j == 0 {
                log_v1.push((new_tup, tup.1));
            } else {
                log_v2.push((new_tup, tup.1));
            }
        }
    }

    let larger;
    let smaller;

    if log_v1.len() >= log_v2.len(){
        larger = &log_v1;
        smaller = &log_v2;
    }
    //TODO doesn't use larger/smaller
    else{
        larger = &log_v2;
        smaller = &log_v1;
    }

    let quantiles_larger = quants.iter().map(|x| median_weight_multi(&larger, *x)).collect::<Vec<_>>();
    let quantiles_smaller = quants.iter().map(|x| median_weight_multi(&smaller, *x)).collect::<Vec<_>>();

    //median log ratios
    let median_larger_3 = &quantiles_larger[1];
    let median_smaller_3 = &quantiles_smaller[1];

    let mut distances = vec![];
    let mut weights = vec![];

    let min_cov_threshold = (MIN_COV_READ as f64 + PSEUDOCOUNT).ln();
    let specificity_threshold = (SPECIFICITY_THRESHOLD as f64 + PSEUDOCOUNT).ln();
    assert!(specificity_threshold > min_cov_threshold);

    for i in 0..ID_THRESHOLD_ITERS{

        //Skip if both are below threshold
        if median_larger_3.unwrap()[i] < min_cov_threshold && median_smaller_3.unwrap()[i] < min_cov_threshold && i != 0{
            continue;
        }

        // Only use the highest  resolution if both are > 30 ; unlikely to be noisy.
        // if median_larger_3.unwrap()[ID_THRESHOLD_ITERS - 1] > specificity_threshold || median_smaller_3.unwrap()[ID_THRESHOLD_ITERS - 1] > specificity_threshold{
        //     if i != ID_THRESHOLD_ITERS - 1 {
        //         continue;
        //     }
        // }

        let d = (median_larger_3.unwrap()[i] - median_smaller_3.unwrap()[i]).abs();

        let upper_dist = (quantiles_larger[ID_THRESHOLD_ITERS - 1].unwrap()[i] + PSEUDOCOUNT).ln() - 
                        (quantiles_smaller[ID_THRESHOLD_ITERS - 1].unwrap()[i] + PSEUDOCOUNT).ln();

        let lower_dist = (quantiles_larger[0].unwrap()[i] + PSEUDOCOUNT).ln() -
                        (quantiles_smaller[0].unwrap()[i] + PSEUDOCOUNT).ln();
        
        let upper_lower_interval = (upper_dist - lower_dist).abs();
        distances.push(d + upper_lower_interval);

        weights.push(1.);
    }

    // // ----- SHAPE distribution distance ------  Not working so well; editing out. 
    // let mut shape_distribution_distances = vec![];
    // let num_shapes = ID_THRESHOLD_ITERS * (ID_THRESHOLD_ITERS - 1) / 2;
    // for i in 0..ID_THRESHOLD_ITERS{
    //     for j in i + 1..ID_THRESHOLD_ITERS{
    //         let mut shape_larger = vec![];
    //         let mut shape_smaller = vec![];
    //         for (log_ratio, _) in larger.iter(){
    //             shape_larger.push(log_ratio[i] - log_ratio[j]);
    //         }
    //         for (log_ratio, _) in smaller.iter(){
    //             shape_smaller.push(log_ratio[i] - log_ratio[j]);
    //         }

    //         let mut shape_distances_dist = vec![];

    //         for i in 0..shape_larger.len(){
    //             shape_distances_dist.push((shape_larger[i] - shape_smaller[i % shape_smaller.len()]).abs());
    //         }

    //         shape_distances_dist.sort_by(|a, b| a.partial_cmp(b).unwrap());
    //         let lower_quartile = shape_distances_dist[distances.len() / 4];
    //         shape_distribution_distances.push(lower_quartile);
    //     }
    // }

    // let shape_distribution_distance;
    // // distribution "shapes" have huge variance, limit this. 
    // if larger.len() < 2 {
    //     shape_distribution_distance = 0.;
    // }
    // else{
    //     shape_distribution_distance = shape_distribution_distances.iter().sum::<f64>() / num_shapes as f64;
    // }

    //let weighted_distance = distances.iter().zip(weights.iter()).map(|(x, y)| x * y).sum::<f64>();
    let weighted_distance = ID_THRESHOLD_ITERS as f64 * distances.iter().zip(weights.iter()).map(|(x, y)| x * y).min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    //return Some(weighted_distance + shape_distribution_distance);
    return Some(weighted_distance);
}

#[inline]
pub fn square_root_poisson_kl(cov1: f64, cov2: f64) -> f64 {
    let p_cov1 = cov1 + 1.;
    let p_cov2 = cov2 + 1.;
    let kl_1 = p_cov1 * (p_cov1.ln() - p_cov2.ln());
    let kl_2 = p_cov2 * (p_cov2.ln() - p_cov1.ln());
    return (kl_1 + kl_2) / (p_cov1 + p_cov2).sqrt();
}

#[inline]
pub fn pseudocount_cov(cov1: f64, cov2: f64) -> f64 {
    let pcount = PSEUDOCOUNT;
    let numerator = cov1.max(cov2) + pcount;
    let denominator = cov1.min(cov2) + pcount;
    return numerator / denominator;
}

#[inline]
pub fn pseudocount_cov_multi(cov1: [f64;ID_THRESHOLD_ITERS], cov2: [f64;ID_THRESHOLD_ITERS]) -> f64 {
    let weights = COV_MULTI_WEIGHTS;
    assert!(weights.len() == ID_THRESHOLD_ITERS);
    let mut cumulative_pcov = 0.;
    for (i,weight) in weights.iter().enumerate(){
        cumulative_pcov += pseudocount_cov(cov1[i], cov2[i]) * weight;
    }
    return cumulative_pcov;
}

pub fn rescue_id_est_nanopore_strict(edge: &UnitigEdge, c: u64) -> bool {
    let id_est = edge.edge_id_est(c as usize);

    let max_allowable_diff_snpmers = (edge.overlap.shared_snpmers / MAX_ALLOWABLE_SNPMER_ERROR_DIVIDER).max(MAX_ALLOWABLE_SNPMER_ERROR_MISC);
    // 0.05, but increased by 0.01 for each 200 SNPmers.
    let additional_fsv_hetero = (edge.overlap.shared_snpmers as f64 / MAX_ALLOWABLE_SNPMER_ERROR_DIVIDER as f64) / 100. / 4.0;
    let additional_fsv_hetero = additional_fsv_hetero.min(0.005);

    if id_est > (1.0 - 0.0050 / 100. - additional_fsv_hetero) && edge.overlap.diff_snpmers < max_allowable_diff_snpmers {
        return true;
    }
    else{
        return false;
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_weight() {
        let mut v = vec![(1.0, 1), (2.0, 2), (3.0, 3)];
        assert_eq!(median_weight(&mut v, 0.5), Some(2.0));

        let mut v = vec![(1.0, 1), (2.0, 5), (3.0, 15)];
        assert_eq!(median_weight(&mut v, 0.5), Some(3.0));

        let mut v = vec![(1.0, 1), (2.0, 5), (3.0, 1)];
        assert_eq!(median_weight(&mut v, 0.5), Some(2.0));

        let mut v = vec![];
        assert_eq!(median_weight(&mut v, 0.5), None);
    }

    #[test]
    fn test_log_dist() {
        let v1 = vec![([1.0, 1.0, 1.0], 1), ([2.0, 2.0, 2.0], 2), ([3.0, 3.0, 3.0], 3)];
        let v2 = vec![([1.0, 1.0, 1.0], 1), ([2.0, 2.0, 2.0], 2), ([3.0, 3.0, 3.0], 3)];
        assert_eq!(log_distribution_distance(&v1, &v2), Some(0.0));

        let v1 = vec![([1.0, 1.0, 1.0], 1), ([2.0, 2.0, 2.0], 2), ([3.0, 3.0, 3.0], 3)];
        let v2 = vec![([2.0, 2.0, 2.0], 1), ([2.0, 2.0, 2.0], 2), ([3.0, 3.0, 3.0], 2)];
        assert_eq!(log_distribution_distance(&v1, &v2), Some(0.0));

        // should be approximately log (2/1)
        let v1 = vec![([10., 5.0, 5.0], 2)];
        let v2 = vec![([10., 10.0, 10.0], 2)];
        let num = 5.0 + PSEUDOCOUNT;
        let denom = 10.0 + PSEUDOCOUNT;
        dbg!(log_distribution_distance(&v1, &v2).unwrap());
        assert!(log_distribution_distance(&v1, &v2).unwrap() - 2. * (num/denom).ln().abs() < 0.0001);

    }

    #[test]
    fn test_log_dist_iqrs() {
        let v1 = vec![([10.;3], 1), ([5.;3], 1), ([5.;3], 1), ([5.;3], 1), ([10.;3], 1)];
        let v2 = vec![([10.;3], 1), ([10.;3], 1), ([10.;3], 1), ([10.;3], 1), ([10.;3], 1)];
        let dist1 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist1);

        let v1 = vec![([5.;3], 1), ([5.;3], 1), ([5.;3], 1), ([5.;3], 1), ([5.;3], 1)];
        let v2 = vec![([10.;3], 1), ([10.;3], 1), ([10.;3], 1), ([10.;3], 1), ([10.;3], 1)];
        let dist2 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist2);

        //Wider distribution -> less distance
        //assert!(dist2 > dist1);
    }

    #[test]
    fn test_log_dist_sample_size_consideration() {
        let v1 = vec![([10.;3], 1), ([5.;3], 1), ([5.;3], 1), ([5.;3], 1), ([10.;3], 1)];
        let v2 = vec![([10.;3], 1), ([10.;3], 1), ([10.;3], 1), ([10.;3], 1), ([10.;3], 1)];
        let dist1 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist1);

        let v1 = vec![([10.;3], 1), ([5.;3], 1), ([5.;3], 1), ([5.;3], 1), ([10.;3], 1)];
        let v2 = vec![([10.;3], 1)];
        let dist2 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist2);

        //Wider distribution -> less distance
        assert!(dist1 == dist2);
    }

    #[test]
    fn test_log_dist_sample_size_variance_2() {
        let v1 = vec![([100.;3], 1), ([5.;3], 1), ([100.;3], 1), ([5.;3], 1), ([100.;3], 1), ([5.;3], 1), ([100.;3], 1), ([5.;3], 1), ([100.;3], 1), ([5.;3], 1)];
        let v2 = vec![([30.;3], 1), ([30.;3], 1), ([30.;3], 1), ([30.;3], 1), ([30.;3], 1)];
        let dist1 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist1);

        let v1 = vec![([5.;3], 1), ([5.;3], 1), ([5.;3], 1), ([5.;3], 1), ([5.;3], 1)];
        let v2 = vec![([10.;3], 1)];
        let dist2 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist2);

        //Wider distribution -> less distance
        assert!(dist1 > dist2);
    }

    #[test]
    fn test_log_dist_interspecies_repeats() {
        let v1 = vec![([100., 50., 20.], 1), ([20., 20., 20.,], 1), ([100., 50., 20.], 1), ([20., 20., 20.,], 1), ([100., 50., 20.], 1), ([20., 20., 20.,], 1), ([100., 50., 20.], 1), ([20., 20., 20.,], 1), ([100., 50., 20.], 1), ([20., 20., 20.,], 1)];
        let v2 = vec![([20., 20., 20.,], 1), ([20., 20., 20.,], 1), ([20., 20., 20.,], 1), ([20., 20., 20.,], 1), ([20., 20., 20.,], 1)];
        let dist1 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist1);

        //Coord 3 is not variable, low distance

        let v1 = vec![([100., 50., 1.], 1), ([90., 20., 2.,], 1), ([110., 20., 5.], 1), ([110., 20., 1.,], 1), ([100., 5., 2.], 1), ([90., 20., 2.,], 1), ([110., 20., 10.], 1), ([110., 20., 1.,], 1), ([100., 5., 2.], 1), ([90., 20., 2.,], 1)];
        let v2 = vec![([100., 50., 10.,], 1), ([100., 20., 1.,], 1)];
        let dist2 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist2);
        

        let v1 = vec![([100., 90., 1.], 1), ([90., 20., 2.,], 1), ([110., 20., 5.], 1), ([110., 20., 1.,], 1), ([100., 5., 2.], 1), ([90., 20., 2.,], 1), ([110., 20., 10.], 1), ([110., 20., 1.,], 1), ([100., 10., 2.], 1), ([90., 20., 5.,], 1)];
        let v2 = vec![([100., 23., 5.,], 1), ([100., 20., 1.,], 1)];
        let dist3 = log_distribution_distance_new(&v1, &v2).unwrap();
        dbg!(dist3);
        assert!(dist1 < 1.);
        assert!(dist2 < 1.);
        assert!(dist3 < 1.);
    }

    #[test]
    fn test_log_dist_shape() {
        let v1 = vec![([90., 90., 1.], 1), ([90., 20., 2.,], 1), ([110., 100., 5.], 1), ([20., 20., 1.,], 1), ([100., 50., 2.], 1), ([90., 20., 2.,], 1), ([110., 100., 10.], 1), ([110., 20., 1.,], 1), ([20., 10., 2.], 1), ([90., 20., 5.,], 1)];
        let mut v1_double = v1.clone();
        v1_double.extend(v1.clone());
        let v2 = vec![([90., 90., 90.,], 1), ([90., 80., 80.,], 1), ([90., 70., 70.,], 1), ([90., 60., 60.,], 1), ([90., 50., 50.,], 1)];
        let dist1 = log_shape_distance(&v1_double, &v2);
        dbg!(dist1);
        assert!(dist1 > 0.5);
    }

    #[test]
    fn test_log_dist_shape_2() {

        fn utup_to_ftup(ut: &[usize; 3]) -> [f64; 3] {
            [ut[0] as f64, ut[1] as f64, ut[2] as f64]
        }
        let dists1 = vec![
            [11,11,8], [11,11,8], [11,11,8], [11,11,8], [11,11,8], 
            [2,2,2], [2,2,2], [2,2,2], [2,2,2], [2,2,2], 
            [20,20,10], [20,20,10], [20,20,10], [20,20,10], [20,20,10], 
            [7,7,5], [7,7,5], [7,7,5], [7,7,5], [7,7,5]];

        let dists2 = vec![
            [9,9,5],[10,10,5],[6,6,6],[6,6,6],[6,6,6],[16,16,8],[5,5,3],[4,4,3],[14,14,8],[7,1,1]];

        let dists3 = vec![
            [1,1,1],[161,3,7],[35,17,7],[70,64,3],[196,192,7],[4,4,4],[5,5,2],[60,54,5],[7,7,7]
        ];

        let v1 = dists1.iter().map(|x| (utup_to_ftup(x), 1)).collect::<Vec<_>>();
        let v2 = dists2.iter().map(|x| (utup_to_ftup(x), 1)).collect::<Vec<_>>();
        let v3 = dists3.iter().map(|x| (utup_to_ftup(x), 1)).collect::<Vec<_>>();

        let dist1 = log_shape_distance(&v1, &v2);
        let dist2 = log_shape_distance(&v1, &v3);
        dbg!(dist1);
        dbg!(dist2);

        let dist1_noshape = log_distribution_distance_new(&v1, &v2).unwrap();
        let dist2_noshape = log_distribution_distance_new(&v1, &v3).unwrap();

        dbg!(dist1_noshape);
        dbg!(dist2_noshape);
    }

}
