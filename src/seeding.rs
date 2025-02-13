use crate::constants::MID_BASE_THRESHOLD_READ;
use crate::constants::MID_BASE_THRESHOLD_INITIAL;
use crate::types::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use std::collections::VecDeque;
use bio_seq::seq::Seq;
use bio_seq::prelude::*;

//create new alias kmer = u64
pub type Kmer64 = u64;
pub type Kmer32 = u32;
pub type KmerHash64 = u64;
pub type KmerHash32 = u32;

#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

#[inline]
pub fn rev_hash_64(hashed_key: u64) -> u64 {
    let mut key = hashed_key;

    // Invert h_key = h_key.wrapping_add(h_key << 31)
    let mut tmp: u64 = key.wrapping_sub(key << 31);
    key = key.wrapping_sub(tmp << 31);

    // Invert h_key = h_key ^ h_key >> 28;
    tmp = key ^ key >> 28;
    key = key ^ tmp >> 28;

    // Invert h_key = h_key.wrapping_add(h_key << 2).wrapping_add(h_key << 4)
    key = key.wrapping_mul(14933078535860113213u64);

    // Invert h_key = h_key ^ h_key >> 14;
    tmp = key ^ key >> 14;
    tmp = key ^ tmp >> 14;
    tmp = key ^ tmp >> 14;
    key = key ^ tmp >> 14;

    // Invert h_key = h_key.wrapping_add(h_key << 3).wrapping_add(h_key << 8)
    key = key.wrapping_mul(15244667743933553977u64);

    // Invert h_key = h_key ^ h_key >> 24
    tmp = key ^ key >> 24;
    key = key ^ tmp >> 24;

    // Invert h_key = (!h_key).wrapping_add(h_key << 21)
    tmp = !key;
    tmp = !(key.wrapping_sub(tmp << 21));
    tmp = !(key.wrapping_sub(tmp << 21));
    key = !(key.wrapping_sub(tmp << 21));

    key
}

pub fn decode(byte: u64) -> u8 {
    if byte == 0 {
        return b'A';
    } else if byte == 1 {
        return b'C';
    } else if byte == 2 {
        return b'G';
    } else if byte == 3 {
        return b'T';
    } else {
        panic!("decoding failed")
    }
}
pub fn print_string(kmer: u64, k: usize) {
    let mut bytes = vec![];
    let mask = 3;
    for i in 0..k {
        let val = kmer >> 2 * i;
        let val = val & mask;
        bytes.push(decode(val));
    }
    dbg!(std::str::from_utf8(&bytes.into_iter().rev().collect::<Vec<u8>>()).unwrap());
}
#[inline]
fn position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn minimizer_seeds_positions(
    string: &[u8],
    kmer_vec: &mut Vec<u64>,
    positions: &mut Vec<u64>,
    w: usize,
    k: usize,
) {

    if string.len() < k + w - 1 {
        return;
    }

    let mut rolling_kmer_f: Kmer64 = 0;
    let mut rolling_kmer_r: Kmer64 = 0;
    let mut canonical_kmer: Kmer64 = 0;

    let reverse_shift_dist = 2 * (k - 1);
    let max_mask = Kmer64::MAX >> (std::mem::size_of::<Kmer64>() * 8 - 2 * k);
    let rev_mask = !(3 << (2 * k - 2));
    let len = string.len();

    let rolling_window = &mut vec![u64::MAX; w];

    // populate the bit representation of the first kmer
    for i in 0..k + w - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as Kmer64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_r >>= 2;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;

        if i >= k - 1 {
            let canonical = rolling_kmer_f < rolling_kmer_r;
            canonical_kmer = if canonical {
                rolling_kmer_f
            } else {
                rolling_kmer_r
            };
            let hash = mm_hash64(canonical_kmer);
            rolling_window[i + 1 - k] = hash;
        }
    }

    let mut min_pos = position_min(rolling_window).unwrap();
    let mut min_val = rolling_window[min_pos];
    kmer_vec.push(canonical_kmer);
    positions.push(min_pos as u64);

    for i in k + w - 1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as Kmer64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= max_mask;
        rolling_kmer_r >>= 2;
        rolling_kmer_r &= rev_mask;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;

        let canonical = rolling_kmer_f < rolling_kmer_r;
        let canonical_kmer = if canonical {
            rolling_kmer_f
        } else {
            rolling_kmer_r
        };

        let hash = mm_hash64(canonical_kmer);
        let kmer_pos_global = i + 1 - k;
        rolling_window[kmer_pos_global % w] = hash;

        if hash < min_val{
            min_val = hash;
            let min_pos_global = i - k + 1;
            min_pos = min_pos_global % w;
            kmer_vec.push(hash);
            positions.push(min_pos_global as u64);
        }

        else if min_pos == (i - k + 1) % w {
            min_pos = position_min(rolling_window).unwrap();
            min_val = rolling_window[min_pos];
            let offset = (((i - k + 1) % w) as i64 - min_pos as i64).rem_euclid(w as i64);
            let min_pos_global = i - k + 1 - offset as usize;
            positions.push(min_pos_global as u64);
            kmer_vec.push(min_val);
        }
    }
}


pub fn fmh_seeds(
    string: &[u8],
    kmer_vec: &mut Vec<u64>,
    c: usize,
    k: usize
) {
    type MarkerBits = u64;
    if string.len() < k {
        return;
    }

    let marker_k = k;
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;

    let threshold_marker = u64::MAX / (c as u64);
    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }
    for i in marker_k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //

        let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
        let canonical_kmer_marker = if canonical_marker {
            rolling_kmer_f_marker
        } else {
            rolling_kmer_r_marker
        };
        let hash_marker = mm_hash64(canonical_kmer_marker);

        if hash_marker < threshold_marker {
            kmer_vec.push(hash_marker as u64);
        }
    }
}

pub fn fmh_seeds_positions(
    string: &[u8],
    kmer_vec: &mut Vec<u64>,
    positions: &mut Vec<u64>,
    c: usize,
    k: usize,
) {
    type MarkerBits = u64;
    if string.len() < k {
        return;
    }

    let marker_k = k;
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;

    let threshold_marker = u64::MAX / (c as u64);
    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }
    for i in marker_k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //

        let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
        let canonical_kmer_marker = if canonical_marker {
            rolling_kmer_f_marker
        } else {
            rolling_kmer_r_marker
        };
        let hash_marker = mm_hash64(canonical_kmer_marker);

        if hash_marker < threshold_marker {
            kmer_vec.push(canonical_kmer_marker as u64);
            positions.push((i + 1 - k) as u64);
        }
    }
}


pub fn get_twin_read_syncmer(
    string: Vec<u8>,
    qualities: Option<Vec<u8>>,
    k: usize,
    c: usize,
    snpmer_set: &FxHashSet<u64>,
    id: String,
) -> Option<TwinRead> {
    let mut snpmers_in_read = vec![];
    let mut minimizers_in_read = vec![];
    let mut dedup_snpmers = FxHashMap::default();
    let marker_k = k;

    type MarkerBits = u64;
    if string.len() < k {
        return None;
    }

    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let split_mask = !(3 << (k-1));
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    let mid_k = k / 2;

    // New syncmer-related variables
    let s = k - c + 1;  // length of syncmers
    let s_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * s);
    let s_rev_mask = !(3 << (2 * s - 2));
    let s_reverse_shift_dist = 2 * (s - 1);
    
    let mut s_mer_hashes = VecDeque::with_capacity(k - s + 1);
    let mut rolling_s_mer_f: MarkerBits = 0;
    let mut rolling_s_mer_r: MarkerBits = 0;

    let mut read_with_all_equal_qualities = false;
    if let Some(qualities) = qualities.as_ref(){
        //Ensure that not all qualities are the same value. If they are, possibly it is an old pacbio run... ignore them
        let mut q_iter = qualities.iter();
        let first_q = q_iter.next().unwrap();
        if q_iter.all(|q| q == first_q){
            read_with_all_equal_qualities = true;
        }
    }

    // Initialize first k-1 bases for k-mer
    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;

        // Also initialize s-mer if within first s-1 bases
        if i < s - 1 {
            rolling_s_mer_f <<= 2;
            rolling_s_mer_f |= nuc_f;
            rolling_s_mer_r >>= 2;
            rolling_s_mer_r |= nuc_r << s_reverse_shift_dist;
        }
    }

    for i in marker_k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;

        // Update k-mers
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;

        let split_f = rolling_kmer_f_marker & split_mask;
        let split_r = rolling_kmer_r_marker & split_mask;
    
        let canonical_marker = split_f < split_r;
        let canonical_kmer_marker = if canonical_marker {
            rolling_kmer_f_marker
        } else {
            rolling_kmer_r_marker
        };

        // Update s-mers
        rolling_s_mer_f <<= 2;
        rolling_s_mer_f |= nuc_f;
        rolling_s_mer_f &= s_mask;
        
        rolling_s_mer_r >>= 2;
        rolling_s_mer_r &= s_rev_mask;
        rolling_s_mer_r |= nuc_r << s_reverse_shift_dist;

        // Get canonical s-mer and its hash
        let canonical_s_mer = if rolling_s_mer_f < rolling_s_mer_r {
            rolling_s_mer_f
        } else {
            rolling_s_mer_r
        };

        let hash = mm_hash64(canonical_s_mer);
        
        // Add to our window of s-mer hashes
        s_mer_hashes.push_back(hash);
        if s_mer_hashes.len() > k - s + 1 {
            s_mer_hashes.pop_front();
        }
        
        // Check SNPmer
        if snpmer_set.contains(&canonical_kmer_marker) {
            let mid_base_qval = if let Some(qualities) = qualities.as_ref() {
                let mid = i + 1 + mid_k - k;
                qualities[mid] - 33
            } else {
                60
            };
            
            if mid_base_qval > MID_BASE_THRESHOLD_READ || read_with_all_equal_qualities {
                snpmers_in_read.push((i + 1 - k, canonical_kmer_marker));
            }
            *dedup_snpmers.entry(canonical_kmer_marker & split_mask).or_insert(0) += 1;
        } 
        // Check for minimizer using syncmer method
       if i >= k - 1 && s_mer_hashes.len() == k - s + 1 {
            let middle_idx = (k - s) / 2;
            let middle_hash = s_mer_hashes[middle_idx];
            
            // Check if middle s-mer has minimum hash
            if s_mer_hashes.iter().all(|h| *h >= middle_hash) {
                minimizers_in_read.push((i + 1 - k, canonical_kmer_marker));
            }
        }
    }

    let mut no_dup_snpmers_in_read = vec![];
    for (pos, kmer) in snpmers_in_read.iter_mut() {
        if dedup_snpmers[&(*kmer & split_mask)] == 1 {
            no_dup_snpmers_in_read.push((*pos, *kmer));
        }
    }

    let seq_id;
    if read_with_all_equal_qualities{
        seq_id = None;
    }
    else{
        seq_id = estimate_sequence_identity(qualities.as_ref());
    }

    let mut qual_seq : Option<Seq<QualCompact3>> = None;
    if let Some(qualities) = qualities{
        qual_seq = Some(qualities.try_into().unwrap());
    }

    Some(TwinRead{
        snpmers: no_dup_snpmers_in_read,
        minimizers: minimizers_in_read,
        base_id: id.clone(),
        id,
        k: k as u8,
        base_length: len,
        dna_seq: string.try_into().unwrap(),
        qual_seq: qual_seq,
        est_id: seq_id,
        outer: false,
        median_depth: None,
        min_depth_multi: None,
        split_chimera: false,
        split_start: 0,
        snpmer_id_threshold: None,
    })

}

pub fn get_twin_read(
    string: Vec<u8>,
    qualities: Option<Vec<u8>>,
    k: usize,
    c: usize,
    snpmer_set: &FxHashSet<u64>,
    id: String,
) -> Option<TwinRead> {

    let mut snpmers_in_read = vec![];
    let mut minimizers_in_read = vec![];
    let mut dedup_snpmers = FxHashMap::default();
    let marker_k = k;

    type MarkerBits = u64;
    if string.len() < k {
        return None;
    }

    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);

    let split_mask = !(3 << (k-1));
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    let threshold = u64::MAX / (c as u64);
    let mid_k = k / 2;

    let mut read_with_all_equal_qualities = false;
    if let Some(qualities) = qualities.as_ref(){
        //Ensure that not all qualities are the same value. If they are, possibly it is an old pacbio run... ignore them
        let mut q_iter = qualities.iter();
        let first_q = q_iter.next().unwrap();
        if q_iter.all(|q| q == first_q){
            read_with_all_equal_qualities = true;
        }
    }

    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }

    for i in marker_k-1..len {

        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;

        let split_f = rolling_kmer_f_marker & split_mask;
        let split_r = rolling_kmer_r_marker & split_mask;
    
        let canonical_marker = split_f < split_r;
        let canonical_kmer_marker; 
        if canonical_marker {
            canonical_kmer_marker = rolling_kmer_f_marker;
        } else {
            canonical_kmer_marker = rolling_kmer_r_marker;
        };
        
        if snpmer_set.contains(&canonical_kmer_marker){
            //Estimate mid base quality
            let mid_base_qval;
            if let Some(qualities) = qualities.as_ref(){
                // --xxoxx
                // pos = 2, k = 5, i = 6, mid_pos = 4
                // We want mid = pos + k/2
                // So mid = i - k + 1 + k/2
                // The middle quality val will be at k/2 + i. 
                let mid = i + 1 + mid_k - k;
                mid_base_qval = qualities[mid] - 33;
            }
            else{
                mid_base_qval = 60;
            }
            if mid_base_qval > MID_BASE_THRESHOLD_READ || read_with_all_equal_qualities{
                snpmers_in_read.push((i + 1 - k, canonical_kmer_marker));
            }
            *dedup_snpmers.entry(canonical_kmer_marker & split_mask).or_insert(0) += 1;
        }
        if mm_hash64(canonical_kmer_marker) < threshold {
            minimizers_in_read.push((i + 1 - k, canonical_kmer_marker));
        }
    }

    let mut no_dup_snpmers_in_read = vec![];
    for (pos, kmer) in snpmers_in_read.iter_mut(){
        if dedup_snpmers[&(*kmer & split_mask)] == 1{
            no_dup_snpmers_in_read.push((*pos, *kmer));
        }
    }

    let seq_id;
    if read_with_all_equal_qualities {
        seq_id = None;
    }
    else{
        seq_id = estimate_sequence_identity(qualities.as_ref());
    }

    let mut qual_seq = None;
    if let Some(qualities) = qualities{
        qual_seq = Some(qualities.try_into().unwrap());
    }

    no_dup_snpmers_in_read.shrink_to_fit();
    minimizers_in_read.shrink_to_fit();

    return Some(TwinRead{
        snpmers: no_dup_snpmers_in_read,
        minimizers: minimizers_in_read,
        base_id: id.clone(),
        id,
        k: k as u8,
        base_length: len,
        dna_seq: string.try_into().unwrap(),
        qual_seq,
        est_id: seq_id,
        outer: false,
        median_depth: None,
        min_depth_multi: None,
        split_chimera: false,
        split_start: 0,
        snpmer_id_threshold: None,
    });

}

fn estimate_sequence_identity(qualities: Option<&Vec<u8>>) -> Option<f64> {
    if qualities.is_none() {
        return None;
    }
    let mut sum = 0.0;
    let mut count = 0;
    for q in qualities.unwrap() {
        let q = (*q - 33) as f64;
        let p = 10.0f64.powf(-q / 10.0);
        sum += p;
        count += 1;
    }
    Some(100. - (sum / count as f64 * 100.))
}

pub fn split_kmer_mid(
    string: Vec<u8>,
    qualities: Option<Vec<u8>>,
    k: usize
) -> Vec<u64>{
    type MarkerBits = u64;
    if string.len() < k {
        return vec![];
    }
    let mut split_kmers = Vec::with_capacity(string.len() - k + 3);

    let marker_k = k;
    if marker_k % 2 != 1 || k > 31{
        panic!("k must be odd and <= 31");
    }
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);

    //split representation 11|11|11|00|11|11|11 for k = 6 and marker_k = 7
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let split_mask = !(3 << (k-1));
    let _split_mask_extract = !split_mask;
    let len = string.len();
    let mid_k = k / 2;
    let mut positions_to_skip = FxHashSet::default();
    if let Some(qualities) = qualities.as_ref(){
        //Ensure that not all qualities are the same value. If they are, possibly it is an old pacbio run... ignore them
        let mut q_iter = qualities.iter();
        let first_q = q_iter.next().unwrap();
        if !q_iter.all(|q| q == first_q){
            for i in marker_k-1..qualities.len(){
                let mid_pos = i + 1 + mid_k - k;
                if qualities[mid_pos] - 33 < MID_BASE_THRESHOLD_INITIAL{
                    positions_to_skip.insert(i);
                }
            }
        }
        
    }

    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }

    for i in marker_k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;

        let split_f = rolling_kmer_f_marker & split_mask;
        let split_r = rolling_kmer_r_marker & split_mask;

        //Palindromes can mess things up because the middle base
        //is automatically a SNPmer. 
        if split_f == split_r{
            continue;
        }

        // Skip low-identity mid bases
        if positions_to_skip.contains(&i){
            continue;
        }

        let canonical_marker = split_f < split_r;
        let canonical_kmer_marker; 
        //let mid_base; 
        if canonical_marker {
            canonical_kmer_marker = rolling_kmer_f_marker;
            //mid_base = (rolling_kmer_f_marker & split_mask_extract) >> (k-1) as u64;
        } else {
            canonical_kmer_marker = rolling_kmer_r_marker;
            //mid_base = (rolling_kmer_r_marker & split_mask_extract) >> (k-1) as u64;
        };
        let final_marked_kmer = canonical_kmer_marker | ((canonical_marker as u64) << (63));
        split_kmers.push(final_marked_kmer);
    }

    return split_kmers;
}
