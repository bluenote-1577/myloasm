use crate::types::*;
use fxhash::FxHashSet;

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

pub fn get_twin_read(
    string: Vec<u8>,
    k: usize,
    c: usize,
    snpmer_set: &FxHashSet<u64>,
    id: String,
) -> Option<TwinRead> {

    let mut snpmers_in_read = vec![];
    let mut minimizers_in_read = vec![];
    let marker_k = k;

    type MarkerBits = u64;
    if string.len() < k {
        return None;
    }

    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);

    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    let threshold = u64::MAX / (c as u64);

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

        let split_mask = !(3 << (k-1));
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
            snpmers_in_read.push((i, canonical_kmer_marker));
        }
        else if mm_hash64(canonical_kmer_marker) < threshold {
            minimizers_in_read.push((i, canonical_kmer_marker));
        }
    }

    return Some(TwinRead{
        snpmers: snpmers_in_read,
        minimizers: minimizers_in_read,
        id,
        k: k as u8,
        base_length: len 
    });

}

pub fn split_kmer_mid(
    string: Vec<u8>,
    k: usize
) -> Vec<SplitKmerInfo>{
    let mut split_kmers = vec![];
    type MarkerBits = u64;
    if string.len() < k {
        return vec![];
    }

    let marker_k = k;
    if marker_k % 2 != 1{
        panic!("k must be odd");
    }
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);

    //split representation 11|11|11|00|11|11|11 for k = 6 and marker_k = 7
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let split_mask = !(3 << (k-1));
    let split_mask_extract = !split_mask;
    let len = string.len();

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
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //

        let canonical_marker = split_f < split_r;
        let canonical_kmer_marker; 
        let mid_base; 
        if canonical_marker {
            canonical_kmer_marker = split_f;
            mid_base = (rolling_kmer_f_marker & split_mask_extract) >> (k-1) as u64;
        } else {
            canonical_kmer_marker = split_r;
            mid_base = (rolling_kmer_r_marker & split_mask_extract) >> (k-1) as u64;
        };
        let split_kmer_info = SplitKmerInfo{
            split_kmer: canonical_kmer_marker,
            mid_base : mid_base as u8,
            canonical: canonical_marker,
            k: k as u8
        };
        split_kmers.push(split_kmer_info);
    }

    return split_kmers;
}
