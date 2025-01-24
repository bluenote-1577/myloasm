use crate::{constants, types::*};
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
                let pseudocov_ratio = pseudocount_cov(uni1.min_read_depth.unwrap(), uni2.min_read_depth.unwrap());
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
                - overhangs[0] as i64
                - overhangs[1] as i64;
            if end < 0 {
                end = 0;
            }
            range = (carryover, end as usize);
            hang_ind += 2;
        } else if i == node.read_indices_ori.len() - 1 {
            range = (
                carryover,
                reads[ind].base_length - right_cut.min(reads[ind].base_length),
            );
            hang_ind += 2;
        } else {
            let mut end = reads[ind].base_length as i64 + 1
                - internal_ol_len[hang_ind] as i64
                - overhangs[hang_ind] as i64
                - overhangs[hang_ind + 1] as i64;
            if end < 0 {
                end = 0;
            }
            range = (carryover, end as usize);
            hang_ind += 2;
        }
        if range.0 >= range.1 {
            carryover -= range.0 - range.1;
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
                base_seq.append(&reads[ind].dna_seq.revcomp()[range.0..range.1]);
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

pub fn get_base_info_mapchunks(
    _left_cut: usize,
    _right_cut: usize,
    _node: &UnitigNode,
    _reads: &[TwinRead],
    _dna_seq_info: bool,
) -> BaseInfo {
    return BaseInfo::default();
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

pub fn quantile_dist(v1: &[(f64, usize)], v2: &[(f64, usize)]) -> Option<f64> {
    let quants = [0.25, 0.5, 0.75];
    let quantiles_v1 = quants.iter().map(|x| median_weight(v1, *x)).collect::<Vec<_>>();
    let quantiles_v2 = quants.iter().map(|x| median_weight(v2, *x)).collect::<Vec<_>>();
    let mut minimum_distance = vec![];
    for i in 0..quants.len() {
        for j in 0..quants.len() {
            if quantiles_v1[i].is_none() || quantiles_v2[j].is_none() {
                continue;
            }
            //minimum_distance.push((quantiles_v2[j].unwrap() - quantiles_v1[i].unwrap()).abs());
            let qvi = quantiles_v1[i].unwrap();
            let qvj = quantiles_v2[j].unwrap();
            minimum_distance.push(pseudocount_cov(qvi, qvj).abs() - 1.);
        }
    }
    if minimum_distance.len() == 0 {
        return None;
    }
    Some(*minimum_distance.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap())
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
    let pcount = constants::PSEUDOCOUNT;
    let numerator = cov1.max(cov2) + pcount;
    let denominator = cov1.min(cov2) + pcount;
    return numerator / denominator;
}


#[cfg(test)]
mod tests {
    use constants::PSEUDOCOUNT;

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
    fn test_quantile_dist() {
        let v1 = vec![(1.0, 1), (2.0, 2), (3.0, 3)];
        let v2 = vec![(1.0, 1), (2.0, 2), (3.0, 3)];
        assert_eq!(quantile_dist(&v1, &v2), Some(0.0));

        let v1 = vec![(1.0, 1), (2.0, 1), (3.0, 1)];
        let v2 = vec![(2.0, 10), (3.0, 10), (4.0, 10)];
        assert_eq!(quantile_dist(&v1, &v2), Some(0.0));

        let v1 = vec![(10.0, 5), (10.0, 6), (20.0, 5)];
        let v2 = vec![(5.0, 1), (5.0, 3), (5.0, 5)];
        assert_eq!(quantile_dist(&v1, &v2), Some((10. + PSEUDOCOUNT) / (5.0 + PSEUDOCOUNT) - 1.));
    }
}
