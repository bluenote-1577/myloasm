/// Minimal-perfect-hash based minimizer index.
///
/// Construction uses NUM_PARTITIONS passes over the input data.  On each pass
/// only k-mers whose `to_u64() % NUM_PARTITIONS == partition` are collected
/// into a partial FxHashMap.  High-multiplicity k-mers are removed from that
/// partial map (same top-100_000 filter as `get_minimizer_index`), and then the
/// map is drained into the key/value vectors for MPHF construction.  This caps
/// peak memory to roughly 1/NUM_PARTITIONS of what a single-pass build would
/// require.

use boomphf::hashmap::BoomHashMap;
use fxhash::FxHashMap;

use crate::mapping::{Anchors, HitInfo};
use crate::types::{Anchor, AnchorBuilder, FlagKmer48, Kmer48, TwinRead};
use crate::utils;

const NUM_PARTITIONS: u64 = 3;

pub struct MphMinIndex {
    inner: BoomHashMap<Kmer48, Vec<HitInfo>>,
}

impl MphMinIndex {
    #[inline]
    pub fn get(&self, kmer: &Kmer48) -> Option<&Vec<HitInfo>> {
        self.inner.get(kmer)
    }
}

/// Collect one partition bucket from the given twin-read map into `partial`.
fn fill_partition(
    tr: impl Iterator<Item = (usize, Vec<(u32, FlagKmer48)>)>,
    partition: u64,
    partial: &mut FxHashMap<Kmer48, Vec<HitInfo>>,
) {
    for (id, minimizers) in tr {
        for (pos, mini) in minimizers {
            let kmer = mini.kmer();
            if kmer.to_u64() % NUM_PARTITIONS != partition {
                continue;
            }
            let hit = HitInfo {
                contig_id_strand: (mini.strand() as u32) << 31 | id as u32,
                position: pos,
            };
            partial.entry(kmer).or_insert_with(Vec::new).push(hit);
        }
    }
}

/// Apply the top-100_000 count filter to `partial` in-place, mirroring the
/// threshold logic in `get_minimizer_index`.
fn threshold_partial(partial: &mut FxHashMap<Kmer48, Vec<HitInfo>>) {
    if partial.len() <= 500_000 {
        return;
    }
    let mut counts: Vec<usize> = partial.values().map(|v| v.len()).collect();
    counts.sort_unstable_by(|a, b| b.cmp(a));
    let threshold = counts[counts.len() / 100_000];
    partial.retain(|_, v| v.len() < threshold);
}

/// Build a minimal-perfect-hash minimizer index.
///
/// Mirrors `get_minimizer_index` in mapping.rs but uses a `BoomHashMap`
/// instead of an `FxHashMap`, constructed with NUM_PARTITIONS passes to
/// limit peak memory.
pub fn get_minimizer_index_mph(
    tr_owned: Option<&FxHashMap<usize, TwinRead>>,
    tr_ref: Option<&FxHashMap<usize, &TwinRead>>,
) -> MphMinIndex {
    let mut all_keys: Vec<Kmer48> = Vec::new();
    let mut all_values: Vec<Vec<HitInfo>> = Vec::new();

    // Pre-sort IDs once so partitions are built deterministically.
    let sorted_ids: Vec<usize> = if let Some(m) = tr_owned {
        let mut v: Vec<usize> = m.keys().cloned().collect();
        v.sort_unstable();
        v
    } else if let Some(m) = tr_ref {
        let mut v: Vec<usize> = m.keys().cloned().collect();
        v.sort_unstable();
        v
    } else {
        panic!("No minimizer data provided");
    };

    for partition in 0..NUM_PARTITIONS {
        let mut partial: FxHashMap<Kmer48, Vec<HitInfo>> = FxHashMap::default();

        if let Some(twinreads) = tr_owned {
            fill_partition(
                sorted_ids.iter().map(|&id| (id, twinreads[&id].minimizers_vec_strand())),
                partition,
                &mut partial,
            );
        } else if let Some(twinreads) = tr_ref {
            fill_partition(
                sorted_ids.iter().map(|&id| (id, twinreads[&id].minimizers_vec_strand())),
                partition,
                &mut partial,
            );
        }

        // Filter high-multiplicity k-mers, then drain into the flat vecs.
        threshold_partial(&mut partial);
        all_keys.reserve(partial.len());
        all_values.reserve(partial.len());
        for (k, v) in partial.into_iter() {
            all_keys.push(k);
            all_values.push(v);
        }
        utils::log_memory_usage(false, &format!("MPH after partition {}", partition));
    }

    all_keys.shrink_to_fit();
    all_values.shrink_to_fit();

    log::debug!("MPH index: {} distinct k-mers", all_keys.len());

    if all_keys.is_empty() {
        // BoomHashMap panics on empty input; return a trivially empty wrapper.
        return MphMinIndex {
            inner: BoomHashMap::new(vec![], vec![]),
        };
    }

    MphMinIndex {
        inner: BoomHashMap::new_parallel(all_keys, all_values),
    }
}

/// Query an `MphMinIndex` for anchors between `seq1` and the indexed reference.
///
/// Mirrors `find_exact_matches_with_full_index` in mapping.rs exactly.
pub fn find_exact_matches_with_full_index_mph(
    seq1: &[(u32, FlagKmer48)],
    index: &MphMinIndex,
) -> FxHashMap<u32, Anchors> {
    let mut max_mult = 0;
    let mut matches: FxHashMap<u32, Vec<AnchorBuilder>> = FxHashMap::default();

    for (pos, flag_kmer) in seq1.iter() {
        let s1 = flag_kmer.strand() as u32;
        if let Some(indices) = index.get(&flag_kmer.kmer()) {
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for hit in indices {
                let s2 = hit.contig_id_strand >> 31;
                let contig = hit.contig_id_strand & 0x7FFF_FFFF;
                let rel_strand = s1 ^ s2;
                let anchor = AnchorBuilder {
                    pos1: (rel_strand << 31) | *pos,
                    pos2: hit.position,
                };
                matches.entry(contig).or_insert_with(Vec::new).push(anchor);
            }
        }
    }

    matches
        .into_iter()
        .map(|(k, v)| {
            let mut anchors = v
                .into_iter()
                .map(|a| Anchor {
                    i: None,
                    j: None,
                    pos1: a.pos1,
                    pos2: a.pos2,
                })
                .collect::<Vec<_>>();
            anchors.sort_by_key(|a| (a.pos1, a.pos2));
            (k, Anchors { anchors, max_mult })
        })
        .collect()
}
