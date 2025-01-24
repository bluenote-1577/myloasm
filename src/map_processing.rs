use crate::cli::Cli;
use crate::constants::ENDPOINT_MAPPING_FUZZ;
use crate::constants::MINIMIZER_END_NTH_COV;
use crate::constants::MIN_READ_LENGTH;
use crate::constants::SAMPLING_RATE_COV;
use crate::twin_graph;
use std::io::Write;
use std::io::BufWriter;
use std::path::Path;
use crate::types::*;
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_lapper::Lapper;
use std::sync::Mutex;

pub struct TwinReadMapping {
    pub tr_index: usize,
    pub mapping_info: MappingInfo,
    pub lapper_strain_max: Lapper<u32, bool>,
}

impl NodeMapping for TwinReadMapping {
    fn median_mapping_depth(&self) -> f64 {
        self.mapping_info.median_depth
    }
    fn mapping_boundaries(&self) -> &Lapper<u32, SmallTwinOl> {
        &self.mapping_info.mapping_boundaries
    }
    fn min_mapping_depth(&self) -> f64 {
        self.mapping_info.minimum_depth
    }
    fn set_mapping_info(&mut self, mapping_info: MappingInfo) {
        self.mapping_info = mapping_info;
    }
    fn mapping_info_present(&self) -> bool {
        self.mapping_info.present
    }
    fn reference_length(&self) -> usize {
        self.mapping_info.length
    }
    fn mapped_indices(&self) -> Vec<usize> {
        self.mapping_boundaries().iter().map(|x| x.val.query_id as usize).collect()
    }
}

pub fn median_and_min_depth_from_lapper<T>(lapper: &Lapper<u32, T>, sampling: usize, seq_start: usize, seq_length: usize) -> Option<(f64,f64)> 
where
    T: Eq + Clone + Send + Sync
{
    //TODO the block size should depend on read length
    let block_size = 20_000;
    let mut min_blocks = vec![];
    let mut median_blocks = vec![];
    let mut depths = vec![];
    if sampling > sampling{
        return None;
    }
    
    //Sample depth every 'sampling' bases, pad by block_size.
    let mut next_block = block_size;
    for pos in ((seq_start + sampling)..(seq_length - sampling)).step_by(sampling)
    {
        let depth = lapper.count(pos as u32, pos as u32 + 1);
        depths.push(depth);
        if pos > next_block && depths.len() > 0{
            depths.sort();
            next_block += block_size;
            let min_ind = 3.min(depths.len()-1);
            let min = depths[min_ind];
            let median = depths[depths.len() / 2];
            for _ in 0..block_size / sampling{
                min_blocks.push(min);
                median_blocks.push(median);
            }
            depths = vec![];
        }
    }

    //Leftover blocks
    if depths.len() > 0{
        depths.sort();
        let min_ind = 3.min(depths.len()-1);
        let min = depths[min_ind];
        let median = depths[depths.len() / 2];
        for _ in 0..depths.len(){
            min_blocks.push(min);
            median_blocks.push(median);
        }
    }

    if min_blocks.len() == 0{
        log::debug!("No blocks found for depth calculation, start {} end {}", seq_start, seq_length);
        return Some((0.,0.))
    }
    
    min_blocks.sort();
    median_blocks.sort();
    let median_over_min_blocks = min_blocks[min_blocks.len() / 2];
    let median_over_median_blocks = median_blocks[median_blocks.len() / 2];
    return Some((median_over_min_blocks as f64, median_over_median_blocks as f64));
}

//Inefficient implemtation. Prob don't need to optimize though
pub fn sliding_window_kmer_coverages(cov_vec: &Vec<(u32,u32)>, window_size: usize, k: usize) -> f64{
    let mut max_covs = vec![];
    let mut rolling_window = vec![0;window_size];
    let mut curr_pos = 0;
    for i in 0..cov_vec.len() {
        let pos = cov_vec[i].0;
        let count = cov_vec[i].1;
        if pos > k as u32 + curr_pos || curr_pos == 0{
            curr_pos = pos;
            rolling_window[i % window_size] = count;
        }
        if i >= window_size{
            let max = *rolling_window.iter().max().unwrap();
            max_covs.push(max);
        }
    }
    log::trace!("Max covs: {:?}", max_covs);
    log::trace!("Raw covs: {:?}", cov_vec);
    let returned_cov = (*max_covs.iter().min().unwrap()) as f64;
    log::trace!("Returned cov: {}", returned_cov);
    return returned_cov;
}

pub fn cov_mapping_breakpoints<T>(mapped: &T) -> Vec<Breakpoints>
where
    T: NodeMapping,
{
    if mapped.mapping_boundaries().intervals.is_empty() {
        return vec![];
    }
    let mut breakpoints = vec![];
    let depths = mapped.mapping_boundaries().depth().collect::<Vec<_>>();
    if depths.len() < 3 {
        return vec![];
    }
    // <ooooo|-------
    if depths[0].start > ENDPOINT_MAPPING_FUZZ{
        breakpoints.push(Breakpoints {
            pos1: 0,
            pos2 : depths[0].start as usize,
            cov: 0,
        });
    }
    //<xxxxx|--------
    //Cut bad left endpoints. 
    let depth_start_right = if mapped.reference_length() > (ENDPOINT_MAPPING_FUZZ as usize) + depths[0].stop as usize{
        mapped
            .mapping_boundaries()
            .count(depths[0].stop + ENDPOINT_MAPPING_FUZZ - 1, depths[0].stop + ENDPOINT_MAPPING_FUZZ)
    } else {
        0
    };
    if depths[0].stop > ENDPOINT_MAPPING_FUZZ && (depths[1].val > 3 || depth_start_right > 3) && depths[0].val == 1 {
        breakpoints.push(Breakpoints {
            pos1: depths[0].start as usize,
            pos2: depths[0].stop as usize,
            cov: depths[0].val as usize,
        });
    }
    // -----|xxxx|----
    for i in 1..depths.len() - 1 {
        let interval = &depths[i];
        let last_cov = depths[i - 1].val;
        let next_cov = depths[i + 1].val;
        let start = interval.start;
        let stop = interval.stop as usize;
        let cov = interval.val as usize;
        let cond1 = last_cov > 3 && next_cov > 3 && cov == 1;
        let cond2;
        if start > 200 && stop + 200 < mapped.reference_length() {
            let left_count = mapped.mapping_boundaries().count(start - 200, start - 198);
            let right_count = mapped.mapping_boundaries().count(start as u32 + 198, start as u32 + 200);
            cond2 = left_count > 3 && right_count > 3 && cov == 1;
        } else {
            cond2 = false;
        }
        let cond3 = (last_cov > (cov as u32 * 5) || next_cov > (cov as u32 * 5)) && cov < 3;
        let cond4 = (last_cov > (cov as u32 * 25) || next_cov > (cov as u32 * 25)) && cov == 3;
        if start > 200
            && start as usize + 200 < mapped.reference_length()
            && (cond1 || cond2 || cond3 || cond4)
        {
            breakpoints.push(Breakpoints {
                pos1: start as usize,
                pos2: stop as usize,
                cov: cov,
            });
        }
    }

    // --------|xxxxx>
    if depths[depths.len() - 1].start > ENDPOINT_MAPPING_FUZZ {
        let depth_stop_left = mapped.mapping_boundaries().count(
            depths[depths.len() - 1].start - ENDPOINT_MAPPING_FUZZ,
            depths[depths.len() - 1].start - ENDPOINT_MAPPING_FUZZ - 1,
        );
        if (depth_stop_left > 3 || depths[depths.len() - 2].val > 3)
            && depths[depths.len() - 1].val == 1
        {
            breakpoints.push(Breakpoints {
                pos1: depths[depths.len() - 1].start as usize,
                pos2: depths[depths.len() - 1].stop as usize,
                cov: 0,
            });
        }
    }

    // -----|ooooo>
    if depths[depths.len() - 1].stop as usize + (ENDPOINT_MAPPING_FUZZ as usize) < mapped.reference_length() {
        breakpoints.push(Breakpoints {
            pos1: depths[depths.len() - 1].stop as usize,
            pos2: mapped.reference_length(),
            cov: depths[depths.len() - 1].val as usize,
        });
    }

    return breakpoints;
}

fn split_read_and_populate_depth(mut twin_read: TwinRead, mapping_info: &TwinReadMapping, mut break_points: Vec<Breakpoints>, args: &Cli) -> Vec<TwinRead>{

    //Return the read, populate the depth from mapping_info
    //TODO this doesn't do anything -- we add a BP and then populate it later down...
    if break_points.len() == 0{
        twin_read.median_depth = Some(mapping_info.median_mapping_depth());
        twin_read.min_depth = Some(mapping_info.min_mapping_depth());
    }

    let mut new_reads = vec![];
    break_points.push(Breakpoints{pos1: twin_read.base_length, pos2: twin_read.base_length, cov: 0});

    let mut last_break = 0;
    let k = twin_read.k as usize;
    
    for (i,break_point) in break_points.iter().enumerate(){
        let bp_start = break_point.pos1;
        let bp_end = break_point.pos2;
        if bp_start - last_break > MIN_READ_LENGTH{
            let mut new_read = TwinRead::default();
            //Repopulate minimizers and snpmers
            new_read.minimizers = twin_read.minimizers.iter().filter(|x| x.0 >= last_break && x.0 + k - 1 < bp_start).copied().map(|x| (x.0 - last_break, x.1)).collect();
            new_read.snpmers = twin_read.snpmers.iter().filter(|x| x.0 >= last_break && x.0 + k - 1 < bp_start).copied().map(|x| (x.0 - last_break, x.1)).collect();
            new_read.id = format!("{}+split{}", &twin_read.id, i);
            log::trace!("Split read {} at {}-{}", &new_read.id, last_break, bp_start);
            new_read.k = twin_read.k;
            new_read.dna_seq = twin_read.dna_seq[last_break..bp_start].to_owned();
            new_read.base_length = new_read.dna_seq.len();
            if break_points.len() > 1 {
                new_read.split_chimera = true;
            }
            else{
                populate_depth_from_map_info(&mut new_read, mapping_info, last_break, bp_start);
            }
            new_reads.push(new_read);
        }
        last_break = bp_end;
    }
    return new_reads;
}

#[inline]
pub fn populate_depth_from_map_info(twin_read: &mut TwinRead, mapping_info: &TwinReadMapping, start: usize, end: usize){
    let (first_mini, last_mini) = first_last_mini_in_range(start, end, twin_read.k as usize, MINIMIZER_END_NTH_COV, &twin_read.minimizers);
    let (min_depth, median_depth) = median_and_min_depth_from_lapper(&mapping_info.lapper_strain_max, SAMPLING_RATE_COV, first_mini, last_mini).unwrap();
    twin_read.min_depth = Some(min_depth);
    twin_read.median_depth = Some(median_depth);
}

pub fn first_last_mini_in_range(start: usize, end: usize, k: usize, nth: usize, minis: &[(usize,u64)]) -> (usize, usize){
    let mut first_mini = start;
    let mut count_first = 0;
    let mut last_mini = end;
    let mut count_last = 0;

    for (mini_pos, _) in minis.iter(){
        if *mini_pos >= start{
            count_first += 1;
            first_mini = *mini_pos;
        }
        if count_first == nth{
            break;
        }
    }

    for (mini_pos, _) in minis.iter().rev(){
        if mini_pos + k - 1 < end{
            last_mini = *mini_pos;
            count_last += 1;
        }
        if count_last == nth{
            break;
        }
    }

    if last_mini <= first_mini{
        return (start, end);
    }

    return (first_mini, last_mini);
}

pub fn split_outer_reads(twin_reads: Vec<TwinRead>, tr_map_info: Vec<TwinReadMapping>, args: &Cli)
-> (Vec<TwinRead>, Vec<usize>){
    let tr_map_info_dict = tr_map_info.iter().map(|x| (x.tr_index, x)).collect::<FxHashMap<usize, &TwinReadMapping>>();
    let new_twin_reads_bools = Mutex::new(vec![]);
    let cov_file = Path::new(args.output_dir.as_str()).join("read_coverages.txt");
    let writer = Mutex::new(BufWriter::new(std::fs::File::create(cov_file).unwrap()));

    twin_reads.into_par_iter().enumerate().for_each(|(i, twin_read)| {
        if tr_map_info_dict.contains_key(&i){
            let map_info = tr_map_info_dict.get(&i).unwrap();
            let breakpoints = cov_mapping_breakpoints(*map_info);

            if log::log_enabled!(log::Level::Trace) {
                let depths = map_info.mapping_boundaries().depth().collect::<Vec<_>>();
                let writer = &mut writer.lock().unwrap();
                for depth in depths{
                    let mut string = format!("{} {}-{} COV:{}, BREAKPOINTS:", twin_read.id, depth.start, depth.stop, depth.val);
                    for breakpoint in breakpoints.iter(){
                        string.push_str(format!("--{} to {}--", breakpoint.pos1, breakpoint.pos2).as_str());
                    }
                    writeln!(writer, "{}", &string).unwrap();
                }
            }
            else{
                let mut read_id_and_breakpoint_string = format!("{} BREAKPOINTS:", twin_read.id);
                for breakpoint in breakpoints.iter(){
                    read_id_and_breakpoint_string.push_str(format!("{}-{},", breakpoint.pos1, breakpoint.pos2).as_str());
                }
                let writer = &mut writer.lock().unwrap();
                writeln!(writer, "{}", &read_id_and_breakpoint_string).unwrap();
            }

            let splitted_reads = split_read_and_populate_depth(twin_read, map_info, breakpoints, args);
            for new_read in splitted_reads{
                new_twin_reads_bools.lock().unwrap().push((new_read, true));
            }
        }
        else{
            new_twin_reads_bools.lock().unwrap().push((twin_read, false));
        }
    });

    let mut ntr_bools = new_twin_reads_bools.into_inner().unwrap();
    ntr_bools.sort_by(|a,b| a.0.id.cmp(&b.0.id));
    let new_outer_indices = ntr_bools.iter().enumerate().filter(|x| x.1.1).map(|x| x.0).collect::<Vec<usize>>();
    let new_twin_reads = ntr_bools.into_iter().map(|x| x.0).collect::<Vec<TwinRead>>();
    return (new_twin_reads, new_outer_indices);
}

pub fn check_maximal_overlap(start1: usize, end1: usize, start2: usize, end2: usize, len1: usize, len2: usize, reverse: bool) -> bool {
    let edge_fuzz = ENDPOINT_MAPPING_FUZZ as usize;

    //Can not extend to the left (cond1) and cannot extend to the right (cond2)
    //  ------->             OR          --------->
    // ---------->                            --------->
    if !reverse{
        if (start1 < edge_fuzz || start2 < edge_fuzz ) && (len1 < edge_fuzz + end1 || len2 < edge_fuzz + end2){
            return true;
        }
    }

    //  ------->             OR          --------->            OR          <----------
    //      <------                    <--------------                            -------->
    else{
        let max1right = len1 < edge_fuzz + end1;
        let max2right = len2 < edge_fuzz + end2;
        let max1left = start1 < edge_fuzz;
        let max2left = start2 < edge_fuzz;

        let ol_plus_minus = max1right && max2right;
        let ol_minus_plus = max1left && max2left;
        let contained1 = max1left && max1right;
        let contained2 = max2left && max2right;

        if ol_plus_minus || ol_minus_plus || contained1 || contained2{
            return true;
        }
    }

    return false;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_maximal_overlap(){
        let len1 = 3000;
        let len2 = 3000;

        //  ------>
        // ---------->
        let start1 = 0;
        let end1 = 3000;
        let start2 = 0;
        let end2 = 3000;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), true);

        //  ------>
        // <----------
        let reverse = true;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), true);

        // ------------>
        //    ----->
        let len1 = 3000;
        let len2 = 2000;
        let start1 = 500;
        let end1 = 2500;
        let start2 = 0;
        let end2 = 2000;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), true);

        let len1 = 3000;
        let len2 = 3000;

        // ----->
        //    ------->
        let start1 = 1500;
        let end1 = 2900;
        let start2 = 100;
        let end2 = 1500;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), true);

        //     ------>
        // ------->
        let start1 = 100;
        let end1 = 1500;
        let start2 = 1400;
        let end2 = 2900;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), true);

        //   ---xxxx>
        // -----x>
        let start1 = 100;
        let end1 = 500;
        let start2 = 2000;
        let end2 = 2500;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), false);

        // <-----xxxx 
        //  -----xxxx->
        let start1 = 1500;
        let end1 = 3000;
        let start2 = 0;
        let end2 = 1500;
        let reverse = true;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), false);

        // <------
        //    -------->
        let start1 = 0;
        let end1 = 1500;
        let start2 = 0;
        let end2 = 1500;
        let reverse = true;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), true);

        // <--|||xx
        //    |||xx--->
        let start1 = 500;
        let end1 = 1500;
        let start2 = 0;
        let end2 = 1000;
        let reverse = true;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse), false);
    }
}
