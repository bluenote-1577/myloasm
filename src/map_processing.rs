use crate::cli::Cli;
use crate::constants::READ_BLOCK_SIZE_FOR_COVERAGE;
use crate::constants::ENDPOINT_MAPPING_FUZZ;
use crate::constants::IDENTITY_THRESHOLDS;
use crate::constants::ID_THRESHOLD_ITERS;
use crate::constants::MINIMIZER_END_NTH_COV;
use crate::constants::MIN_COV_READ;
use crate::constants::MIN_READ_LENGTH;
use crate::constants::SAMPLING_RATE_COV;
use std::io::Write;
use flate2::write::GzEncoder;
use rust_lapper::Interval;
use std::io::BufWriter;
use std::path::Path;
use crate::types::*;
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_lapper::Lapper;
use std::sync::Mutex;
use std::collections::BTreeMap;
use fxhash::FxHashSet;
use flate2::Compression;
use std::path::PathBuf;


pub struct TwinReadMapping {
    pub tr_index: usize,
    pub mapping_info: MappingInfo,
    pub all_intervals: Vec<BareInterval>,
}


pub fn median_and_min_depth_from_lapper_new(
    lapper: &Lapper<u32, BareMappingOverlap>, 
    sampling: usize, 
    seq_start: usize, 
    seq_length: usize, 
    snpmer_identity_cutoff: f64
) -> Option<(f64, f64)> {
    let block_size = READ_BLOCK_SIZE_FOR_COVERAGE;
    
    // Filter intervals based on snpmer_identity and maximal_overlap
    let intervals_cutoff = lapper
        .iter()
        .filter(|x| (x.val.snpmer_identity as f64 >= snpmer_identity_cutoff))
        .map(|x| x.clone())
        .collect::<Vec<_>>();
    let lapper_cutoff = Lapper::new(intervals_cutoff);
    
    // If no qualifying intervals, return early
    if lapper_cutoff.intervals.is_empty() {
        log::trace!("No qualifying intervals found for depth calculation, start {} end {}", seq_start, seq_length);
        return Some((0., 0.));
    }
    
    // Get depths at regular intervals using our new function
    let depths_at_points = depths_at_points(
        &lapper_cutoff,
        (seq_start + sampling) as u32,
        (seq_length - sampling) as u32,
        sampling as u32
    );
    
    if depths_at_points.is_empty() {
        log::trace!("No points found for depth calculation, start {} end {}", seq_start, seq_length);
        return Some((0., 0.));
    }
    
    // Process depths in blocks
    let mut min_blocks = Vec::new();
    let mut median_blocks = Vec::new();
    let mut current_block = Vec::new();
    let mut next_block_boundary = block_size + seq_start;
    
    for (pos, depth) in depths_at_points {
        current_block.push(depth as u32);
        
        if pos as usize >= next_block_boundary && !current_block.is_empty() {
            // Process the completed block
            current_block.sort_unstable();
            let min = current_block.first().unwrap();
            let median = current_block[current_block.len() / 2];
            
            // Add block_size/sampling copies of the min and median
            let copies = block_size / sampling;
            min_blocks.extend(std::iter::repeat(min).take(copies));
            median_blocks.extend(std::iter::repeat(median).take(copies));
            
            // Reset for next block
            current_block.clear();
            next_block_boundary += block_size;
        }
    }
    
    // Process any remaining points in the last block
    if !current_block.is_empty() {
        current_block.sort_unstable();
        let min = current_block.first().unwrap();
        let median = current_block[current_block.len() / 2];
        
        let remaining = current_block.len();
        min_blocks.extend(std::iter::repeat(min).take(remaining));
        median_blocks.extend(std::iter::repeat(median).take(remaining));
    }
    
    if min_blocks.is_empty() {
        log::trace!("No blocks processed for depth calculation, start {} end {}", seq_start, seq_length);
        return Some((0., 0.));
    }
    
    // Calculate final medians
    min_blocks.sort_unstable();
    median_blocks.sort_unstable();
    let median_over_min_blocks: u32 = min_blocks[min_blocks.len() / 2];
    let median_over_median_blocks : u32 = median_blocks[median_blocks.len() / 2];
    
    Some((median_over_min_blocks as f64, median_over_median_blocks as f64))
}

pub fn median_and_min_depth_from_lapper(lapper: &Lapper<u32, BareMappingOverlap>, sampling: usize, seq_start: usize, seq_length: usize, snpmer_identity_cutoff: f64) -> Option<(f64,f64)> 
{
    //TODO the block size should depend on read length
    let block_size = READ_BLOCK_SIZE_FOR_COVERAGE;
    let mut min_blocks = vec![];
    let mut median_blocks = vec![];
    let mut depths = vec![];
    let intervals_cutoff = lapper
        .iter()
        .filter(|x| (x.val.snpmer_identity as f64 >= snpmer_identity_cutoff))
        .map (|x| x.clone())
        .collect::<Vec<_>>();
    let lapper_cutoff = Lapper::new(intervals_cutoff);
    
    //Sample depth every 'sampling' bases, pad by block_size.
    let mut next_block = block_size;
    for pos in ((seq_start + sampling)..(seq_length - sampling)).step_by(sampling)
    {
        let depth = lapper_cutoff.count(pos as u32, pos as u32 + 1);
        depths.push(depth);
        if pos > next_block && depths.len() > 0{
            depths.sort();
            next_block += block_size;
            //let min_ind = 3.min(depths.len()-1);
            //let min = depths[min_ind];
            let min = depths[0];
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
        //let min_ind = 3.min(depths.len()-1);
        //let min = depths[min_ind];
        let min = depths[0];
        let median = depths[depths.len() / 2];
        for _ in 0..depths.len(){
            min_blocks.push(min);
            median_blocks.push(median);
        }
    }

    if min_blocks.len() == 0{
        log::trace!("No blocks found for depth calculation, start {} end {}", seq_start, seq_length);
        return Some((0.,0.))
    }
    
    min_blocks.sort();
    median_blocks.sort();
    let median_over_min_blocks = min_blocks[min_blocks.len() / 2];
    let median_over_median_blocks = median_blocks[median_blocks.len() / 2];
    return Some((median_over_min_blocks as f64, median_over_median_blocks as f64));
}

pub fn cov_mapping_breakpoints(intervals: &Vec<BareInterval>, reference_length: u32) -> Vec<Breakpoints>
{
    if intervals.is_empty() {
        return vec![];
    }

    let mut breakpoints = vec![];
    let sampling = 5;
    let coverages_at_sampling = depths_at_points_interval(
        intervals,
        0,
        reference_length,
        sampling,
    );
    let intervals = intervals.iter().map(|x| Interval{start: x.start, stop: x.stop, val: false}).collect::<Vec<Interval<u32,bool>>>();
    let lapper = Lapper::new(intervals);
    let depths = lapper.depth().collect::<Vec<_>>();
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
    let depth_start_right = if reference_length as usize > (ENDPOINT_MAPPING_FUZZ as usize) + depths[0].stop as usize{
        lapper 
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
        if start > 200 && stop + 200 < reference_length as usize {
            // let left_count = mapped.mapping_boundaries().count(start - 200, start - 198);
            // let right_count = mapped.mapping_boundaries().count(start as u32 + 198, start as u32 + 200);

            let left_count = coverages_at_sampling[(start as usize - 200) / sampling as usize].1;
            let right_count = coverages_at_sampling[(start as usize + 200) / sampling as usize].1;
            cond2 = left_count > 3 && right_count > 3 && cov == 1;
        } else {
            cond2 = false;
        }
        let cond3 = (last_cov > (cov as u32 * 5) || next_cov > (cov as u32 * 5)) && cov < 3;
        let cond4 = (last_cov > (cov as u32 * 25) || next_cov > (cov as u32 * 25)) && cov == 3;
        if start > 200
            && start as usize + 200 < reference_length as usize
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
    if depths[depths.len() - 1].start < reference_length - ENDPOINT_MAPPING_FUZZ {
        let depth_stop_left = lapper.count(
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
    if depths[depths.len() - 1].stop as usize + (ENDPOINT_MAPPING_FUZZ as usize) < reference_length as usize {
        breakpoints.push(Breakpoints {
            pos1: depths[depths.len() - 1].stop as usize,
            pos2: reference_length as usize ,
            cov: depths[depths.len() - 1].val as usize,
        });
    }

    return breakpoints;
}

fn split_read_and_populate_depth(mut twin_read: TwinRead, mapping_info: &TwinReadMapping, mut break_points: Vec<Breakpoints>, _args: &Cli) -> Vec<TwinRead>{

    let mut new_reads = vec![];
    break_points.push(Breakpoints{pos1: twin_read.base_length, pos2: twin_read.base_length, cov: 0});

    let mut last_break = 0;
    let k = twin_read.k as usize;
    
    for (i,break_point) in break_points.iter().enumerate(){

        let bp_start = break_point.pos1;
        let bp_end = break_point.pos2;

        if bp_start - last_break > MIN_READ_LENGTH{
            let mut new_read = TwinRead::default();
            //Repopulate minimizers and snpmers. This is kind of bad; should have a method that does this...
            new_read.shift_and_retain(&twin_read, last_break, bp_start, k);
            new_read.id = format!("{}+split{}", &twin_read.id, i);
            new_read.split_start = last_break as u32;
            log::trace!("Split read {} at {}-{}", &new_read.id, last_break, bp_start);
            new_read.k = twin_read.k;
            new_read.dna_seq = twin_read.dna_seq[last_break..bp_start].to_owned();
            new_read.est_id = twin_read.est_id;

            if let Some(qual_seq) = twin_read.qual_seq.as_ref(){
                new_read.qual_seq = Some(qual_seq[last_break..bp_start].to_owned());
            }

            new_read.base_length = new_read.dna_seq.len();
            new_read.base_id = twin_read.base_id.clone();

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
    twin_read.clear();
    return new_reads;
}

pub fn populate_depth_from_map_info(twin_read: &mut TwinRead, mapping_info: &TwinReadMapping, start: usize, end: usize){
    let (first_mini, last_mini) = first_last_mini_in_range(start, end, twin_read.k as usize, MINIMIZER_END_NTH_COV, &twin_read.minimizer_positions);
    let mut min_depths = [0.; ID_THRESHOLD_ITERS];
    let mut median_depth = 0.;

    for (i,id) in IDENTITY_THRESHOLDS.iter().enumerate() {
        let (min_depth, median_depth_t) = median_and_min_depth_from_lapper_new(mapping_info.mapping_info.max_mapping_boundaries.as_ref().unwrap(), SAMPLING_RATE_COV, first_mini, last_mini, *id).unwrap();
        median_depth = median_depth_t;
        min_depths[i] = min_depth;
    }

    // Change snpmer id threshold based on coverage
    if twin_read.snpmer_id_threshold.is_none(){

        let step = 0.05 / 100.;

        let sufficient_gap = min_depths[0] > 3.0 * min_depths[ID_THRESHOLD_ITERS - 1];
        let sufficient_depth = min_depths[0] >= MIN_COV_READ as f64;
        let strain_specific_not_enough = min_depths[ID_THRESHOLD_ITERS - 1] < MIN_COV_READ as f64; // Always
        let pass_cond_1 = sufficient_gap;
        let pass_cond_2 = strain_specific_not_enough;

        let mut prev_id = IDENTITY_THRESHOLDS[ID_THRESHOLD_ITERS - 1];

        //TODO
        if (pass_cond_1 && pass_cond_2) && sufficient_depth {
            //0.10% increments; 100 -> 99.90 -> 99.8 ...
            let mut min_depth_prev = min_depths[ID_THRESHOLD_ITERS - 1];
            let mut try_id = IDENTITY_THRESHOLDS[ID_THRESHOLD_ITERS - 1] - step;
            loop{
                let (min_depth_try, _) = median_and_min_depth_from_lapper_new(mapping_info.mapping_info.max_mapping_boundaries.as_ref().unwrap(), SAMPLING_RATE_COV, first_mini, last_mini, try_id).unwrap();

                //if min_depth_try >= MIN_COV_READ as f64  && min_depth_try < min_depth_prev * 1.50 {
                if min_depth_try <= min_depth_prev * 1.50 && min_depth_try >= MIN_COV_READ as f64 {
                    log::trace!("Read {} min_depth_identity:{} cov:{}", twin_read.id, prev_id * 100., min_depth_try);
                    twin_read.snpmer_id_threshold = Some(prev_id * 100.);
                    break;
                }

                prev_id = try_id;
                min_depth_prev = min_depth_try;
                try_id -= step;
            }
        }
        else{
            twin_read.snpmer_id_threshold = Some(IDENTITY_THRESHOLDS[ID_THRESHOLD_ITERS - 1] * 100.);
        }
    }

    twin_read.min_depth_multi = Some(min_depths);
    twin_read.median_depth = Some(median_depth);

    //Percentage, not fractinoal
    if twin_read.snpmer_id_threshold.is_some(){
        assert!(twin_read.snpmer_id_threshold.unwrap() > 1.1);
    }
}

pub fn first_last_mini_in_range(start: usize, end: usize, k: usize, nth: usize, minis: &[u32]) -> (usize, usize){
    let mut first_mini = start;
    let mut count_first = 0;
    let mut last_mini = end;
    let mut count_last = 0;

    for mini_pos in minis.iter(){
        if *mini_pos as usize >= start{
            count_first += 1;
            first_mini = *mini_pos as usize;
        }
        if count_first == nth{
            break;
        }
    }

    for mini_pos in minis.iter().rev(){
        if *mini_pos as usize + k - 1 < end{
            last_mini = *mini_pos as usize;
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

pub fn split_outer_reads(twin_reads: Vec<TwinRead>, tr_map_info: Vec<TwinReadMapping>, temp_dir: &PathBuf, args: &Cli)
-> (Vec<TwinRead>, Vec<usize>){
    let tr_map_info_dict = tr_map_info.iter().map(|x| (x.tr_index, x)).collect::<FxHashMap<usize, &TwinReadMapping>>();
    let new_twin_reads_bools = Mutex::new(vec![]);
    let cov_file = Path::new(temp_dir).join("read_coverages.txt.gz");
    let buf = BufWriter::new(std::fs::File::create(cov_file).unwrap());
    let writer = Mutex::new(GzEncoder::new(buf, Compression::default()));

    twin_reads.into_par_iter().enumerate().for_each(|(i, twin_read)| {
        if tr_map_info_dict.contains_key(&i){
            let map_info = tr_map_info_dict.get(&i).unwrap();
            let breakpoints = cov_mapping_breakpoints(&map_info.all_intervals, map_info.mapping_info.length as u32);

            //Broke right now 
            if log::log_enabled!(log::Level::Trace) && false {
                let depths = map_info.mapping_info.max_mapping_boundaries.as_ref().unwrap().depth().collect::<Vec<_>>();
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
                if breakpoints.len() > 0{
                    let mut read_id_and_breakpoint_string = format!("{} BREAKPOINTS:", twin_read.id);
                    for breakpoint in breakpoints.iter(){
                        read_id_and_breakpoint_string.push_str(format!("{}-{},", breakpoint.pos1, breakpoint.pos2).as_str());
                    }
                    let writer = &mut writer.lock().unwrap();
                    writeln!(writer, "{}", &read_id_and_breakpoint_string).unwrap();
                }
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

pub fn check_maximal_overlap(start1: usize, end1: usize, start2: usize, end2: usize, len1: usize, len2: usize, reverse: bool, endpoint_fuzz: usize) -> bool {
    let edge_fuzz = endpoint_fuzz;

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

pub fn depths_at_points_interval(intervals: &Vec<BareInterval>, start: u32, end: u32, step: u32) -> Vec<(u32, usize)> {
    if intervals.is_empty() || start > end || step == 0 {
        return Vec::new();
    }

    let mut result = Vec::new();
    let points = (end - start) / step + 1;
    result.reserve(points as usize);
    
    let mut pos = start;
    let mut active_intervals = FxHashSet::default();
    let mut next_interval_idx = 0;
    let mut stops_of_active_intervals: BTreeMap<u32, Vec<usize>> = BTreeMap::new();

    while pos <= end {
        let mut toremove_intervals = vec![];
        let mut toremove_stops = vec![];

        // Remove intervals that end before or at current position
        for (stop,idsxs) in stops_of_active_intervals.iter(){
            if *stop <= pos{
                toremove_stops.push(*stop);
                toremove_intervals.extend(idsxs.iter().cloned());
            }
            else{
                break;
            }
        }

        for idx in toremove_intervals.iter() {
            active_intervals.remove(idx);
        }

        for stop in toremove_stops.iter() {
            stops_of_active_intervals.remove(stop);
        }

        // Add new intervals that start before or at current position
        while next_interval_idx < intervals.len() 
            && intervals[next_interval_idx].start <= pos 
        {
            if intervals[next_interval_idx].stop > pos {
                active_intervals.insert(next_interval_idx);
                stops_of_active_intervals.entry(intervals[next_interval_idx].stop).or_insert(vec![]).push(next_interval_idx);
            }
            next_interval_idx += 1;
        }

        // Record depth at current position
        result.push((pos, active_intervals.len()));

        // Move to next position
        pos = pos + step;
    }

    result
}

pub fn depths_at_points<T>(lapper: &Lapper<u32, T>, start: u32, end: u32, step: u32) -> Vec<(u32, usize)> 
where
    T: Eq + Clone + Send + Sync
{
    if lapper.intervals.is_empty() || start > end || step == 0 {
        return Vec::new();
    }

    let mut result = Vec::new();
    let points = (end - start) / step + 1;
    result.reserve(points as usize);
    
    let mut pos = start;
    let mut active_intervals = FxHashSet::default();
    let mut next_interval_idx = 0;
    let mut stops_of_active_intervals: BTreeMap<u32, Vec<usize>> = BTreeMap::new();

    while pos <= end {
        let mut toremove_intervals = vec![];
        let mut toremove_stops = vec![];

        // Remove intervals that end before or at current position
        for (stop,idsxs) in stops_of_active_intervals.iter(){
            if *stop <= pos{
                toremove_stops.push(*stop);
                toremove_intervals.extend(idsxs.iter().cloned());
            }
            else{
                break;
            }
        }

        for idx in toremove_intervals.iter() {
            active_intervals.remove(idx);
        }

        for stop in toremove_stops.iter() {
            stops_of_active_intervals.remove(stop);
        }

        // Add new intervals that start before or at current position
        while next_interval_idx < lapper.intervals.len() 
            && lapper.intervals[next_interval_idx].start <= pos 
        {
            if lapper.intervals[next_interval_idx].stop > pos {
                active_intervals.insert(next_interval_idx);
                stops_of_active_intervals.entry(lapper.intervals[next_interval_idx].stop).or_insert(vec![]).push(next_interval_idx);
            }
            next_interval_idx += 1;
        }

        // Record depth at current position
        result.push((pos, active_intervals.len()));

        // Move to next position
        pos = pos + step;
    }

    result
}


#[cfg(test)]
mod tests {

    use super::*;
    use approx::assert_relative_eq;


    fn create_small_twinol(identity: f32, _maximal: bool) -> BareMappingOverlap{
        BareMappingOverlap {
            snpmer_identity: identity,
            ..Default::default()
        }
    }

    #[test]
    fn test_maximal_overlap(){

        let fuzz = 200; 

        let len1 = 3000;
        let len2 = 3000;

        //  ------>
        // ---------->
        let start1 = 0;
        let end1 = 3000;
        let start2 = 0;
        let end2 = 3000;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), true);

        //  ------>
        // <----------
        let reverse = true;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), true);

        // ------------>
        //    ----->
        let len1 = 3000;
        let len2 = 2000;
        let start1 = 500;
        let end1 = 2500;
        let start2 = 0;
        let end2 = 2000;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), true);

        let len1 = 3000;
        let len2 = 3000;

        // ----->
        //    ------->
        let start1 = 1500;
        let end1 = 2900;
        let start2 = 100;
        let end2 = 1500;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), true);

        //     ------>
        // ------->
        let start1 = 100;
        let end1 = 1500;
        let start2 = 1400;
        let end2 = 2900;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), true);

        //   ---xxxx>
        // -----x>
        let start1 = 100;
        let end1 = 500;
        let start2 = 2000;
        let end2 = 2500;
        let reverse = false;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), false);

        // <-----xxxx 
        //  -----xxxx->
        let start1 = 1500;
        let end1 = 3000;
        let start2 = 0;
        let end2 = 1500;
        let reverse = true;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), false);

        // <------
        //    -------->
        let start1 = 0;
        let end1 = 1500;
        let start2 = 0;
        let end2 = 1500;
        let reverse = true;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), true);

        // <--|||xx
        //    |||xx--->
        let start1 = 500;
        let end1 = 1500;
        let start2 = 0;
        let end2 = 1000;
        let reverse = true;
        assert_eq!(check_maximal_overlap(start1, end1, start2, end2, len1, len2, reverse, fuzz), false);
    }

    // fn _test_get_min_depth_various_thresholds(){
    //     let mut intervals = vec![];
    //     for _ in 0..100{
    //         intervals.push((0 as u32, 0 as u32 + 100, OrderedFloat(0.99)));
    //     }

    //     for _ in 0..100{
    //         intervals.push((0 as u32, 0 as u32 + 100, OrderedFloat(1.00)));
    //     }

    //     for _ in 0..100{
    //         intervals.push((0 as u32, 0 as u32 + 100, OrderedFloat(0.98)));
    //     }
    //     let intervals = intervals.into_iter().map(|x| Interval{start: x.0, stop: x.1, val: x.2}).collect::<Vec<_>>();
    //     let lapper = Lapper::new(intervals);

    //     let min_depths = _get_min_depth_various_thresholds(&lapper, 10, 50, &[0.97, 0.98, 0.99, 1.0]);
    //     dbg!(&min_depths);
    //     assert!(min_depths[0] == 300.);
    //     assert!(min_depths[1] == 300.);
    //     assert!(min_depths[2] == 200.);
    //     assert!(min_depths[3] == 100.);

    //     let min_depths = _get_min_depth_various_thresholds(&lapper, 101, 120, &[0.97, 0.98, 0.99, 1.0]);
    //     dbg!(&min_depths);
    //     assert!(min_depths[0] == 0.);
    //     assert!(min_depths[1] == 0.);
    //     assert!(min_depths[2] == 0.);
    //     assert!(min_depths[3] == 0.);
    // }

    //#[test]
    // fn test_get_min_depth_various_thresholds_staggered(){
    //     let mut intervals = vec![];
    //     for i in 0..100{
    //         intervals.push((i as u32, i as u32 + 100, OrderedFloat(0.995)));
    //     }

    //     for j in 0..100{
    //         intervals.push((j as u32 + 500 , j as u32 + 500 + 100, OrderedFloat(1.00)));
    //     }

    //     for k in 0..100{
    //         intervals.push((k as u32, k as u32 + 100, OrderedFloat(0.97)));
    //     }
    //     let intervals = intervals.into_iter().map(|x| Interval{start: x.0, stop: x.1, val: x.2}).collect::<Vec<_>>();
    //     let lapper = Lapper::new(intervals);

    //     let min_depths = _get_min_depth_various_thresholds(&lapper, 20, 30, &[0.99]);
    //     dbg!(&min_depths);
    //     assert!(min_depths[0] == 21.);

    //     let min_depths = _get_min_depth_various_thresholds(&lapper, 0, 1, &[0.99]);
    //     dbg!(&min_depths);
    //     assert!(min_depths[0] == 1.);

    //     //Closed/open [) intervals
    //     let min_depths = _get_min_depth_various_thresholds(&lapper, 100, 100, &[0.99]);
    //     dbg!(&min_depths);
    //     assert!(min_depths[0] == 99.);

    //     let min_depths = _get_min_depth_various_thresholds(&lapper, 50, 5000, &[0.99]);
    //     dbg!(&min_depths);
    //     assert!(min_depths[0] == 0.);

    //     let min_depths = _get_min_depth_various_thresholds(&lapper, 50, 120, &[0.99]);
    //     dbg!(&min_depths);
    //     assert!(min_depths[0] == 51.);

    //     let min_depths = _get_min_depth_various_thresholds(&lapper, 50, 120, &[0.95]);
    //     dbg!(&min_depths);
    //     assert!(min_depths[0] == 102.);

    // }

    #[test]
    fn test_depths_at_points() {
        // Test case with overlapping intervals
        let data = vec![
            Interval { start: 0, stop: 10, val: SmallTwinOl::default() },
            Interval { start: 5, stop: 15, val: SmallTwinOl::default() },
            Interval { start: 10, stop: 20, val: SmallTwinOl::default() },
        ];
        let lapper = Lapper::new(data);
        
        let depths = depths_at_points(&lapper, 0, 20, 5);
        assert_eq!(depths, vec![
            (0, 1),   // Only first interval
            (5, 2),   // First and second intervals
            (10, 2),  // Second and third intervals
            (15, 1),  // Only third interval
            (20, 0),  // No intervals
        ]);
    }

    #[test]
    fn test_depths_at_points_empty() {
        let data: Vec<Interval<u32, SmallTwinOl>> = vec![];
        let lapper = Lapper::new(data);
        
        let depths = depths_at_points(&lapper, 0, 10, 5);
        assert_eq!(depths, vec![]);
    }

    #[test]
    fn test_depths_at_points_single_interval() {
        let data = vec![
            Interval { start: 5, stop: 15, val: SmallTwinOl::default() },
        ];
        let lapper = Lapper::new(data);
        
        let depths = depths_at_points(&lapper, 0, 20, 5);
        assert_eq!(depths, vec![
            (0, 0),
            (5, 1),
            (10, 1),
            (15, 0),
            (20, 0),
        ]);
    }


        #[test]
    fn test_depths_at_points_dense_overlap() {
        let data = vec![
            Interval { start: 0, stop: 10, val: SmallTwinOl::default() },
            Interval { start: 0, stop: 5, val: SmallTwinOl::default() },
            Interval { start: 0, stop: 3, val: SmallTwinOl::default() },
        ];
        let lapper = Lapper::new(data);
        
        let depths = depths_at_points(&lapper, 0, 10, 1);
        assert_eq!(depths, vec![
            (0, 3),  // All three intervals
            (1, 3),  // All three intervals
            (2, 3),  // All three intervals
            (3, 2),  // Two intervals remain
            (4, 2),  // Two intervals remain
            (5, 1),  // Only the longest interval
            (6, 1),  // Only the longest interval
            (7, 1),  // Only the longest interval
            (8, 1),  // Only the longest interval
            (9, 1),  // Only the longest interval
            (10, 0), // No intervals
        ]);
    }

    #[test]
    fn test_depths_at_points_staggered_overlap() {
        let data = vec![
            Interval { start: 0, stop: 10, val: SmallTwinOl::default() },
            Interval { start: 2, stop: 8, val: SmallTwinOl::default() },
            Interval { start: 4, stop: 6, val: SmallTwinOl::default() },
        ];
        let lapper = Lapper::new(data);
        
        let depths = depths_at_points(&lapper, 0, 10, 2);
        assert_eq!(depths, vec![
            (0, 1),  // First interval only
            (2, 2),  // First and second intervals
            (4, 3),  // All three intervals
            (6, 2),  // Back to two intervals
            (8, 1),  // Back to one interval
            (10, 0), // No intervals
        ]);
    }

    #[test]
    fn test_depths_at_points_exact_boundaries() {
        let data = vec![
            Interval { start: 5, stop: 10, val: SmallTwinOl::default() },
            Interval { start: 10, stop: 15, val: SmallTwinOl::default() },
            Interval { start: 15, stop: 20, val: SmallTwinOl::default() },
        ];
        let lapper = Lapper::new(data);
        
        let depths = depths_at_points(&lapper, 5, 20, 5);
        assert_eq!(depths, vec![
            (5, 1),   // First interval
            (10, 1),  // Second interval (boundary)
            (15, 1),  // Third interval (boundary)
            (20, 0),  // No intervals
        ]);
    }

    #[test]
    fn test_depths_at_points_non_standard_step() {
        let data = vec![
            Interval { start: 0, stop: 20, val: SmallTwinOl::default() },
            Interval { start: 5, stop: 15, val: SmallTwinOl::default() },
        ];
        let lapper = Lapper::new(data);
        
        let depths = depths_at_points(&lapper, 0, 20, 3);
        assert_eq!(depths, vec![
            (0, 1),   // First interval only
            (3, 1),   // First interval only
            (6, 2),   // Both intervals
            (9, 2),   // Both intervals
            (12, 2),  // Both intervals
            (15, 1),  // First interval only
            (18, 1),  // First interval only
        ]);
    }

    #[test]
    fn test_depths_at_points_query_subset() {
        let data = vec![
            Interval { start: 0, stop: 100, val: SmallTwinOl::default() },
            Interval { start: 20, stop: 80, val: SmallTwinOl::default() },
            Interval { start: 40, stop: 60, val: SmallTwinOl::default() },
        ];
        let lapper = Lapper::new(data);
        
        // Query only middle section
        let depths = depths_at_points(&lapper, 30, 70, 10);
        assert_eq!(depths, vec![
            (30, 2),  // Two intervals
            (40, 3),  // Three intervals
            (50, 3),  // Three intervals
            (60, 2),  // Back to two intervals
            (70, 2),  // Two intervals
        ]);
    }

    #[test]
    fn test_depths_at_points_edge_cases() {
        let data = vec![
            Interval { start: 10, stop: 20, val: SmallTwinOl::default() },
        ];
        let lapper = Lapper::new(data);
        
        // Test empty range
        assert_eq!(depths_at_points(&lapper, 5, 4, 1), vec![]);
        
        // Test zero step
        assert_eq!(depths_at_points(&lapper, 0, 20, 0), vec![]);
        
        // Test single point
        assert_eq!(
            depths_at_points(&lapper, 15, 15, 5),
            vec![(15, 1)]
        );
        
        // Test points beyond interval
        assert_eq!(
            depths_at_points(&lapper, 0, 30, 10),
            vec![(0, 0), (10, 1), (20, 0), (30, 0)]
        );
    }

     #[test]
    fn test_median_and_min_depth_basic() {
        let intervals = vec![
            Interval {
                start: 0,
                stop: 100,
                val: BareMappingOverlap {
                    snpmer_identity: 1.0,
                    ..Default::default()
                }
            },
            Interval {
                start: 25,
                stop: 75,
                val: BareMappingOverlap {
                    snpmer_identity: 1.0,
                    ..Default::default()
                }
            },
        ];
        let lapper = Lapper::new(intervals);
        
        let result = median_and_min_depth_from_lapper(&lapper, 10, 0, 100, 0.9);
        assert!(result.is_some());
        let (min_depth, median_depth) = result.unwrap();
        assert!(min_depth > 0.0);
        assert!(median_depth >= min_depth);

        let result = median_and_min_depth_from_lapper_new(&lapper, 10, 0, 100, 0.9);
        assert!(result.is_some());
        let (min_depth, median_depth) = result.unwrap();
        assert!(min_depth > 0.0);
        assert!(median_depth >= min_depth);

    }

    #[test]
    fn test_median_and_min_depth_empty() {
        let intervals = vec![];
        let lapper = Lapper::new(intervals);
        
        let result = median_and_min_depth_from_lapper(&lapper, 10, 0, 100, 0.9);
        assert_eq!(result, Some((0.0, 0.0)));

        let result = median_and_min_depth_from_lapper_new(&lapper, 10, 0, 100, 0.9);
        assert_eq!(result, Some((0.0, 0.0)));
    }

    #[test]
    fn test_median_and_min_depth_filtered_out() {
        let intervals = vec![
            Interval {
                start: 0,
                stop: 100,
                val: BareMappingOverlap {
                    snpmer_identity: 0.5, // Below cutoff
                    ..Default::default()
                }
            },
        ];
        let lapper = Lapper::new(intervals);
        
        let result = median_and_min_depth_from_lapper(&lapper, 10, 0, 100, 0.9);
        assert_eq!(result, Some((0.0, 0.0)));


        let result = median_and_min_depth_from_lapper_new(&lapper, 10, 0, 100, 0.9);
        assert_eq!(result, Some((0.0, 0.0)));
    }

      #[test]
    fn test_sparse_coverage() {
        let intervals = vec![
            Interval {
                start: 0,
                stop: 20,
                val: create_small_twinol(1.0, true)
            },
            Interval {
                start: 50,
                stop: 70,
                val: create_small_twinol(1.0, true)
            },
            Interval {
                start: 90,
                stop: 100,
                val: create_small_twinol(1.0, true)
            },
        ];
        let lapper = Lapper::new(intervals);
        
        let old_result = median_and_min_depth_from_lapper(&lapper, 5, 0, 100, 0.9);
        let new_result = median_and_min_depth_from_lapper_new(&lapper, 5, 0, 100, 0.9);
        
        assert!(old_result.is_some() && new_result.is_some());
        let (old_min, old_median) = old_result.unwrap();
        let (new_min, new_median) = new_result.unwrap();
        
        assert_relative_eq!(old_min, new_min, epsilon = 1e-10);
        assert_relative_eq!(old_median, new_median, epsilon = 1e-10);
    }

    #[test]
    fn test_maximal_overlap_filtering() {
        let intervals = vec![
            Interval {
                start: 0,
                stop: 100,
                val: create_small_twinol(1.0, true)
            },
        ];
        let lapper = Lapper::new(intervals);
        
        let old_result = median_and_min_depth_from_lapper(&lapper, 10, 0, 100, 0.9);
        let new_result = median_and_min_depth_from_lapper_new(&lapper, 10, 0, 100, 0.9);
        
        assert!(old_result.is_some() && new_result.is_some());
        let (old_min, old_median) = old_result.unwrap();
        let (new_min, new_median) = new_result.unwrap();
        
        assert_relative_eq!(old_min, new_min, epsilon = 1e-10);
        assert_relative_eq!(old_median, new_median, epsilon = 1e-10);
        assert_relative_eq!(old_min, 1.0, epsilon = 1e-10);  // Should only count maximal intervals
    }

    #[test]
    fn test_edge_cases() {
        // Empty intervals
        let empty_lapper = Lapper::new(vec![]);
        let old_empty = median_and_min_depth_from_lapper(&empty_lapper, 10, 0, 100, 0.9);
        let new_empty = median_and_min_depth_from_lapper_new(&empty_lapper, 10, 0, 100, 0.9);
        assert_eq!(old_empty, new_empty);
        assert_eq!(old_empty, Some((0.0, 0.0)));

        // Single point interval
        let point_intervals = vec![
            Interval {
                start: 50,
                stop: 51,
                val: create_small_twinol(1.0, true)
            },
        ];
        let point_lapper = Lapper::new(point_intervals);
        let old_point = median_and_min_depth_from_lapper(&point_lapper, 1, 0, 100, 0.9);
        let new_point = median_and_min_depth_from_lapper_new(&point_lapper, 1, 0, 100, 0.9);
        assert!(old_point.is_some() && new_point.is_some());
        assert_eq!(old_point, new_point);
    }

    #[test]
    fn test_varying_sampling_rates() {
        let intervals = vec![
            Interval {
                start: 0,
                stop: 100,
                val: create_small_twinol(1.0, true)
            },
            Interval {
                start: 25,
                stop: 75,
                val: create_small_twinol(1.0, true)
            },
        ];
        let lapper = Lapper::new(intervals);
        
        for sampling in [1, 5, 10].iter() {
            let old_result = median_and_min_depth_from_lapper(&lapper, *sampling, 0, 100, 0.9);
            let new_result = median_and_min_depth_from_lapper_new(&lapper, *sampling, 0, 100, 0.9);
            
            dbg!(sampling);
            assert!(old_result.is_some() && new_result.is_some());
            let (old_min, old_median) = old_result.unwrap();
            let (new_min, new_median) = new_result.unwrap();
            
            assert_relative_eq!(old_min, new_min, epsilon = 1e-10);
            assert_relative_eq!(old_median, new_median, epsilon = 1e-10);
        }
    }

    #[test]
    fn timing_test_depth_over_points() {
        let mut intervals = vec![];
        for i in 0..1000{
            intervals.push((i as u32, i as u32 + 100, create_small_twinol(1.0, true)));
        }

        for j in 0..1000{
            intervals.push((j as u32 + 500 , j as u32 + 500 + 100, create_small_twinol(1.0, true)));

        }

        for k in 0..1000{
            intervals.push((k as u32, k as u32 + 100, create_small_twinol(1.0, true)));
        }

        let intervals = intervals.into_iter().map(|x| Interval{start: x.0, stop: x.1, val: x.2}).collect::<Vec<_>>();
        let lapper = Lapper::new(intervals);

        let start = std::time::Instant::now();
        let _depths = depths_at_points(&lapper, 0, 1000, 1);
        dbg!(start.elapsed());
    }

    /// Creates a basic TwinRead with evenly spaced minimizers
    fn make_test_read() -> TwinRead {
        TwinRead {
            minimizer_positions: vec![10, 20, 30, 40, 50],
            minimizer_kmers: vec![0.into(), 0.into(), 0.into(), 0.into(), 0.into()],
            k: 21,
            ..Default::default()
        }
    }

    /// Creates a mapping interval with sensible defaults
    fn make_interval(start: u32, stop: u32) -> (u32, u32, f32, bool) {
        (start, stop, 1.0, true)  // Default to 100% identity and maximal overlap
    }

    /// Creates a TwinReadMapping from a list of (start, stop) positions
    /// Optionally override identity and maximal_overlap with full tuples
    fn make_test_mapping<T>(intervals: T) -> TwinReadMapping 
    where 
        T: IntoIterator<Item = (u32, u32, f32, bool)> + Clone
    {

        let lapper_intervals_map = intervals.clone()
            .into_iter()
            .enumerate()
            .map(|(i, (start, stop, identity, _))| Interval {
                start,
                stop,
                val: BareMappingOverlap {
                    snpmer_identity: identity,
                },
            })
            .collect();

        let lapper_intervals = intervals
            .into_iter()
            .enumerate()
            .map(|(i, (start, stop, identity, _))| Interval {
                start,
                stop,
                val: SmallTwinOl {
                    query_id: i as u32,
                    snpmer_identity: identity,
                    ..Default::default()
                },
            })
            .collect();

        TwinReadMapping {
            tr_index: 0,
            mapping_info: MappingInfo {
                median_depth: 0.0,
                minimum_depth: 0.0,
                max_mapping_boundaries: Some(Lapper::new(lapper_intervals_map)),
                max_alignment_boundaries: Some(Lapper::new(lapper_intervals)),
                present: true,
                length: 1000,
            },
            all_intervals: vec![],
        }
    }

    #[test]
    fn test_basic_coverage() {
        let mut read = make_test_read();
        // Single interval covering whole region
        let mapping = make_test_mapping(vec![make_interval(0, 100)]);
        
        populate_depth_from_map_info(&mut read, &mapping, 0, 100);

        assert!(read.min_depth_multi.is_some());
        assert!(read.median_depth.is_some());
        if let Some(min_depths) = read.min_depth_multi {
            assert_relative_eq!(min_depths[0], 1.0, epsilon = 0.001);
        }
    }

    #[test]
    fn test_varying_identity() {
        let mut read = make_test_read();
        let mapping = make_test_mapping(vec![
            (0, 100, 1.0, true),    // 100% identity
            (0, 100, 0.998, true),   // 95% identity
            (0, 100, 0.991, true),   // 90% identity
        ]);
        
        populate_depth_from_map_info(&mut read, &mapping, 0, 100);

        if let Some(min_depths) = read.min_depth_multi {
            assert_relative_eq!(min_depths[0], 3.0, epsilon = 0.001);
            assert_relative_eq!(min_depths[1], 2.0, epsilon = 0.001);
            assert_relative_eq!(min_depths[2], 1.0, epsilon = 0.001);
        }
    }

    #[test]
    fn test_varying_identity_adaptive() {
        let mut read = make_test_read();
        let mapping = make_test_mapping(vec![
            (0, 100, 1.0, true),    // 100% identity
            (0, 100, 0.998, true),   // 95% identity
            (0, 100, 0.997, true),   // 95% identity
            (0, 100, 0.996, true),   // 95% identity
            (0, 100, 0.995, true),   // 95% identity
            (0, 100, 0.994, true),   // 95% identity
            (0, 100, 0.993, true),   // 95% identity
            (0, 100, 0.992, true),   // 95% identity
            (0, 100, 0.991, true),   // 90% identity
            (0, 100, 0.991, true),   // 90% identity
            (0, 100, 0.991, true),   // 90% identity
            (0, 100, 0.991, true),   // 90% identity
            (0, 100, 0.991, true),   // 90% identity
            (0, 100, 0.991, true),   // 90% identity
            (0, 100, 0.991, true),   // 90% identity
        ]);
        
        populate_depth_from_map_info(&mut read, &mapping, 0, 100);

        dbg!(&read.snpmer_id_threshold);
        assert!(read.snpmer_id_threshold.unwrap() < 99.9);
        assert!(read.snpmer_id_threshold.unwrap() > 1.0);
    }

    #[test]
    fn test_gaps_in_coverage() {
        let mut read = make_test_read();
        let mapping = make_test_mapping(vec![
            make_interval(0, 25),     // First quarter
            make_interval(75, 100),   // Last quarter
        ]);
        
        populate_depth_from_map_info(&mut read, &mapping, 0, 100);

        if let Some(min_depths) = read.min_depth_multi {
            assert!(min_depths[0] < 1.0);  // Should be less than 1 due to gaps
        }
    }
}
