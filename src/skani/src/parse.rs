use crate::cmd_line::*;
use crate::params::*;
use crate::regression;
use clap::parser::ArgMatches;
use log::LevelFilter;
use log::*;
use std::fs;
use std::fs::File;
use std::io::{prelude::*, BufReader};

pub fn parse_params(matches: &ArgMatches) -> (SketchParams, CommandParams) {
    let mode;
    let matches_subc;
    match matches.subcommand_name() {
        Some(SKETCH_STRING) => {
            mode = Mode::Sketch;
            matches_subc = matches.subcommand_matches(SKETCH_STRING).unwrap()
        }
        Some(DIST_STRING) => {
            mode = Mode::Dist;
            matches_subc = matches.subcommand_matches(DIST_STRING).unwrap()
        }
        Some(TRIANGLE_STRING) => {
            mode = Mode::Triangle;
            matches_subc = matches.subcommand_matches(TRIANGLE_STRING).unwrap()
        }
        Some(SEARCH_STRING) => {
            mode = Mode::Search;
            matches_subc = matches.subcommand_matches(SEARCH_STRING).unwrap();
        }
        _ => {
            panic!()
        } // Either no subcommand or one not tested for...
    }

    let threads = matches_subc.value_of("t").unwrap();
    let threads = threads.parse::<usize>().unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    simple_logging::log_to_stderr(LevelFilter::Info);
    if matches_subc.is_present("v") {
        simple_logging::log_to_stderr(LevelFilter::Debug);
    }
    if matches_subc.is_present("trace") {
        simple_logging::log_to_stderr(LevelFilter::Trace);
    }

    if mode == Mode::Search {
        return parse_params_search(matches_subc);
    }

    let amino_acid;
    if matches_subc.is_present("aai") {
        warn!("Amino acid mode (AAI) detected. This mode is not stable.");
        amino_acid = true;
    } else {
        amino_acid = false;
    }

    let rescue_small;
    if mode == Mode::Triangle || mode == Mode::Dist{
        if matches_subc.is_present(FAST_SMALL) || matches_subc.is_present(MODE_SMALL_GENOMES) {
            rescue_small = false;
        }
        else{
            rescue_small = true;
        }
    }
    else{
        rescue_small = false;
    }

    let mut ref_files: Vec<String>;
    let mut ref_file_list = None;
    let mut sparse = false;
    if mode == Mode::Triangle {
        sparse = matches_subc.is_present("sparse");
    }
    if mode == Mode::Triangle || mode == Mode::Sketch {
        if let Some(values) = matches_subc.values_of("fasta_files") {
            ref_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.value_of("fasta_list") {
            ref_files = vec![];
            ref_file_list = Some(values);
        } else {
            error!("No reference inputs found.");
            std::process::exit(1);
        }
    } else if mode == Mode::Dist {
        if let Some(values) = matches_subc.values_of("reference") {
            ref_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.value_of("reference list file") {
            ref_files = vec![];
            ref_file_list = Some(values);
        } else if let Some(values) = matches_subc.values_of("references") {
            ref_files = values.map(|x| x.to_string()).collect();
        } else {
            error!("No reference inputs found.");
            std::process::exit(1);
        }
    } else {
        panic!("PATH TODO");
    }

    if ref_file_list.is_some() {
        let ref_file_list = ref_file_list.unwrap();
        let file = File::open(ref_file_list).expect("-l specified file could not be opened properly. Make sure this file exists. Exiting.");
        let reader = BufReader::new(file);
        let mut temp_vec = vec![];

        for line in reader.lines() {
            temp_vec.push(line.unwrap().trim().to_string());
        }
        ref_files = temp_vec;
    }

    let mut query_files = vec![];
    let mut query_file_list = None;
    let mut max_results = usize::MAX;

    if mode == Mode::Dist {
        max_results = matches_subc
            .value_of("n")
            .unwrap_or("1000000000000")
            .parse::<usize>()
            .unwrap();
        if let Some(values) = matches_subc.values_of("query") {
            query_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.values_of("queries") {
            query_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.value_of("query list file") {
            query_file_list = Some(values)
        }
    }

    if query_file_list.is_some() {
        let query_file_list = query_file_list.unwrap();
        let file = File::open(query_file_list).unwrap();
        let reader = BufReader::new(file);
        let mut temp_vec = vec![];

        for line in reader.lines() {
            temp_vec.push(line.unwrap().trim().to_string());
        }
        query_files = temp_vec;
    }

    let def_k = if amino_acid { DEFAULT_K_AAI } else { DEFAULT_K };
    let def_c = if amino_acid { DEFAULT_C_AAI } else { DEFAULT_C };
    let k = matches_subc
        .value_of("k")
        .unwrap_or(def_k)
        .parse::<usize>()
        .unwrap();

    let use_syncs = false;
    let mut c = matches_subc
        .value_of("c")
        .unwrap_or(def_c)
        .parse::<usize>()
        .unwrap();

    let mut marker_c = matches_subc
        .value_of("marker_c")
        .unwrap_or(MARKER_C_DEFAULT)
        .parse::<usize>()
        .unwrap();

    if matches_subc.is_present(MODE_FAST) &&
        matches_subc.is_present(MODE_SLOW){
            panic!("Both --slow and --fast were set. This is not allowed.");
    }
    if matches_subc.is_present(MODE_FAST){
        if matches_subc.is_present("c"){
            warn!("-c value is set but --fast is also set. Using --fast mode instead (-c 200)");
        }
        c = FAST_C
    }
    if matches_subc.is_present(MODE_SLOW){
        if matches_subc.is_present("c"){
            warn!("-c value is set but --slow is also set. Using --slow mode instead (-c 30)");
        }
        c = SLOW_C
    }
    if matches_subc.is_present(MODE_MEDIUM){
        if matches_subc.is_present("c"){
            warn!("-c value is set but --fast is also set. Using --medium mode instead (-c 70)");
        }
        c = MEDIUM_C
    }
    if mode == Mode::Triangle || mode == Mode::Dist{
        if matches_subc.is_present(MODE_SMALL_GENOMES){
            if matches_subc.is_present("c") || matches_subc.is_present("marker_c") {
                warn!("-c or -m value is set but --small-genomes is also set. Using -c 30 and -m 200 instead.");
            }
            c = SLOW_C;
            marker_c = SMALL_M;
        }
    }

    let min_aligned_frac;
    let est_ci;
    let detailed_out;
    if mode != Mode::Sketch {
        let def_maf = if amino_acid {
            D_FRAC_COVER_CUTOFF_AA
        } else {
            D_FRAC_COVER_CUTOFF
        };
        min_aligned_frac = matches_subc
            .value_of(MIN_ALIGN_FRAC)
            .unwrap_or(def_maf)
            .parse::<f64>()
            .unwrap()
            / 100.;
        est_ci = matches_subc.is_present(CONF_INTERVAL);
        detailed_out = matches_subc.is_present(DETAIL_OUT);
    } else {
        min_aligned_frac = 0.;
        est_ci = false;
        detailed_out = false;
    }

    let out_file_name;
    if mode == Mode::Triangle {
        out_file_name = matches_subc.value_of("output").unwrap_or("").to_string();
    } else if mode == Mode::Sketch {
        out_file_name = matches_subc
            .value_of("output sketch folder")
            .unwrap_or("")
            .to_string();
    } else if mode == Mode::Dist {
        out_file_name = matches_subc.value_of("output").unwrap_or("").to_string();
    } else {
        panic!("Mode doesn't exist");
    }

    let mut screen_val = 0.;
    let mut robust = false;
    let mut median = false;
    if mode == Mode::Triangle || mode == Mode::Dist {
        screen_val = matches_subc
            .value_of("s")
            .unwrap_or("0.00")
            .parse::<f64>()
            .unwrap()
            / 100.;
    }
    if mode == Mode::Triangle || mode == Mode::Search || mode == Mode::Dist {
        robust = matches_subc.is_present("robust");
        median = matches_subc.is_present("median");
    }

    let sketch_params = SketchParams::new(marker_c, c, k, use_syncs, amino_acid);

    let mut refs_are_sketch = !ref_files.is_empty();
    for ref_file in ref_files.iter() {
        if !ref_file.contains(".sketch")
            && !ref_file.contains(".marker")
            && !ref_file.contains("markers.bin")
        {
            refs_are_sketch = false;
            break;
        }
    }

    let mut queries_are_sketch = !query_files.is_empty();
    for query_file in query_files.iter() {
        if !query_file.contains(".sketch") && !query_file.contains("markers.bin") {
            queries_are_sketch = false;
            break;
        }
    }

    let individual_contig_q;
    let individual_contig_r;

    if mode == Mode::Triangle {
        if matches_subc.is_present("individual contig") {
            individual_contig_q = true;
            individual_contig_r = true;
        } else {
            individual_contig_q = false;
            individual_contig_r = false;
        }
    } else if mode == Mode::Dist {
        individual_contig_q = matches_subc.is_present("individual contig query");
        individual_contig_r = matches_subc.is_present("individual contig ref");
    } else if mode == Mode::Sketch{
        individual_contig_q = false;
        individual_contig_r = if matches_subc.is_present("individual contig"){ true } else{ false };
    }
    else{
        individual_contig_q = false;
        individual_contig_r = false;

    }

    let full_matrix;
    let diagonal;
    if mode == Mode::Triangle {
        full_matrix = matches_subc.is_present(FULL_MAT);
        diagonal = matches_subc.is_present(DIAG);
    } else {
        full_matrix = false;
        diagonal = false;
    }

    let screen;
    if mode == Mode::Dist {
        if (query_files.len() > FULL_INDEX_THRESH || individual_contig_q) && !matches_subc.is_present(NO_FULL_INDEX) {
            screen = true;
        } else {
            screen = false;
        }
    } else if mode == Mode::Triangle {
        screen = true;
    } else {
        screen = false;
    }

    let learned_ani;
    if mode == Mode::Sketch{
        learned_ani = false;
    }
    else if matches_subc.is_present(NO_LEARNED_ANI){
        learned_ani = false;
    }
    else{
        learned_ani = regression::use_learned_ani(c, individual_contig_q, individual_contig_r, median);
    }

    let mut distance = false;
    if mode == Mode::Triangle{
        if matches_subc.is_present(DISTANCE_OUT){
            distance = true;
        }
    }

    let command_params = CommandParams {
        screen,
        screen_val,
        mode,
        out_file_name,
        ref_files,
        query_files,
        refs_are_sketch,
        queries_are_sketch,
        robust,
        median,
        sparse,
        full_matrix,
        diagonal,
        max_results,
        individual_contig_q,
        individual_contig_r,
        min_aligned_frac,
        keep_refs: false,
        est_ci,
        learned_ani,
        detailed_out,
        distance,
        rescue_small,
    };

    (sketch_params, command_params)
}

pub fn parse_params_search(matches_subc: &ArgMatches) -> (SketchParams, CommandParams) {
    let mode = Mode::Search;
    let out_file_name = matches_subc.value_of("output").unwrap_or("").to_string();

    let mut query_files = vec![];
    let mut query_file_list = None;
    let max_results = matches_subc
        .value_of("n")
        .unwrap_or("10000000")
        .parse::<usize>()
        .unwrap();
    if let Some(values) = matches_subc.values_of("query") {
        query_files = values.map(|x| x.to_string()).collect();
    } else if let Some(values) = matches_subc.values_of("queries") {
        query_files = values.map(|x| x.to_string()).collect();
    } else if let Some(values) = matches_subc.value_of("query list file") {
        query_file_list = Some(values);
    }
    if query_file_list.is_some() {
        let query_file_list = query_file_list.unwrap();
        let file = File::open(query_file_list).unwrap();
        let reader = BufReader::new(file);
        let mut temp_vec = vec![];

        for line in reader.lines() {
            temp_vec.push(line.unwrap().trim().to_string());
        }
        query_files = temp_vec;
    }
    let ref_folder = matches_subc.value_of("sketched database folder").unwrap();
    let paths =
        fs::read_dir(ref_folder).expect("Issue with folder specified by -d option; exiting");
    let ref_files = paths
        .into_iter()
        .map(|x| x.unwrap().path().to_str().unwrap().to_string())
        .collect();
    let refs_are_sketch = true;

    let mut queries_are_sketch = !query_files.is_empty();
    for query_file in query_files.iter() {
        if !query_file.contains(".sketch") && !query_file.contains("markers.bin") {
            queries_are_sketch = false;
            break;
        }
    }

    let robust = matches_subc.is_present("robust");
    let median = matches_subc.is_present("median");
    let sparse = false;

    let screen_val = matches_subc
        .value_of("s")
        .unwrap_or("0.00")
        .parse::<f64>()
        .unwrap()
        / 100.;
    let screen;
    let individual_contig_q = matches_subc.is_present("individual contig query");
    if (query_files.len() > FULL_INDEX_THRESH || individual_contig_q) && !matches_subc.is_present(NO_FULL_INDEX) {
        screen = true;
    } else {
        screen = false;
    }

    let min_aligned_frac = matches_subc
        .value_of(MIN_ALIGN_FRAC)
        .unwrap_or("-100.0")
        .parse::<f64>()
        .unwrap()
        / 100.;
    let keep_refs = matches_subc.is_present(KEEP_REFS);
    let est_ci = matches_subc.is_present(CONF_INTERVAL);
    let detailed_out = matches_subc.is_present(DETAIL_OUT);


    let learned_ani;
    if matches_subc.is_present(NO_LEARNED_ANI){
        learned_ani = false;
    }
    else{
        learned_ani = true;
    }

    let command_params = CommandParams {
        screen,
        screen_val,
        mode,
        out_file_name,
        ref_files,
        query_files,
        refs_are_sketch,
        queries_are_sketch,
        robust,
        median,
        sparse,
        full_matrix: false,
        diagonal: false,
        max_results,
        individual_contig_q,
        individual_contig_r: false,
        min_aligned_frac,
        keep_refs,
        est_ci,
        learned_ani,
        detailed_out,
        distance: false,
        rescue_small: false
    };

    if command_params.ref_files.is_empty() {
        error!("No valid reference fastas or sketches found.");
        std::process::exit(1)
    }



    (SketchParams::default(), command_params)
}
