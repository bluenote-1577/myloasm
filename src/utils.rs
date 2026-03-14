use memory_stats::memory_stats;

pub fn log_memory_usage(info: bool, message: &str) {
    if let Some(usage) = memory_stats() {
        if info{
            log::info!(
                "{} --- Memory usage: {:.2} GB",
                message,
                usage.physical_mem as f64 / 1_000_000_000.
            );
        }
        else{
            log::debug!(
                "{} --- Memory usage: {:.2} GB",
                message,
                usage.physical_mem as f64 / 1_000_000_000.
            );
        }
    }
    else{
        log::info!("Memory usage: unknown (WARNING)");
    }
}

#[inline]
pub fn div_rounded(a: usize, b: usize) -> usize {
    (a + b / 2) / b
}



pub fn first_word(s: &str) -> String{
    s.split_whitespace().next().unwrap_or(s).to_string()
}

pub fn get_nx_from_vec(vec: &[usize], n_vec: &[usize]) -> String{
    let mut sorted_vec = vec.to_vec();
    sorted_vec.sort_unstable_by(|a, b| b.cmp(a)); // Sort in descending order
    let total: usize = sorted_vec.iter().map(|x| *x as usize).sum();
    let mut nx_values = Vec::new();

    for &n in n_vec {
        let mut cumulative_sum = 0;
        let threshold = total * n / 100;
        for value in &sorted_vec {
            cumulative_sum += value;
            if cumulative_sum >= threshold {
                nx_values.push(value);
                break;
            }
        }
    }

    format!("N{}: {}, N{}: {}, N{}: {}", n_vec[0], nx_values[0], n_vec[1], nx_values[1], n_vec[2], nx_values[2])
}