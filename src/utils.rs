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

