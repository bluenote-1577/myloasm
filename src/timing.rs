use memory_stats::memory_stats;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::{Duration, Instant};

struct StepRecord {
    name: String,
    duration: Duration,
    ram_gb: Option<f64>,
}

pub struct PipelineTimer {
    total_start: Instant,
    steps: Vec<StepRecord>,
}

impl PipelineTimer {
    pub fn new() -> Self {
        PipelineTimer {
            total_start: Instant::now(),
            steps: Vec::new(),
        }
    }

    pub fn measure<F: FnOnce() -> R, R>(&mut self, name: &str, f: F) -> R {
        let t = Instant::now();
        let result = f();
        let duration = t.elapsed();
        let ram_gb = memory_stats().map(|m| m.physical_mem as f64 / 1_073_741_824.0);
        self.steps.push(StepRecord {
            name: name.to_string(),
            duration,
            ram_gb,
        });
        result
    }

    pub fn write_tsv(&self, path: &Path) {
        let total = self.total_start.elapsed();
        let sum: Duration = self.steps.iter().map(|s| s.duration).sum();
        let other = total.saturating_sub(sum);

        match File::create(path) {
            Ok(f) => {
                let mut w = BufWriter::new(f);
                let _ = writeln!(w, "step\tduration_s\tduration_percent\tram_after_gb");
                for s in &self.steps {
                    let ram_str = s.ram_gb.map_or("-".to_string(), |r| format!("{:.2}", r));
                    let percent = s.duration.as_secs_f64() / total.as_secs_f64() * 100.0;
                    let _ = writeln!(
                        w,
                        "{}\t{:.3}\t{:.2}\t{}",
                        s.name,
                        s.duration.as_secs_f64(),
                        percent,
                        ram_str
                    );
                }
                let other_percent = other.as_secs_f64() / total.as_secs_f64() * 100.0;
                let _ = writeln!(
                    w,
                    "OTHER\t{:.3}\t{:.2}\t-",
                    other.as_secs_f64(),
                    other_percent
                );
                let _ = writeln!(w, "TOTAL\t{:.3}\t-", total.as_secs_f64());
            }
            Err(e) => log::warn!("Could not write timing TSV to {:?}: {}", path, e),
        }

        // log::debug!("Pipeline timing summary:");
        // for s in &self.steps {
        //     let percent = s.duration.as_secs_f64() / total.as_secs_f64() * 100.0;
        //     log::debug!("  {:<45} {:.3}s {:.2}%", s.name, s.duration.as_secs_f64(), percent);
        // }
        // let other_percent = other.as_secs_f64() / total.as_secs_f64() * 100.0;
        // log::debug!("  {:<45} {:.3}s {:.2}%", "OTHER", other.as_secs_f64(), other_percent);
        // log::debug!("  {:<45} {:.3}s {:.2}%", "TOTAL", total.as_secs_f64(), 100.0);
    }
}
