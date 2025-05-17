use std::process::Command;
use assert_cmd::prelude::*;
use std::fs;
use std::path::Path;
use tempfile::TempDir;


#[test]
fn test_missing_input() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("myloasm")?;
    cmd.assert()
        .failure();
    Ok(())
}

#[test]
fn test_nanopore_basic_run_2kb_plasmid() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("output");
    
    let test_input = Path::new("tests/reads/2kb_plas.fq");

    let mut cmd = Command::cargo_bin("myloasm")?;
    cmd.arg(test_input.to_str().unwrap())
        .arg("-o")
        .arg(output_dir.to_str().unwrap())
        .arg("-t")
        .arg("20"); // Use fewer threads for testing
    
    cmd.assert().success();
    
    // Check output directory exists and contains expected files
    assert!(output_dir.exists());
    assert!(output_dir.join("assembly_primary.fa").exists());

    //Check the assembly file has > 1.8kb and <2kb and has only one contig
    let assembly_file = output_dir.join("assembly_primary.fa");
    let assembly_content = fs::read_to_string(&assembly_file)?;
    let contigs: Vec<&str> = assembly_content.lines().filter(|line| line.starts_with('>')).collect();
    println!("Contigs found: {:?}", contigs);
    assert_eq!(contigs.len(), 1, "Expected 1 contig, found {}", contigs.len());
    let contig_length = assembly_content.lines().filter(|line| !line.starts_with('>')).collect::<String>().len();
    assert!(contig_length > 1800 && contig_length < 2000, "Expected contig length between 1800 and 2000, found {}", contig_length);
    
    Ok(())
}

#[test]
fn test_nanopore_basic_run_48kb_plasmid() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("output");
    
    let test_input = Path::new("tests/reads/40kb_plas.fq");

    let mut cmd = Command::cargo_bin("myloasm")?;
    cmd.arg(test_input.to_str().unwrap())
        .arg("-o")
        .arg(output_dir.to_str().unwrap())
        .arg("-t")
        .arg("20"); // Use fewer threads for testing
    
    cmd.assert().success();
    
    // Check output directory exists and contains expected files
    assert!(output_dir.exists());
    assert!(output_dir.join("assembly_primary.fa").exists());

    //Check the assembly file has > 1.8kb and <2kb and has only one contig
    let assembly_file = output_dir.join("assembly_primary.fa");
    let assembly_content = fs::read_to_string(&assembly_file)?;
    let contigs: Vec<&str> = assembly_content.lines().filter(|line| line.starts_with('>')).collect();
    println!("Contigs found: {:?}", contigs);
    assert_eq!(contigs.len(), 1, "Expected 1 contig, found {}", contigs.len());
    let contig_length = assembly_content.lines().filter(|line| !line.starts_with('>')).collect::<String>().len();
    assert!(contig_length > 40000 && contig_length < 60000, "Expected contig length between 40000 and 60000, found {}", contig_length);
    
    Ok(())
}
