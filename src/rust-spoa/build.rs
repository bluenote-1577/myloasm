// build.rs
extern crate cc;
extern crate cmake;
use cmake::Config;
use std::env;

fn main() {
    // First build zlib
    // Then build spoa
    let spoa_dst = Config::new("src/spoa")
        .define("CMAKE_BUILD_TYPE", "Release")
        .build();

    let out_dir = env::var("OUT_DIR").unwrap();


    cc::Build::new()
        .cpp(true)
        .shared_flag(false)
        .static_flag(true)
        .flag_if_supported("-O3")
        .flag_if_supported("-D_GNU_SOURCE")
        .flag_if_supported("-Wall")
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-Isrc/spoa/include")
        .flag_if_supported(&format!("-L{}/lib64 -L{}/lib", &out_dir, &out_dir))
        .file("src/poa_func.cpp")
        .compile("poa_func");


    // Explicitly specify the library paths
    println!("cargo:rustc-link-search=native={}/lib", spoa_dst.display());
    println!("cargo:rustc-link-search=native={}/lib64", spoa_dst.display());
    
    // Specify the libraries to link
    println!("cargo:rustc-link-lib=static=poa_func");
    println!("cargo:rustc-link-lib=static=spoa");
    println!("cargo:rustc-link-lib=stdc++");

    println!("cargo:rustc-flags=-L {}/lib64/ -L {}/lib/", &out_dir, &out_dir);
}