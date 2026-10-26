//! Write the rainbon operator in the standard exchange format
//! (n / weights / magnetic=potential / "i j z" electric entries, j >= i), so
//! `cluster`, `cheb`, `jumpto` and `sturmrs` run on it unchanged.
//! usage: cloudb L harmonic|linear OUT_MATRIX     (json=0)
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic)]

use rainbon::{hops, potential, Lattice, Potential};
use std::collections::BTreeMap;
use std::io::Write;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: cloudb L harmonic|linear OUT_MATRIX");
        std::process::exit(1);
    }
    let l: i64 = a[1].parse().unwrap_or(0);
    let kind = match a[2].as_str() {
        "harmonic" => Potential::Harmonic,
        "linear" => Potential::Linear,
        _ => {
            println!("CLOUDB|error=potential_must_be_harmonic_or_linear|json=0");
            std::process::exit(1);
        }
    };
    if !(3..=12).contains(&l) {
        println!("CLOUDB|error=L_out_of_range|json=0");
        std::process::exit(1);
    }
    let t0 = std::time::Instant::now();
    let lat = Lattice { l };
    let n = lat.dim();
    let f = std::fs::File::create(&a[3]).unwrap_or_else(|e| {
        println!("CLOUDB|error=write:{e}|json=0");
        std::process::exit(1)
    });
    let mut out = std::io::BufWriter::new(f);
    let mut ok = writeln!(out, "{n}").is_ok();
    ok &= writeln!(out, "{}", vec!["1"; n].join(" ")).is_ok();
    let pots: Vec<String> = (0..n).map(|i| potential(&lat, &lat.decode(i), kind).to_string()).collect();
    ok &= writeln!(out, "{}", pots.join(" ")).is_ok();
    let mut nnz = 0usize;
    for i in 0..n {
        let v = lat.decode(i);
        let mut row: BTreeMap<usize, i64> = BTreeMap::new();
        row.insert(i, 18);
        for j in hops(&lat, &v) {
            *row.entry(j).or_insert(0) -= 1;
        }
        for (&j, &z) in &row {
            if z != 0 {
                nnz += 1;
                if j >= i {
                    ok &= writeln!(out, "{i} {j} {z}").is_ok();
                }
            }
        }
    }
    ok &= out.flush().is_ok();
    println!(
        "CLOUDB|L={l}|potential={}|states={n}|nnz={nnz}|write={}|ms={}|json=0",
        a[2],
        if ok { "ok" } else { "failed" },
        t0.elapsed().as_millis()
    );
}
