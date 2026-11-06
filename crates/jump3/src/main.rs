//! jump3: the power-of-three jump self-test and a jump of any ring stencil.
//!
//! usage: jump3 --self-test [--timing]
//!        jump3 --ring D --stencil s:c,s:c --ks k,k [--direct-max K] [--case NAME] [--step TEXT]
//!
//! Prints the same HBP rows, byte for byte, as python/matrix_tricks/jump3_ntt.py
//! (timing rows excepted: they appear only under --timing). Output json=0.
//! Integer-only, Rust 1.81.
#![deny(clippy::float_arithmetic)]

use jump3::{
    cases, circulant, cube_ladder, eigenvalues, head_row, matvec, ntt_jump, parse_stencil,
    rows_text, run_case, start_state, step, verify_row, Case, Flags, Ops, Ring,
};
use std::time::Instant;

fn timing_rows() -> Result<Vec<String>, String> {
    let case = &cases()[0];
    let ring = Ring::new(case.d, &case.stencil)?;
    let lam = eigenvalues(&ring, &mut Ops::default());
    let x0 = start_state(case.d);
    let steps = 3u128.pow(case.direct_max);
    let t0 = Instant::now();
    let mut x = x0.clone();
    for _ in 0..steps {
        x = step(&ring, &x, &mut Ops::default());
    }
    let direct_ns = t0.elapsed().as_nanos();
    let mut rows = Vec::new();
    for &k in &case.ks {
        let n = 3u128.pow(k);
        let t0 = Instant::now();
        let ladder = cube_ladder(circulant(&ring), k, &mut Ops::default());
        matvec(&ladder[k as usize].0, &x0, &mut Ops::default());
        let cube_ns = t0.elapsed().as_nanos();
        let t0 = Instant::now();
        ntt_jump(&ring, &lam, &x0, n, &mut Ops::default());
        let ntt_ns = t0.elapsed().as_nanos();
        let need = direct_ns.checked_mul(n).map(|v| v / steps);
        let (need_s, years_s) = match need {
            Some(v) => (
                v.to_string(),
                (v / (365 * 24 * 3600 * 1_000_000_000)).to_string(),
            ),
            None => ("EXCEEDS_U128".to_string(), "EXCEEDS_U128".to_string()),
        };
        rows.push(format!(
            "JUMP3_TIMING|lane=rust|case={}|k={k}|n={n}|cube_ns={cube_ns}|ntt_ns={ntt_ns}|direct_ns_for_3^{}={direct_ns}|\
             direct_ns_needed={need_s}|direct_years_needed={years_s}|basis=direct_rate_at_3^{}_steps|json=0",
            case.name, case.direct_max, case.direct_max
        ));
    }
    Ok(rows)
}

fn arg_value(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1).cloned())
}

fn run() -> Result<bool, String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut out = vec![head_row()];
    let mut flags = Flags::default();
    if args.iter().any(|a| a == "--self-test") {
        let all = cases();
        for case in &all {
            run_case(case, &mut out, &mut flags)?;
        }
        let (row, ok) = verify_row(&flags, all.len());
        out.push(row);
        if args.iter().any(|a| a == "--timing") {
            out.extend(timing_rows()?);
        }
        print!("{}", rows_text(&out));
        return Ok(ok);
    }
    let (Some(ring), Some(stencil), Some(ks)) = (
        arg_value(&args, "--ring"),
        arg_value(&args, "--stencil"),
        arg_value(&args, "--ks"),
    ) else {
        return Err(
            "usage: jump3 --self-test [--timing] | --ring D --stencil s:c,s:c --ks k,k \
                    [--direct-max K] [--case NAME] [--step TEXT]"
                .to_string(),
        );
    };
    let case = Case {
        name: arg_value(&args, "--case").unwrap_or_else(|| "custom".to_string()),
        d: ring.parse().map_err(|e| format!("--ring: {e}"))?,
        stencil: parse_stencil(&stencil)?,
        label: arg_value(&args, "--step").unwrap_or_else(|| "stencil".to_string()),
        ks: ks
            .split(',')
            .map(|k| k.parse::<u32>().map_err(|e| format!("--ks: {e}")))
            .collect::<Result<_, _>>()?,
        direct_max: arg_value(&args, "--direct-max").map_or(Ok(0), |v| {
            v.parse().map_err(|e| format!("--direct-max: {e}"))
        })?,
    };
    run_case(&case, &mut out, &mut flags)?;
    let (row, ok) = verify_row(&flags, 1);
    out.push(row);
    print!("{}", rows_text(&out));
    Ok(ok)
}

fn main() {
    match run() {
        Ok(true) => {}
        Ok(false) => std::process::exit(1),
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    }
}
