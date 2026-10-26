//! Ground state of the rainbon (three-preon cloud) and its charge radius.
//! Builds H = JE*E + JB*M in memory, finds the ground vector by an integer
//! Chebyshev filter (doubling blocks; b = Gershgorin, a restarted at the
//! current Rayleigh quotient + GUESS/2), then reports, in exact big integers:
//!   E0 (Rayleigh quotient) and a rigorous residual radius eps
//!   (some eigenvalue lies in [E0 - eps, E0 + eps]),
//!   <r^2> = <R>/9 of that vector (charge radius squared, lattice units^2).
//! The radius is a measurement of the computed vector (its error is controlled
//! by convergence across step counts), not a certificate.
//! usage: cloud L harmonic|linear JE JB STEPS GUESS_MICRO    (micro = 1e-6 of
//! the energy unit in which JE, JB are given)                  json=0
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic)]

use num_bigint::BigInt;
use rainbon::{hops, potential, Lattice, Potential};
use std::collections::BTreeMap;

const VB: u32 = 40;

fn top_bits(v: &[i128]) -> u32 {
    128 - v.iter().map(|x| x.unsigned_abs()).max().unwrap_or(0).leading_zeros()
}
fn shr(v: &mut [i128], s: u32) {
    for x in v.iter_mut() {
        *x >>= s;
    }
}

struct Op {
    rows: Vec<Vec<(u32, i64)>>,
}
impl Op {
    fn apply(&self, v: &[i128], out: &mut [i128]) {
        for (o, row) in out.iter_mut().zip(self.rows.iter()) {
            let mut s = 0i128;
            for &(j, c) in row {
                s += i128::from(c) * v[j as usize];
            }
            *o = s;
        }
    }
}

fn micro(num: &BigInt, den: &BigInt) -> BigInt {
    num * BigInt::from(1_000_000) / den
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 7 {
        eprintln!("usage: cloud L harmonic|linear JE JB STEPS GUESS_MICRO");
        std::process::exit(1);
    }
    let t0 = std::time::Instant::now();
    let l: i64 = a[1].parse().unwrap_or(0);
    let kind = if a[2] == "linear" { Potential::Linear } else { Potential::Harmonic };
    let je: i64 = a[3].parse().unwrap_or(1);
    let jb: i64 = a[4].parse().unwrap_or(0);
    let steps: usize = a[5].parse().unwrap_or(2000);
    let guess: i128 = a[6].parse().unwrap_or(1_000_000);
    if !(3..=12).contains(&l) {
        println!("CLOUD|error=L_out_of_range|json=0");
        std::process::exit(1);
    }
    let lat = Lattice { l };
    let n = lat.dim();
    let mut rows = Vec::with_capacity(n);
    let mut robs = Vec::with_capacity(n);
    for i in 0..n {
        let v = lat.decode(i);
        let mut row: BTreeMap<usize, i64> = BTreeMap::new();
        row.insert(i, 18 * je + jb * potential(&lat, &v, kind));
        for j in hops(&lat, &v) {
            *row.entry(j).or_insert(0) -= je;
        }
        rows.push(row.into_iter().filter(|&(_, c)| c != 0).map(|(j, c)| (u32::try_from(j).unwrap_or(0), c)).collect::<Vec<_>>());
        robs.push(lat.radius_obs(&v));
    }
    let op = Op { rows };
    let beta: i128 = op.rows.iter().map(|r| r.iter().map(|&(_, c)| i128::from(c.abs())).sum::<i128>()).max().unwrap_or(0) + 1;
    let gs = filter(&op, n, beta, guess, steps, None);
    let ex = filter(&op, n, beta, guess, steps, Some(&gs));
    let (e0, eps0, r20) = verify(&op, &gs, &robs);
    let (e1, eps1, r21) = verify(&op, &ex, &robs);
    println!(
        "CLOUD|L={l}|potential={}|JE={je}|JB={jb}|states={n}|steps={steps}|E0={e0}/1000000|eps0_up={eps0}/1000000|r2_ground={r20}/1000000|E1={e1}/1000000|eps1_up={eps1}/1000000|r2_excited={r21}/1000000|gap={}/1000000|kind=MEASUREMENT_WITH_RESIDUAL_BOUNDS|ordering=UNVERIFIED|ms={}|json=0",
        a[2],
        &e1 - &e0,
        t0.elapsed().as_millis()
    );
}

fn project_out(v: &mut [i128], g: &[i128]) {
    // unweighted (all weights 1); coefficient in fixed point 2^30 via big integers
    let gv: BigInt = v.iter().zip(g.iter()).map(|(x, y)| BigInt::from(*x) * BigInt::from(*y)).sum();
    let gg: BigInt = g.iter().map(|y| BigInt::from(*y) * BigInt::from(*y)).sum();
    let c = i128::try_from((gv << 30u32) / gg.max(BigInt::from(1))).unwrap_or(0);
    for (x, y) in v.iter_mut().zip(g.iter()) {
        *x -= (c * y) >> 30;
    }
}

fn filter(op: &Op, n: usize, beta: i128, guess: i128, steps: usize, defl: Option<&[i128]>) -> Vec<i128> {
    // start: positive smooth vector plus a small deterministic perturbation
    // (for the deflated run the ground component is projected out at once)
    let mut v: Vec<i128> = (0..n).map(|i| (1i128 << 38) + i128::from((i as u64).wrapping_mul(2654435761) % 1_000_003) * 1024).collect();
    if let Some(g) = defl {
        project_out(&mut v, g);
    }
    let mut av = vec![0i128; n];
    let guess_e = guess / 1_000_000;
    let mut done = 0usize;
    let mut block = 50usize;
    while done < steps {
        block = block.min(steps - done).max(2);
        let mut u = v.clone();
        let b = top_bits(&u);
        if b > 36 {
            shr(&mut u, b - 36);
        }
        op.apply(&u, &mut av);
        let num: i128 = u.iter().zip(av.iter()).map(|(x, y)| x * y).sum();
        let den: i128 = u.iter().map(|x| x * x).sum();
        let alpha = (num / den.max(1) + guess_e.max(1) / 2 + 1).min(beta - 1);
        let width = beta - alpha;
        let mut y0 = v.clone();
        op.apply(&y0, &mut av);
        let mut y1: Vec<i128> = y0.iter().zip(av.iter()).map(|(y, x)| ((alpha + beta) * y - 2 * x) / width).collect();
        for _ in 1..block {
            if let Some(g) = defl {
                project_out(&mut y1, g);
                project_out(&mut y0, g);
            }
            op.apply(&y1, &mut av);
            let y2: Vec<i128> = (0..n).map(|i| (2 * ((alpha + beta) * y1[i] - 2 * av[i]) - width * y0[i]) / width).collect();
            y0 = y1;
            y1 = y2;
            let b = top_bits(&y1);
            if b > VB {
                shr(&mut y0, b - VB);
                shr(&mut y1, b - VB);
            }
        }
        v = y1;
        if let Some(g) = defl {
            project_out(&mut v, g);
        }
        done += block;
        block *= 2;
    }
    v
}

/// exact E (micro), residual radius (micro, rounded up), <r^2> = <R>/9 (micro)
fn verify(op: &Op, v0: &[i128], robs: &[i64]) -> (BigInt, BigInt, BigInt) {
    let n = v0.len();
    let mut v = v0.to_vec();
    let mut av = vec![0i128; n];
    // exact verification
    let b = top_bits(&v);
    if b > 36 {
        shr(&mut v, b - 36);
    }
    op.apply(&v, &mut av);
    let bv: Vec<BigInt> = v.iter().map(|&x| BigInt::from(x)).collect();
    let ba: Vec<BigInt> = av.iter().map(|&x| BigInt::from(x)).collect();
    let mut nm = BigInt::from(0);
    let mut dn = BigInt::from(0);
    let mut rr_obs = BigInt::from(0);
    for i in 0..n {
        nm += &bv[i] * &ba[i];
        dn += &bv[i] * &bv[i];
        rr_obs += &bv[i] * &bv[i] * BigInt::from(robs[i]);
    }
    let mut res = BigInt::from(0);
    for i in 0..n {
        let r = &dn * &ba[i] - &nm * &bv[i];
        res += &r * &r;
    }
    // eps^2 = res / dn^3 ; eps_micro = ceil(1e6 sqrt)
    let d3 = &dn * &dn * &dn;
    let q = (res * BigInt::from(1_000_000_000_000u64) + &d3 - BigInt::from(1)) / &d3;
    let mut e = q.sqrt();
    if &e * &e < q {
        e += 1;
    }
    let e0 = micro(&nm, &dn);
    // <r^2> = <R>/9 in micro lattice units^2
    let r2 = &rr_obs * BigInt::from(1_000_000) / (&dn * BigInt::from(9));
    (e0, e, r2)
}
