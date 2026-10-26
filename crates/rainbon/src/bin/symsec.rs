//! Symmetry-sector solver for the rainbon (three-preon cloud): compress first,
//! then filter. Integer-only, Rust 1.81.
//!
//! Symmetries of H (equal-mass preons, symmetric potential, cubic periodic
//! lattice): the 48 cube operations (axis permutations x sign flips) acting on
//! both relative vectors, and the 6 preon permutations (positions (0, d2, d3)
//! permuted, then re-referenced to the new preon 1). 288 operations in all.
//!
//! sector A  : trivial character of all 288 operations. Contains the ground
//!             state (the electron).
//! sector Pz : subgroup K = {cube operations fixing the z axis up to sign
//!             (16)} x {preon identity, swap 2<->3} (32 elements), character
//!             chi = sign of the z flip. Contains the P-wave level (lambda mode
//!             polarised along z), part of the first excited multiplet.
//! Orbit basis: one representative per orbit (smallest index); a vector
//! constant up to chi on each orbit. Orbits whose stabiliser contains an
//! element with chi = -1 carry no amplitude and are dropped. The sector
//! operator A(i, c) = sum over hops i -> j of H(i,j) chi(g_j), g_j mapping j to
//! its representative c, is W-self-adjoint with w = orbit size (checked
//! exactly). Lowest state by the integer Chebyshev filter (W inner products),
//! exact residual radius, and the charge radius <r^2> = <R>/9 (R summed over
//! each orbit, which equals the full-lattice expectation).
//!
//! usage: symsec L harmonic|linear A|Pz JE JB STEPS GUESS_MICRO      json=0
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic)]

use num_bigint::BigInt;
use rainbon::{hops, potential, Lattice, Potential};
use std::collections::BTreeMap;

const VB: u32 = 40;

/// one group element: axis permutation, axis signs, preon permutation, character
#[derive(Clone, Copy)]
struct Elem {
    perm: [usize; 3],
    sign: [i64; 3],
    pre: [usize; 3],
    chi: i64,
}

fn group(sector: &str) -> Vec<Elem> {
    let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let pres_all = perms;
    let mut out = Vec::new();
    for p in perms {
        for s in 0..8u32 {
            let sign = [if s & 1 == 1 { -1 } else { 1 }, if s & 2 == 2 { -1 } else { 1 }, if s & 4 == 4 { -1 } else { 1 }];
            for pre in pres_all {
                if sector == "A" {
                    out.push(Elem { perm: p, sign, pre, chi: 1 });
                } else {
                    // Pz: cube ops with z -> +-z (perm fixes index 2), preon ops {id, (2 3)}
                    if p[2] != 2 || !(pre == [0, 1, 2] || pre == [0, 2, 1]) {
                        continue;
                    }
                    out.push(Elem { perm: p, sign, pre, chi: sign[2] });
                }
            }
        }
    }
    out
}

fn act(lat: &Lattice, e: &Elem, v: &[i64; 6]) -> [i64; 6] {
    // positions of the three preons: x1 = 0, x2 = d2, x3 = d3
    let pos = [[0i64; 3], [v[0], v[1], v[2]], [v[3], v[4], v[5]]];
    // cube operation on each position: new component k = sign[k] * old component perm[k]
    let rot = |x: &[i64; 3]| [e.sign[0] * x[e.perm[0]], e.sign[1] * x[e.perm[1]], e.sign[2] * x[e.perm[2]]];
    let r = [rot(&pos[0]), rot(&pos[1]), rot(&pos[2])];
    // preon permutation: new preon k is old preon pre[k]; re-reference to new preon 0
    let p = [r[e.pre[0]], r[e.pre[1]], r[e.pre[2]]];
    let mut out = [0i64; 6];
    for k in 0..3 {
        out[k] = (p[1][k] - p[0][k]).rem_euclid(lat.l);
        out[3 + k] = (p[2][k] - p[0][k]).rem_euclid(lat.l);
    }
    out
}

fn top_bits(v: &[i128]) -> u32 {
    128 - v.iter().map(|x| x.unsigned_abs()).max().unwrap_or(0).leading_zeros()
}
fn shr(v: &mut [i128], s: u32) {
    for x in v.iter_mut() {
        *x >>= s;
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 8 {
        eprintln!("usage: symsec L harmonic|linear A|Pz JE JB STEPS GUESS_MICRO");
        std::process::exit(1);
    }
    let t0 = std::time::Instant::now();
    let l: i64 = a[1].parse().unwrap_or(0);
    let kind = if a[2] == "linear" { Potential::Linear } else { Potential::Harmonic };
    let sector = a[3].as_str();
    let je: i64 = a[4].parse().unwrap_or(1);
    let jb: i64 = a[5].parse().unwrap_or(0);
    let steps: usize = a[6].parse().unwrap_or(2000);
    let guess: i128 = a[7].parse().unwrap_or(1_000_000);
    if !(3..=16).contains(&l) || !(sector == "A" || sector == "Pz") {
        println!("SYMSEC|error=bad_arguments|json=0");
        std::process::exit(1);
    }
    let lat = Lattice { l };
    let n = lat.dim();
    let g = group(sector);
    // representative and sign of every state; u32::MAX = dropped (chi-odd stabiliser)
    let mut rep = vec![u32::MAX; n];
    let mut sgn = vec![0i8; n];
    let mut reps: Vec<usize> = Vec::new();
    let mut wts: Vec<i64> = Vec::new();
    // sum of the radius observable over each orbit: R is not exactly invariant
    // under the symmetries when a separation sits at L/2 (minimum image), so the
    // orbit sum reproduces the full-lattice expectation exactly
    let mut rsum: Vec<i64> = Vec::new();
    let mut idx_of: BTreeMap<usize, u32> = BTreeMap::new();
    for i in 0..n {
        if rep[i] != u32::MAX || sgn[i] != 0 {
            continue;
        }
        let v = lat.decode(i);
        // orbit with characters: image -> chi; a conflicting chi means dropped
        let mut img: BTreeMap<usize, i64> = BTreeMap::new();
        let mut dropped = false;
        for e in &g {
            let j = lat.encode(&act(&lat, e, &v));
            match img.get(&j) {
                Some(&c) if c != e.chi => dropped = true,
                Some(_) => {}
                None => {
                    img.insert(j, e.chi);
                }
            }
        }
        let first = *img.keys().next().unwrap_or(&i);
        if dropped {
            for &j in img.keys() {
                sgn[j] = 2; // marks visited + dropped
            }
            continue;
        }
        // the representative is `first`; chi(g) maps j to first: v(j) = chi_j * v(first) / chi_first
        let chi_first = img[&first];
        let r = u32::try_from(reps.len()).unwrap_or(u32::MAX);
        for (&j, &c) in &img {
            rep[j] = r;
            sgn[j] = i8::try_from(c * chi_first).unwrap_or(0);
        }
        idx_of.insert(first, r);
        reps.push(first);
        wts.push(i64::try_from(img.len()).unwrap_or(0));
        rsum.push(img.keys().map(|&j| lat.radius_obs(&lat.decode(j))).sum());
    }
    let m = reps.len();
    // sector operator rows
    let mut rows: Vec<Vec<(u32, i64)>> = Vec::with_capacity(m);
    for &i in &reps {
        let v = lat.decode(i);
        let mut row: BTreeMap<u32, i64> = BTreeMap::new();
        *row.entry(rep[i]).or_insert(0) += 18 * je + jb * potential(&lat, &v, kind);
        for j in hops(&lat, &v) {
            if rep[j] == u32::MAX {
                continue;
            }
            *row.entry(rep[j]).or_insert(0) -= je * i64::from(sgn[j]);
        }
        rows.push(row.into_iter().filter(|&(_, c)| c != 0).collect());
    }
    // exact W-self-adjointness
    let mut sa = true;
    for (i, row) in rows.iter().enumerate() {
        for &(c, z) in row {
            let back = rows[c as usize].iter().find(|&&(k, _)| k as usize == i).map_or(0, |&(_, b)| b);
            if wts[i] * z != wts[c as usize] * back {
                sa = false;
            }
        }
    }
    let w: Vec<i128> = wts.iter().map(|&x| i128::from(x)).collect();
    let apply = |v: &[i128], out: &mut [i128]| {
        for (o, row) in out.iter_mut().zip(rows.iter()) {
            let mut s = 0i128;
            for &(c, z) in row {
                s += i128::from(z) * v[c as usize];
            }
            *o = s;
        }
    };
    let beta: i128 = rows.iter().map(|r| r.iter().map(|&(_, c)| i128::from(c.abs())).sum::<i128>()).max().unwrap_or(0) + 1;
    let mut v: Vec<i128> = if sector == "A" {
        vec![1i128 << 38; m]
    } else {
        let mut s = 777u64;
        (0..m)
            .map(|_| {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                i128::from(s >> 25) - (1i128 << 38)
            })
            .collect()
    };
    let mut av = vec![0i128; m];
    let guess_e = guess / 1_000_000;
    let mut done = 0usize;
    let mut block = 50usize;
    while done < steps {
        block = block.min(steps - done).max(2);
        let mut u = v.clone();
        let b = top_bits(&u);
        if b > 34 {
            shr(&mut u, b - 34);
        }
        apply(&u, &mut av);
        let num: i128 = (0..m).map(|i| w[i] * u[i] * av[i]).sum();
        let den: i128 = (0..m).map(|i| w[i] * u[i] * u[i]).sum();
        let alpha = (num / den.max(1) + guess_e.max(1) / 2 + 1).min(beta - 1);
        let width = beta - alpha;
        let mut y0 = v.clone();
        apply(&y0, &mut av);
        let mut y1: Vec<i128> = y0.iter().zip(av.iter()).map(|(y, x)| ((alpha + beta) * y - 2 * x) / width).collect();
        for _ in 1..block {
            apply(&y1, &mut av);
            let y2: Vec<i128> = (0..m).map(|i| (2 * ((alpha + beta) * y1[i] - 2 * av[i]) - width * y0[i]) / width).collect();
            y0 = y1;
            y1 = y2;
            let b = top_bits(&y1);
            if b > VB {
                shr(&mut y0, b - VB);
                shr(&mut y1, b - VB);
            }
        }
        v = y1;
        done += block;
        block *= 2;
    }
    // exact verification (W inner products)
    let b = top_bits(&v);
    if b > 34 {
        shr(&mut v, b - 34);
    }
    apply(&v, &mut av);
    let mut nm = BigInt::from(0);
    let mut dn = BigInt::from(0);
    let mut ro = BigInt::from(0);
    for i in 0..m {
        let wi = BigInt::from(w[i]);
        let vi = BigInt::from(v[i]);
        nm += &wi * &vi * BigInt::from(av[i]);
        dn += &wi * &vi * &vi;
        ro += &vi * &vi * BigInt::from(rsum[i]);
    }
    let mut res = BigInt::from(0);
    for i in 0..m {
        let r = &dn * BigInt::from(av[i]) - &nm * BigInt::from(v[i]);
        res += BigInt::from(w[i]) * &r * &r;
    }
    let d3 = &dn * &dn * &dn;
    let q = (res * BigInt::from(1_000_000_000_000u64) + &d3 - BigInt::from(1)) / &d3;
    let mut eps = q.sqrt();
    if &eps * &eps < q {
        eps += 1;
    }
    let e = &nm * BigInt::from(1_000_000) / &dn;
    let r2 = &ro * BigInt::from(1_000_000) / (&dn * BigInt::from(9));
    println!(
        "SYMSEC|L={l}|potential={}|sector={sector}|group={}|full_states={n}|sector_states={m}|w_self_adjoint={}|JE={je}|JB={jb}|steps={steps}|E_lowest={e}/1000000|eps_up={eps}/1000000|r2={r2}/1000000|kind=MEASUREMENT_WITH_RESIDUAL_BOUND|ms={}|json=0",
        a[2],
        g.len(),
        if sa { "PASS" } else { "FAIL" },
        t0.elapsed().as_millis()
    );
    let _ = idx_of;
}
