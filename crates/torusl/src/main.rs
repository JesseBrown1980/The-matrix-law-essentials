//! Gauge-orbit Hamiltonian of the periodic Lx x Ly S3 lattice (Kogut-Susskind
//! form, finite group), generalising native/src/torus2.rs (Stage 2C) from the
//! 2x2 torus. Integer-only (HEB-C rule), Rust 1.81, no external crates.
//!
//! Sites (x,y), index x*Ly+y. Link a = x-link (x,y)->(x+1,y); link N+a =
//! y-link (x,y)->(x,y+1). Spanning tree ("comb"): y-links on column x=0 for
//! y < Ly-1, and x-links (x,y) for x < Lx-1, every y. For 2x2 this is exactly
//! the Stage 2C tree, so FREE = [2,3,5,6,7] and the encoding is identical.
//! The k = N+1 free links are gauge-fixed loop variables; the residual gauge
//! freedom is global conjugation, and states are its orbits (weights = orbit
//! sizes 1,2,3,6). Electric term per link: sum over the 3 transpositions t of
//! (1 - L_t); magnetic term per plaquette: 2 - chi_std(U_p) in {0,2,3}.
//!
//! usage: torusl LX LY OUT_MATRIX
//! Writes the Sturm/evolver exchange format (n / weights / magnetic /
//! "i j z" for j >= i) and prints one TORUS line, json=0.
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use std::collections::{BTreeMap, VecDeque};
use std::io::Write;

const PERM: [[u8; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
const TRANS: [u8; 3] = [1, 2, 5];

fn mul_raw(a: u8, b: u8) -> u8 {
    let p = PERM[a as usize];
    let q = PERM[b as usize];
    let t = [p[q[0] as usize], p[q[1] as usize], p[q[2] as usize]];
    PERM.iter().position(|x| *x == t).map_or(255, |v| v as u8)
}

struct G {
    mul: [[u8; 6]; 6],
    inv: [u8; 6],
}

impl G {
    fn new() -> G {
        let mut mul = [[0u8; 6]; 6];
        for a in 0..6u8 {
            for b in 0..6u8 {
                mul[a as usize][b as usize] = mul_raw(a, b);
            }
        }
        let mut inv = [0u8; 6];
        for a in 0..6usize {
            inv[a] = (0..6u8).find(|&b| mul[a][b as usize] == 0).unwrap_or(255);
        }
        G { mul, inv }
    }
    fn m(&self, a: u8, b: u8) -> u8 {
        self.mul[a as usize][b as usize]
    }
    fn i(&self, a: u8) -> u8 {
        self.inv[a as usize]
    }
}

fn character(a: u8) -> i64 {
    if a == 0 {
        2
    } else if a == 3 || a == 4 {
        -1
    } else {
        0
    }
}

struct Lat {
    lx: usize,
    ly: usize,
    n: usize,
    edges: Vec<(usize, usize)>,
    tree: Vec<usize>, // in the order h is propagated
    free: Vec<usize>,
}

impl Lat {
    fn new(lx: usize, ly: usize) -> Lat {
        let n = lx * ly;
        let site = |x: usize, y: usize| (x % lx) * ly + (y % ly);
        let mut edges = vec![(0, 0); 2 * n];
        for x in 0..lx {
            for y in 0..ly {
                let a = site(x, y);
                edges[a] = (a, site(x + 1, y));
                edges[n + a] = (a, site(x, y + 1));
            }
        }
        let mut tree = Vec::new();
        for y in 0..ly - 1 {
            tree.push(n + site(0, y));
        }
        for x in 0..lx - 1 {
            for y in 0..ly {
                tree.push(site(x, y));
            }
        }
        let free: Vec<usize> = (0..2 * n).filter(|e| !tree.contains(e)).collect();
        Lat { lx, ly, n, edges, tree, free }
    }
    fn site(&self, x: usize, y: usize) -> usize {
        (x % self.lx) * self.ly + (y % self.ly)
    }
    fn encode(&self, q: &[u8]) -> usize {
        let mut id = 0;
        let mut p = 1;
        for &v in q {
            id += usize::from(v) * p;
            p *= 6;
        }
        id
    }
    fn decode(&self, mut id: usize) -> Vec<u8> {
        (0..self.free.len())
            .map(|_| {
                let v = (id % 6) as u8;
                id /= 6;
                v
            })
            .collect()
    }
    fn canonical(&self, g: &G, q: &[u8]) -> usize {
        (0..6u8)
            .map(|c| self.encode(&q.iter().map(|&x| g.m(g.m(c, x), g.i(c))).collect::<Vec<_>>()))
            .min()
            .unwrap_or(usize::MAX)
    }
    fn representative(&self, q: &[u8]) -> Vec<u8> {
        let mut links = vec![0u8; 2 * self.n];
        for (i, &e) in self.free.iter().enumerate() {
            links[e] = q[i];
        }
        links
    }
    fn gauge_fix(&self, g: &G, links: &[u8]) -> Vec<u8> {
        let mut h = vec![0u8; self.n];
        let mut set = vec![false; self.n];
        set[0] = true;
        for &e in &self.tree {
            let (t, hd) = self.edges[e];
            debug_assert!(set[t] && !set[hd]);
            h[hd] = g.m(h[t], links[e]);
            set[hd] = true;
        }
        self.free
            .iter()
            .map(|&e| {
                let (t, hd) = self.edges[e];
                g.m(g.m(h[t], links[e]), g.i(h[hd]))
            })
            .collect()
    }
    fn magnetic(&self, g: &G, s: &[u8]) -> i64 {
        let mut act = 0;
        for x in 0..self.lx {
            for y in 0..self.ly {
                let a = self.site(x, y);
                let right = self.site(x + 1, y);
                let up = self.site(x, y + 1);
                let p = g.m(g.m(g.m(s[a], s[self.n + right]), g.i(s[up])), g.i(s[self.n + a]));
                act += 2 - character(p);
            }
        }
        act
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: torusl LX LY OUT_MATRIX");
        std::process::exit(1);
    }
    let lx: usize = a[1].parse().unwrap_or(0);
    let ly: usize = a[2].parse().unwrap_or(0);
    if !(2..=4).contains(&lx) || !(2..=4).contains(&ly) || lx * ly > 9 {
        println!("TORUS|error=size_out_of_range|json=0");
        std::process::exit(1);
    }
    let t0 = std::time::Instant::now();
    let g = G::new();
    let lat = Lat::new(lx, ly);
    let k = lat.free.len();
    let tree_count = 6usize.pow(u32::try_from(k).unwrap_or(0));
    // weights by canonical key
    let mut w_by_key = vec![0u8; tree_count];
    for id in 0..tree_count {
        let c = lat.canonical(&g, &lat.decode(id));
        w_by_key[c] += 1;
    }
    let keys: Vec<usize> = (0..tree_count).filter(|&i| w_by_key[i] != 0).collect();
    let dim = keys.len();
    // Burnside: (6^k + 3*2^k + 2*3^k)/6
    let burn = (tree_count + 3 * (1usize << k) + 2 * 3usize.pow(u32::try_from(k).unwrap_or(0))) / 6;
    let weights: Vec<i64> = keys.iter().map(|&key| i64::from(w_by_key[key])).collect();
    let mut idx = vec![u32::MAX; tree_count];
    for (i, &key) in keys.iter().enumerate() {
        idx[key] = u32::try_from(i).unwrap_or(u32::MAX);
    }
    drop(w_by_key);
    let nlinks = 2 * lat.n;
    let diag = 3 * i64::try_from(nlinks).unwrap_or(0);
    let mut rows: Vec<Vec<(usize, i64)>> = Vec::with_capacity(dim);
    let mut mag = Vec::with_capacity(dim);
    for &key in &keys {
        let q = lat.decode(key);
        let raw = lat.representative(&q);
        mag.push(lat.magnetic(&g, &raw));
        let mut row = BTreeMap::<usize, i64>::new();
        *row.entry(idx[key] as usize).or_default() += diag;
        for e in 0..nlinks {
            for &t in &TRANS {
                let mut next = raw.clone();
                next[e] = g.m(t, next[e]);
                let j = idx[lat.canonical(&g, &lat.gauge_fix(&g, &next))];
                if j == u32::MAX {
                    println!("TORUS|error=electric_escaped|json=0");
                    std::process::exit(1);
                }
                *row.entry(j as usize).or_default() -= 1;
            }
        }
        rows.push(row.into_iter().filter(|&(_, v)| v != 0).collect());
    }
    // integer consistency checks
    let mut errors = Vec::new();
    if dim != burn {
        errors.push("burnside_mismatch");
    }
    if weights.iter().map(|&w| usize::try_from(w).unwrap_or(0)).sum::<usize>() != tree_count {
        errors.push("weight_sum");
    }
    if weights.iter().any(|w| ![1, 2, 3, 6].contains(w)) {
        errors.push("weight_value");
    }
    let maxmag = 3 * i64::try_from(lat.n).unwrap_or(0);
    let mut rowsum_bad = false;
    let mut sym_bad = false;
    let mut offdiag_bad = false;
    for (i, row) in rows.iter().enumerate() {
        if row.iter().map(|&(_, v)| v).sum::<i64>() != 0 {
            rowsum_bad = true;
        }
        if !(0..=maxmag).contains(&mag[i]) {
            errors.push("magnetic_range");
        }
        for &(j, v) in row {
            if j != i && v >= 0 {
                offdiag_bad = true;
            }
            let back = rows[j].iter().find(|&&(c, _)| c == i).map_or(0, |&(_, b)| b);
            if weights[i] * v != weights[j] * back {
                sym_bad = true;
            }
        }
    }
    if rowsum_bad {
        errors.push("row_sum");
    }
    if sym_bad {
        errors.push("weighted_self_adjointness");
    }
    if offdiag_bad {
        errors.push("offdiag_sign");
    }
    // connectivity
    let mut seen = vec![false; dim];
    let mut qd = VecDeque::from([0usize]);
    seen[0] = true;
    while let Some(i) = qd.pop_front() {
        for &(j, _) in &rows[i] {
            if !seen[j] {
                seen[j] = true;
                qd.push_back(j);
            }
        }
    }
    if !seen.iter().all(|&s| s) {
        errors.push("disconnected");
    }
    let flat = mag.iter().filter(|&&m| m == 0).count();
    let nnz: usize = rows.iter().map(Vec::len).sum();
    // write exchange format
    let file = std::fs::File::create(&a[3]).unwrap_or_else(|e| {
        println!("TORUS|error=write:{e}|json=0");
        std::process::exit(1)
    });
    let mut out = std::io::BufWriter::new(file);
    let join = |v: &[i64]| v.iter().map(i64::to_string).collect::<Vec<_>>().join(" ");
    let mut ok = writeln!(out, "{dim}").is_ok();
    ok &= writeln!(out, "{}", join(&weights)).is_ok();
    ok &= writeln!(out, "{}", join(&mag)).is_ok();
    for (i, row) in rows.iter().enumerate() {
        for &(j, v) in row {
            if j >= i {
                ok &= writeln!(out, "{i} {j} {v}").is_ok();
            }
        }
    }
    ok &= out.flush().is_ok();
    if !ok {
        errors.push("write_failed");
    }
    println!(
        "TORUS|L={lx}x{ly}|sites={}|links={nlinks}|free_loops={k}|tree_states={tree_count}|dim={dim}|burnside={burn}|nnz={nnz}|flat_orbits={flat}|electric_diag={diag}|magnetic_max={maxmag}|free={:?}|checks={}|ms={}|json=0",
        lat.n,
        lat.free,
        if errors.is_empty() { "PASS".to_string() } else { errors.join(",") },
        t0.elapsed().as_millis()
    );
    if !errors.is_empty() {
        std::process::exit(2);
    }
}
