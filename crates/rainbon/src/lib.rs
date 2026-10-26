//! Finite composite-electron toy model ("rainbon"): three preons on a periodic
//! L x L x L lattice, centre of mass removed (total momentum zero), bound by an
//! integer potential. Integer-only (HEB-C rule), Rust 1.81.
//!
//! State: relative positions d2 = x2 - x1, d3 = x3 - x1, each in Z_L^3, so the
//! space has L^6 states (weights all 1, the operator is symmetric).
//! H = JE * E + JB * M:
//!   E  kinetic lattice Laplacian of the three preons: diagonal 18 (3 preons x
//!      6 hops), -1 for each hop. Preon 2 moves d2, preon 3 moves d3, preon 1
//!      moves both d2 and d3 the opposite way. In units JE = hbar^2/(2 m a^2).
//!   M  binding potential on minimum-image pair vectors r12 = d2, r13 = d3,
//!      r23 = d3 - d2:
//!        harmonic : |r12|^2 + |r13|^2 + |r23|^2   (JB = k a^2 / 2)
//!        linear   : floor(10|r12|) + floor(10|r13|) + floor(10|r23|)
//!                   (JB = sigma a / 10; floor of 10x the Euclidean length by
//!                   integer square root)
//! Charge cloud: equal charges q/3 on equal-mass preons. With x1 = 0,
//!   sum_i |x_i - X|^2 = (2|d2|^2 + 2|d3|^2 - 2 d2.d3) / 3, and the charge
//!   radius squared is <r^2> = <sum_i |x_i - X|^2> / 3 = <R> / 9, R integer.

pub struct Lattice {
    pub l: i64,
}

impl Lattice {
    pub fn dim(&self) -> usize {
        usize::try_from(self.l.pow(6)).unwrap_or(0)
    }
    /// minimum-image component in (-L/2, L/2]
    pub fn mi(&self, x: i64) -> i64 {
        let y = x.rem_euclid(self.l);
        if 2 * y > self.l {
            y - self.l
        } else {
            y
        }
    }
    pub fn decode(&self, idx: usize) -> [i64; 6] {
        let mut v = [0i64; 6];
        let mut r = i64::try_from(idx).unwrap_or(0);
        for k in (0..6).rev() {
            v[k] = r % self.l;
            r /= self.l;
        }
        v
    }
    pub fn encode(&self, v: &[i64; 6]) -> usize {
        let mut r = 0i64;
        for &x in v {
            r = r * self.l + x.rem_euclid(self.l);
        }
        usize::try_from(r).unwrap_or(0)
    }
    /// minimum-image pair vectors r12, r13, r23
    pub fn pairs(&self, v: &[i64; 6]) -> [[i64; 3]; 3] {
        let d2 = [self.mi(v[0]), self.mi(v[1]), self.mi(v[2])];
        let d3 = [self.mi(v[3]), self.mi(v[4]), self.mi(v[5])];
        let d23 = [self.mi(v[3] - v[0]), self.mi(v[4] - v[1]), self.mi(v[5] - v[2])];
        [d2, d3, d23]
    }
    /// integer radius observable R (charge radius^2 = R / 9), using minimum-image
    /// positions x1 = 0, x2 = d2, x3 = d3
    pub fn radius_obs(&self, v: &[i64; 6]) -> i64 {
        let d2 = [self.mi(v[0]), self.mi(v[1]), self.mi(v[2])];
        let d3 = [self.mi(v[3]), self.mi(v[4]), self.mi(v[5])];
        (0..3).map(|k| 2 * d2[k] * d2[k] + 2 * d3[k] * d3[k] - 2 * d2[k] * d3[k]).sum()
    }
}

/// integer floor square root (Newton, integers only)
pub fn isqrt_int(n: i64) -> i64 {
    if n < 2 {
        return n.max(0);
    }
    let mut x = n;
    let mut y = (x + 1) / 2;
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Potential {
    Harmonic,
    Linear,
}

pub fn potential(lat: &Lattice, v: &[i64; 6], kind: Potential) -> i64 {
    let p = lat.pairs(v);
    match kind {
        Potential::Harmonic => p.iter().map(|r| r.iter().map(|c| c * c).sum::<i64>()).sum(),
        Potential::Linear => p.iter().map(|r| isqrt_int(100 * r.iter().map(|c| c * c).sum::<i64>())).sum(),
    }
}

/// neighbours under one preon hop (18 of them), with multiplicity
pub fn hops(lat: &Lattice, v: &[i64; 6]) -> Vec<usize> {
    let mut out = Vec::with_capacity(18);
    for k in 0..3 {
        for s in [-1i64, 1] {
            // preon 2: d2 += s e_k
            let mut a = *v;
            a[k] += s;
            out.push(lat.encode(&a));
            // preon 3: d3 += s e_k
            let mut b = *v;
            b[3 + k] += s;
            out.push(lat.encode(&b));
            // preon 1: x1 += s e_k  =>  d2 -= s e_k, d3 -= s e_k
            let mut c = *v;
            c[k] -= s;
            c[3 + k] -= s;
            out.push(lat.encode(&c));
        }
    }
    out
}
