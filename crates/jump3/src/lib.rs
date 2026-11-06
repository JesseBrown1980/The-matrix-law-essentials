//! Power-of-three jump, exact integers mod p = 1000081 (docs/LAW-POWER-OF-THREE-JUMP.md).
//!
//! Rust twin of python/matrix_tricks/jump3_ntt.py and transposition_memo.py.
//! The binaries `jump3` and `transmemo` print the same HBP rows, byte for
//! byte, as the Python reference (verify.sh compares them with cmp).
//!
//! A step that is the same at every site of a ring of D sites,
//! x'[i] = sum_s c_s x[(i - s) mod D] (mod p), reaches step n = 3^k by
//! k cubings of its D x D matrix (any fixed linear map), and, when D divides
//! p - 1, by the integer NTT: eigenvalues lambda_m = sum_s c_s w^(s m) raised
//! to n. mulmod fields count the modular products performed. SHA-256 is
//! implemented here (FIPS 180-4) so the crate has no dependencies.
#![deny(clippy::float_arithmetic)]

/// The kernel prime; p - 1 = 2^4 * 3^3 * 5 * 463.
pub const P: u64 = 1_000_081;
/// Largest k accepted: 3^80 < 2^127.
pub const K_MAX: u32 = 80;

/// Counter of modular products performed.
#[derive(Default)]
pub struct Ops {
    pub n: u64,
}

pub type Mat = Vec<Vec<u64>>;

pub fn factor(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            let mut e = 0;
            while n % d == 0 {
                n /= d;
                e += 1;
            }
            out.push((d, e));
        }
        d += 1;
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

pub fn factor_text(n: u64) -> String {
    factor(n)
        .iter()
        .map(|&(q, e)| {
            if e > 1 {
                format!("{q}^{e}")
            } else {
                format!("{q}")
            }
        })
        .collect::<Vec<_>>()
        .join("*")
}

/// Uncounted modular power for setup (roots, inverses of D).
pub fn modpow(b: u64, mut e: u128, m: u64) -> u64 {
    let mut r = 1 % m;
    let mut b = b % m;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % m;
        }
        e >>= 1;
        b = b * b % m;
    }
    r
}

/// Smallest primitive root mod P.
pub fn generator() -> u64 {
    let qs: Vec<u64> = factor(P - 1).iter().map(|&(q, _)| q).collect();
    let mut g = 2;
    while qs
        .iter()
        .any(|&q| modpow(g, u128::from((P - 1) / q), P) == 1)
    {
        g += 1;
    }
    g
}

/// Primitive d-th root of unity mod P, or None when d does not divide P - 1.
pub fn root_of_unity(d: usize) -> Option<u64> {
    let d64 = d as u64;
    if d < 2 || (P - 1) % d64 != 0 {
        return None;
    }
    let w = modpow(generator(), u128::from((P - 1) / d64), P);
    let ok = modpow(w, u128::from(d64), P) == 1
        && factor(d64)
            .iter()
            .all(|&(q, _)| modpow(w, u128::from(d64 / q), P) != 1);
    assert!(ok, "root of unity check failed");
    Some(w)
}

/// Right-to-left binary power mod P, every product counted.
pub fn powmod(b: u64, mut e: u128, ops: &mut Ops) -> u64 {
    let mut r = 1;
    let mut b = b % P;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % P;
            ops.n += 1;
        }
        e >>= 1;
        if e > 0 {
            b = b * b % P;
            ops.n += 1;
        }
    }
    r
}

/// lambda^n mod P for n >= 0; a zero eigenvalue gives 0 for every n >= 1.
pub fn lam_power(lam: u64, n: u128, ops: &mut Ops) -> u64 {
    if lam == 0 {
        return u64::from(n == 0);
    }
    powmod(lam, n % u128::from(P - 1), ops)
}

/// A step that is the same at every site of a ring of d sites.
pub struct Ring {
    pub d: usize,
    pub signed: Vec<(i64, i64)>,
    pub stencil: Vec<(usize, u64)>,
    pub w: Option<u64>,
    pub tw: Vec<u64>,
    pub tw_inv: Vec<u64>,
    pub inv_d: u64,
}

impl Ring {
    pub fn new(d: usize, stencil: &[(i64, i64)]) -> Result<Ring, String> {
        if d < 2 || stencil.is_empty() {
            return Err("ring of at least 2 sites and a nonempty stencil required".to_string());
        }
        let di = d as i64;
        let pi = P as i64;
        let reduced = stencil
            .iter()
            .map(|&(s, c)| (s.rem_euclid(di) as usize, c.rem_euclid(pi) as u64))
            .collect();
        let w = root_of_unity(d);
        let (mut tw, mut tw_inv, mut inv_d) = (Vec::new(), Vec::new(), 0);
        if let Some(w) = w {
            tw = vec![1u64; d];
            for e in 1..d {
                tw[e] = tw[e - 1] * w % P;
            }
            tw_inv = (0..d).map(|e| tw[(d - e) % d]).collect();
            inv_d = modpow(d as u64, u128::from(P - 2), P);
        }
        Ok(Ring {
            d,
            signed: stencil.to_vec(),
            stencil: reduced,
            w,
            tw,
            tw_inv,
            inv_d,
        })
    }

    pub fn stencil_text(&self) -> String {
        self.signed
            .iter()
            .map(|(s, c)| format!("{s}:{c}"))
            .collect::<Vec<_>>()
            .join(",")
    }
}

pub fn start_state(d: usize) -> Vec<u64> {
    (0..d as u64).map(|i| (7 * i * i + 3 * i + 1) % P).collect()
}

const K256: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

/// SHA-256 (FIPS 180-4).
pub fn sha256(data: &[u8]) -> [u8; 32] {
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let bit_len = (data.len() as u64).wrapping_mul(8);
    let mut msg = data.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());
    for chunk in msg.chunks(64) {
        let mut w = [0u32; 64];
        for (i, word) in chunk.chunks(4).enumerate() {
            w[i] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let mut a = h;
        for (k, wi) in K256.iter().zip(w.iter()) {
            let s1 = a[4].rotate_right(6) ^ a[4].rotate_right(11) ^ a[4].rotate_right(25);
            let ch = (a[4] & a[5]) ^ (!a[4] & a[6]);
            let t1 = a[7]
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(*k)
                .wrapping_add(*wi);
            let s0 = a[0].rotate_right(2) ^ a[0].rotate_right(13) ^ a[0].rotate_right(22);
            let maj = (a[0] & a[1]) ^ (a[0] & a[2]) ^ (a[1] & a[2]);
            let t2 = s0.wrapping_add(maj);
            a = [
                t1.wrapping_add(t2),
                a[0],
                a[1],
                a[2],
                a[3].wrapping_add(t1),
                a[4],
                a[5],
                a[6],
            ];
        }
        for (x, y) in h.iter_mut().zip(a.iter()) {
            *x = x.wrapping_add(*y);
        }
    }
    let mut out = [0u8; 32];
    for (o, x) in out.chunks_mut(4).zip(h.iter()) {
        o.copy_from_slice(&x.to_be_bytes());
    }
    out
}

pub fn hex(bytes: &[u8]) -> String {
    use std::fmt::Write;
    bytes
        .iter()
        .fold(String::with_capacity(2 * bytes.len()), |mut s, b| {
            let _ = write!(s, "{b:02x}");
            s
        })
}

pub fn sha256_hex(data: &[u8]) -> String {
    hex(&sha256(data))
}

/// SHA-256 of the components as 8-byte big-endian integers (the seat's encoding).
pub fn state_sha(x: &[u64]) -> String {
    let bytes: Vec<u8> = x.iter().flat_map(|v| v.to_be_bytes()).collect();
    sha256_hex(&bytes)
}

pub fn step(ring: &Ring, x: &[u64], ops: &mut Ops) -> Vec<u64> {
    let d = ring.d;
    let mut acc = vec![0u64; d];
    for &(s, c) in &ring.stencil {
        for (i, a) in acc.iter_mut().enumerate() {
            *a += c * x[(i + d - s) % d];
        }
    }
    ops.n += (d * ring.stencil.len()) as u64;
    acc.iter().map(|a| a % P).collect()
}

pub fn circulant(ring: &Ring) -> Mat {
    let d = ring.d;
    let mut m = vec![vec![0u64; d]; d];
    for (i, row) in m.iter_mut().enumerate() {
        for &(s, c) in &ring.stencil {
            let j = (i + d - s) % d;
            row[j] = (row[j] + c) % P;
        }
    }
    m
}

pub fn matmul(a: &Mat, b: &Mat, ops: &mut Ops) -> Mat {
    let d = a.len();
    ops.n += (d * d * d) as u64;
    a.iter()
        .map(|row| {
            (0..d)
                .map(|j| {
                    row.iter()
                        .zip(b.iter())
                        .map(|(u, brow)| u * brow[j])
                        .sum::<u64>()
                        % P
                })
                .collect()
        })
        .collect()
}

pub fn matvec(a: &Mat, x: &[u64], ops: &mut Ops) -> Vec<u64> {
    ops.n += (a.len() * x.len()) as u64;
    a.iter()
        .map(|row| row.iter().zip(x.iter()).map(|(u, v)| u * v).sum::<u64>() % P)
        .collect()
}

/// a^3 mod P: two matrix products.
pub fn cube(a: &Mat, ops: &mut Ops) -> Mat {
    matmul(&matmul(a, a, ops), a, ops)
}

/// (m^(3^r), products so far) for r = 0..kmax: one cache of cubings serves every k.
pub fn cube_ladder(m: Mat, kmax: u32, ops: &mut Ops) -> Vec<(Mat, u64)> {
    let mut out = vec![(m, ops.n)];
    for _ in 0..kmax {
        let next = cube(&out[out.len() - 1].0, ops);
        out.push((next, ops.n));
    }
    out
}

pub fn ntt(x: &[u64], tw: &[u64], ops: &mut Ops) -> Vec<u64> {
    let d = x.len();
    ops.n += (d * d) as u64;
    (0..d)
        .map(|m| {
            x.iter()
                .enumerate()
                .map(|(j, v)| v * tw[j * m % d])
                .sum::<u64>()
                % P
        })
        .collect()
}

pub fn eigenvalues(ring: &Ring, ops: &mut Ops) -> Vec<u64> {
    let d = ring.d;
    ops.n += (d * ring.stencil.len()) as u64;
    (0..d)
        .map(|m| {
            ring.stencil
                .iter()
                .map(|&(s, c)| c * ring.tw[s * m % d])
                .sum::<u64>()
                % P
        })
        .collect()
}

pub fn ntt_apply(ring: &Ring, x: &[u64], factors: &[u64], ops: &mut Ops) -> Vec<u64> {
    let x_hat = ntt(x, &ring.tw, ops);
    let y_hat: Vec<u64> = x_hat
        .iter()
        .zip(factors.iter())
        .map(|(a, f)| a * f % P)
        .collect();
    ops.n += ring.d as u64;
    let y = ntt(&y_hat, &ring.tw_inv, ops);
    ops.n += ring.d as u64;
    y.iter().map(|v| v * ring.inv_d % P).collect()
}

/// M^n x through the eigenvalues: D powers, two transforms.
pub fn ntt_jump(ring: &Ring, lam: &[u64], x: &[u64], n: u128, ops: &mut Ops) -> Vec<u64> {
    let factors: Vec<u64> = lam.iter().map(|&v| lam_power(v, n, ops)).collect();
    ntt_apply(ring, x, &factors, ops)
}

/// M^-1 x through the inverse eigenvalues; None when an eigenvalue is zero.
pub fn ntt_unstep(ring: &Ring, lam: &[u64], x: &[u64], ops: &mut Ops) -> Option<Vec<u64>> {
    if lam.iter().any(|&v| v == 0) {
        return None;
    }
    let factors: Vec<u64> = lam
        .iter()
        .map(|&v| powmod(v, u128::from(P - 2), ops))
        .collect();
    Some(ntt_apply(ring, x, &factors, ops))
}

/// Decimal digits of n * m for n < 2^127 and m < 2^32, without overflow.
pub fn mul_small_decimal(n: u128, m: u64) -> String {
    const BASE: u128 = 10_000_000_000_000_000_000; // 10^19
    let m = u128::from(m);
    let lo = (n % BASE) * m;
    let hi = (n / BASE) * m + lo / BASE;
    let lo = lo % BASE;
    if hi == 0 {
        format!("{lo}")
    } else {
        format!("{hi}{lo:019}")
    }
}

/// The rows joined, each ending in a newline.
pub fn rows_text(rows: &[String]) -> String {
    let mut s = String::new();
    for r in rows {
        s.push_str(r);
        s.push('\n');
    }
    s
}

pub fn mark(v: bool) -> &'static str {
    if v {
        "1"
    } else {
        "0"
    }
}

/// One case of the self-test: a ring, its stencil, the powers k and the largest k stepped directly.
pub struct Case {
    pub name: String,
    pub d: usize,
    pub stencil: Vec<(i64, i64)>,
    pub label: String,
    pub ks: Vec<u32>,
    pub direct_max: u32,
}

pub struct Flags {
    pub cube_ntt: bool,
    pub direct: bool,
    pub leaf_f: bool,
    pub leaf_b: bool,
    pub zero: bool,
    pub shift: bool,
    pub rows: u64,
    pub zero_rows: u64,
}

impl Default for Flags {
    fn default() -> Self {
        Flags {
            cube_ntt: true,
            direct: true,
            leaf_f: true,
            leaf_b: true,
            zero: true,
            shift: true,
            rows: 0,
            zero_rows: 0,
        }
    }
}

pub fn cases() -> Vec<Case> {
    let mk =
        |name: &str, d: usize, stencil: &[(i64, i64)], label: &str, ks: &[u32], direct_max: u32| {
            Case {
                name: name.to_string(),
                d,
                stencil: stencil.to_vec(),
                label: label.to_string(),
                ks: ks.to_vec(),
                direct_max,
            }
        };
    vec![
        mk(
            "signed_return",
            27,
            &[(0, 1), (1, 1), (2, -1)],
            "I+R-R^2",
            &[1, 2, 3, 6, 9, 10, 20, 40, 80],
            10,
        ),
        mk(
            "curved_sum",
            27,
            &[(0, 1), (1, 1), (2, 1)],
            "I+R+R^2",
            &[1, 2, 3, 9, 40],
            9,
        ),
        mk("shift", 27, &[(1, 1)], "R", &[1, 2, 3, 40], 3),
        mk(
            "ring45",
            45,
            &[(0, 2), (1, -1), (44, 3)],
            "2I-R+3R^-1",
            &[1, 2, 3, 6, 9, 20],
            9,
        ),
        mk(
            "ring81",
            81,
            &[(0, 1), (1, 1), (80, -1)],
            "I+R-R^-1",
            &[1, 2, 3],
            3,
        ),
    ]
}

pub fn head_row() -> String {
    format!(
        "JUMP3_HEAD|law=power_of_three_jump|p={P}|p_minus_1_factors={}|generator={}|\
         start=7i^2+3i+1_mod_p|state_hash=sha256_of_8_byte_big_endian_components|json=0",
        factor_text(P - 1),
        generator()
    )
}

/// Append the JUMP3_CASE, JUMP3 and JUMP3_ZERO_MODE rows of one case to out.
pub fn run_case(case: &Case, out: &mut Vec<String>, flags: &mut Flags) -> Result<(), String> {
    if case.ks.is_empty() || case.ks.iter().any(|&k| !(1..=K_MAX).contains(&k)) {
        return Err(format!("k must lie in 1..{K_MAX}"));
    }
    let ring = Ring::new(case.d, &case.stencil)?;
    let d = case.d;
    let x0 = start_state(d);
    let s_len = ring.stencil.len() as u64;
    let mut ops = Ops::default();
    let has_ntt = ring.w.is_some();
    let (lam, zeros, root, ntt_state, zm, zl, backward);
    if let Some(w) = ring.w {
        lam = eigenvalues(&ring, &mut ops);
        zeros = (0..d).filter(|&m| lam[m] == 0).collect::<Vec<_>>();
        root = w.to_string();
        ntt_state = "DEFINED";
        zm = zeros.len().to_string();
        zl = if zeros.is_empty() {
            "none".to_string()
        } else {
            zeros
                .iter()
                .map(|m| m.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };
        backward = if zeros.is_empty() {
            "DEFINED"
        } else {
            "REFUSED_ZERO_MODE"
        };
    } else {
        lam = Vec::new();
        zeros = Vec::new();
        root = "none".to_string();
        ntt_state = "REFUSED_RING_DOES_NOT_DIVIDE_P_MINUS_1";
        zm = "NOT_DEFINED".to_string();
        zl = "none".to_string();
        backward = "NOT_DEFINED";
    }
    out.push(format!(
        "JUMP3_CASE|case={}|ring={d}|step={}|stencil={}|root={root}|ntt={ntt_state}|eigen_mulmods={}|\
         zero_modes={zm}|zero_mode_list={zl}|backward={backward}|start_sha256={}|json=0",
        case.name,
        case.label,
        ring.stencil_text(),
        ops.n,
        state_sha(&x0)
    ));
    // direct stepping, one step at a time, to 3^k and 3^k + 1 for every k <= direct_max
    let mut targets: Vec<u128> = Vec::new();
    for &k in &case.ks {
        if k <= case.direct_max {
            targets.push(3u128.pow(k));
            targets.push(3u128.pow(k) + 1);
        }
    }
    let mut direct: std::collections::BTreeMap<u128, (Vec<u64>, u64)> =
        std::collections::BTreeMap::new();
    if let Some(&last) = targets.iter().max() {
        let mut dops = Ops::default();
        let mut x = x0.clone();
        for t in 1..=last {
            x = step(&ring, &x, &mut dops);
            if targets.contains(&t) {
                direct.insert(t, (x.clone(), dops.n));
            }
        }
    }
    let kmax = *case.ks.iter().max().unwrap_or(&1);
    let ladder = cube_ladder(circulant(&ring), kmax, &mut Ops::default());
    let mut scratch = Ops::default();
    for &k in &case.ks {
        let n = 3u128.pow(k);
        let (mk, lad) = &ladder[k as usize];
        let mut cops = Ops { n: *lad };
        let xc = matvec(mk, &x0, &mut cops);
        let (ntt_mul, c_eq) = if has_ntt {
            let mut nops = Ops::default();
            let xn = ntt_jump(&ring, &lam, &x0, n, &mut nops);
            flags.cube_ntt &= xc == xn;
            (nops.n.to_string(), mark(xc == xn).to_string())
        } else {
            ("NOT_DEFINED".to_string(), "NOT_DEFINED".to_string())
        };
        let (d_run, d_mul, d_eq) = match direct.get(&n) {
            Some((xd, cnt)) => {
                flags.direct &= *xd == xc;
                ("RUN", cnt.to_string(), mark(*xd == xc).to_string())
            }
            None => (
                "NOT_RUN",
                mul_small_decimal(n, d as u64 * s_len),
                "NOT_RUN".to_string(),
            ),
        };
        let fwd = step(&ring, &xc, &mut scratch);
        let leaf_f = if has_ntt {
            let ok = fwd == ntt_jump(&ring, &lam, &x0, n + 1, &mut scratch);
            flags.leaf_f &= ok;
            mark(ok)
        } else if let Some((xd, _)) = direct.get(&(n + 1)) {
            let ok = fwd == *xd;
            flags.leaf_f &= ok;
            mark(ok)
        } else {
            "NOT_DEFINED"
        };
        let leaf_b = if has_ntt && zeros.is_empty() {
            let b = ntt_unstep(&ring, &lam, &xc, &mut scratch).ok_or("zero eigenvalue")?;
            let x1 = step(&ring, &x0, &mut scratch);
            let ok = step(&ring, &b, &mut scratch) == xc
                && b == ntt_jump(&ring, &lam, &x0, n - 1, &mut scratch)
                && ntt_unstep(&ring, &lam, &x1, &mut scratch).as_ref() == Some(&x0);
            flags.leaf_b &= ok;
            mark(ok)
        } else if has_ntt {
            "REFUSED_ZERO_MODE"
        } else {
            "NOT_DEFINED"
        };
        let at_start = xc == x0;
        if case.name == "shift" {
            flags.shift &= at_start == (k >= 3);
        }
        flags.rows += 1;
        out.push(format!(
            "JUMP3|case={}|k={k}|n={n}|cubings={k}|cube_mulmods={}|ntt_mulmods={ntt_mul}|cube_equals_ntt={c_eq}|\
             direct={d_run}|direct_mulmods={d_mul}|direct_equals_jump={d_eq}|equals_start={}|leaf_forward={leaf_f}|\
             leaf_back={leaf_b}|sha256={}|json=0",
            case.name,
            cops.n,
            mark(at_start),
            state_sha(&xc)
        ));
    }
    if has_ntt && !zeros.is_empty() {
        let k = kmax;
        let n = 3u128.pow(k);
        let mk = &ladder[k as usize].0;
        let jump_k = matvec(mk, &x0, &mut scratch);
        for &m in &zeros {
            let v: Vec<u64> = (0..d).map(|j| ring.tw_inv[j * m % d]).collect(); // eigenvector w^(-j*m), eigenvalue 0
            let xb: Vec<u64> = x0.iter().zip(v.iter()).map(|(a, b)| (a + b) % P).collect();
            let image_b = matvec(mk, &xb, &mut scratch);
            let ntt_a = ntt_jump(&ring, &lam, &x0, n, &mut scratch);
            let checks = [
                xb != x0,
                step(&ring, &x0, &mut scratch) == step(&ring, &xb, &mut scratch),
                jump_k == image_b,
                ntt_a == ntt_jump(&ring, &lam, &xb, n, &mut scratch),
                image_b == ntt_a,
            ];
            for c in checks {
                flags.zero &= c;
            }
            flags.zero_rows += 1;
            out.push(format!(
                "JUMP3_ZERO_MODE|case={}|mode={m}|k={k}|preimages_distinct={}|images_equal_after_1_step={}|\
                 images_equal_after_3^{k}_cube={}|images_equal_after_3^{k}_ntt={}|image_equals_case_jump={}|\
                 preimage_b_sha256={}|image_sha256={}|json=0",
                case.name,
                mark(checks[0]),
                mark(checks[1]),
                mark(checks[2]),
                mark(checks[3]),
                mark(checks[4]),
                state_sha(&xb),
                state_sha(&image_b)
            ));
        }
    }
    Ok(())
}

pub fn verify_row(flags: &Flags, cases: usize) -> (String, bool) {
    let ok =
        flags.cube_ntt && flags.direct && flags.leaf_f && flags.leaf_b && flags.zero && flags.shift;
    (
        format!(
            "JUMP3_VERIFY|status={}|cases={cases}|jump_rows={}|cube_equals_ntt_all={}|direct_equals_jump_all={}|\
             leaf_forward_all={}|leaf_back_all={}|zero_mode_rows={}|zero_mode_all={}|shift_returns_at_27={}|json=0",
            if ok { "PASS" } else { "FAIL" },
            flags.rows,
            mark(flags.cube_ntt),
            mark(flags.direct),
            mark(flags.leaf_f),
            mark(flags.leaf_b),
            flags.zero_rows,
            mark(flags.zero),
            mark(flags.shift)
        ),
        ok,
    )
}

pub fn parse_stencil(text: &str) -> Result<Vec<(i64, i64)>, String> {
    text.split(',')
        .map(|part| {
            let (s, c) = part
                .split_once(':')
                .ok_or_else(|| format!("bad stencil entry {part}"))?;
            Ok((
                s.parse::<i64>().map_err(|e| e.to_string())?,
                c.parse::<i64>().map_err(|e| e.to_string())?,
            ))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_vectors() {
        assert_eq!(
            sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(
            sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
    }

    #[test]
    fn root_and_generator() {
        assert_eq!(factor_text(P - 1), "2^4*3^3*5*463");
        assert_eq!(root_of_unity(27), Some(951_846));
        assert_eq!(root_of_unity(81), None);
    }

    #[test]
    fn jump_equals_direct_and_leaves() {
        let ring = Ring::new(27, &[(0, 1), (1, 1), (2, -1)]).unwrap();
        let x0 = start_state(27);
        let mut ops = Ops::default();
        let lam = eigenvalues(&ring, &mut ops);
        let mut x = x0.clone();
        for _ in 0..243 {
            x = step(&ring, &x, &mut ops);
        }
        let ladder = cube_ladder(circulant(&ring), 5, &mut ops);
        assert_eq!(matvec(&ladder[5].0, &x0, &mut ops), x);
        assert_eq!(ntt_jump(&ring, &lam, &x0, 243, &mut ops), x);
        let back = ntt_unstep(&ring, &lam, &x, &mut ops).unwrap();
        assert_eq!(step(&ring, &back, &mut ops), x);
    }

    #[test]
    fn zero_mode_refuses_backward() {
        let ring = Ring::new(27, &[(0, 1), (1, 1), (2, 1)]).unwrap();
        let mut ops = Ops::default();
        let lam = eigenvalues(&ring, &mut ops);
        let zeros: Vec<usize> = (0..27).filter(|&m| lam[m] == 0).collect();
        assert_eq!(zeros, vec![9, 18]);
        assert!(ntt_unstep(&ring, &lam, &start_state(27), &mut ops).is_none());
    }

    #[test]
    fn decimal_product() {
        let n = 3u128.pow(80);
        assert_eq!(
            mul_small_decimal(n, 81),
            "11972515182562019788602740026717047105681"
        );
        assert_eq!(mul_small_decimal(27, 3), "81");
    }
}
