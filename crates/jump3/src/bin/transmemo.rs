//! transmemo: hash-keyed position memo on the 27-state cube, the path-count
//! recurrence C_{n+1}(s) = sum over p -> s of C_n(p), the digit-sum invariant,
//! and the counts after 3^k turns by k cubings (exact mod p).
//!
//! usage: transmemo --self-test [--timing]
//!
//! Prints the same HBP rows, byte for byte, as
//! python/matrix_tricks/transposition_memo.py (timing rows excepted). json=0.
#![deny(clippy::float_arithmetic)]

use jump3::{
    cube_ladder, mark, matvec, modpow, rows_text, sha256, sha256_hex, state_sha, Mat, Ops, P,
};
use std::collections::BTreeMap;
use std::time::Instant;

/// Position index a*9 + b*3 + c.
fn digits(i: usize) -> [u8; 3] {
    [(i / 9) as u8, (i / 3 % 3) as u8, (i % 3) as u8]
}

fn index(p: [u8; 3]) -> usize {
    usize::from(p[0]) * 9 + usize::from(p[1]) * 3 + usize::from(p[2])
}

fn name(i: usize) -> String {
    let p = digits(i);
    format!("{}{}{}", p[0], p[1], p[2])
}

fn turn(i: usize, d: usize) -> usize {
    let mut p = digits(i);
    p[d] = (p[d] + 1) % 3;
    index(p)
}

fn digit_sum(i: usize) -> usize {
    digits(i).iter().map(|&v| usize::from(v)).sum()
}

/// Every one of the 3^L move orders; the end position looked up by its hash.
fn enumerate_orders(l: u32) -> ([u64; 27], BTreeMap<String, u64>, u64, String) {
    let slot: BTreeMap<String, usize> = (0..27).map(|i| (sha256_hex(&digits(i)), i)).collect();
    let mut table: BTreeMap<String, u64> = BTreeMap::new();
    let mut turns = 0u64;
    for code in 0..3u64.pow(l) {
        let mut p = [0u8; 3];
        let mut c = code;
        for _ in 0..l {
            let d = (c % 3) as usize;
            c /= 3;
            p[d] = (p[d] + 1) % 3;
        }
        turns += u64::from(l);
        *table.entry(sha256_hex(&p)).or_insert(0) += 1;
    }
    let mut counts = [0u64; 27];
    for (key, c) in &table {
        counts[slot[key]] = *c;
    }
    let text = rows_text(
        &table
            .iter()
            .map(|(k, c)| format!("{k}={c}"))
            .collect::<Vec<_>>(),
    );
    let table_sha = sha256_hex(text.as_bytes());
    (counts, table, turns, table_sha)
}

/// C_L over all 27 positions; exact (u128, L <= 80), or mod `modulus` when it is nonzero.
fn recurrence(l: u32, modulus: u128) -> ([u128; 27], u64) {
    let mut c = [0u128; 27];
    c[0] = 1;
    let mut updates = 0u64;
    for _ in 0..l {
        let mut nxt = [0u128; 27];
        for (q, &cq) in c.iter().enumerate() {
            for d in 0..3 {
                nxt[turn(q, d)] += cq;
                updates += 1;
            }
        }
        if modulus > 0 {
            for v in nxt.iter_mut() {
                *v %= modulus;
            }
        }
        c = nxt;
    }
    (c, updates)
}

/// The same recurrence over the 9 positions of the current class only.
fn recurrence_invariant(l: u32) -> ([u128; 27], u64) {
    let mut c = [0u128; 27];
    c[0] = 1;
    let mut updates = 0u64;
    for n in 0..l as usize {
        let mut nxt = [0u128; 27];
        for (q, &cq) in c.iter().enumerate() {
            if digit_sum(q) % 3 != n % 3 {
                continue;
            }
            for d in 0..3 {
                nxt[turn(q, d)] += cq;
                updates += 1;
            }
        }
        c = nxt;
    }
    (c, updates)
}

/// A[s][q] = number of turns taking q to s.
fn adjacency() -> Mat {
    let mut a = vec![vec![0u64; 27]; 27];
    for q in 0..27 {
        for d in 0..3 {
            a[turn(q, d)][q] += 1;
        }
    }
    a
}

fn counts_text(counts: &[u128]) -> String {
    counts
        .iter()
        .enumerate()
        .filter(|(_, &v)| v != 0)
        .map(|(i, v)| format!("{}:{v}", name(i)))
        .collect::<Vec<_>>()
        .join(",")
}

fn exact_sha(counts: &[u128]) -> String {
    sha256_hex(
        counts
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",")
            .as_bytes(),
    )
}

const PERM_012: [[u8; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

fn self_test(timing: bool) -> bool {
    let mut out = vec![format!(
        "TRANSMEMO_HEAD|cube=3x3x3|positions=27|turn=add_1_mod_3_to_one_digit|start=000|\
         position_key=sha256_of_3_digit_bytes|invariant=digit_sum_mod_3_equals_turns_mod_3|p={P}|json=0"
    )];
    let mut ok = true;
    let mut rows_t = 0;
    for l in [6u32, 9, 12] {
        let (counts, table, turns, table_sha) = enumerate_orders(l);
        let counts: Vec<u128> = counts.iter().map(|&v| u128::from(v)).collect();
        let paths = 3u64.pow(l);
        let (rec, upd) = recurrence(l, 0);
        let (inv, upd_inv) = recurrence_invariant(l);
        let cls: Vec<usize> = (0..27)
            .filter(|&i| digit_sum(i) % 3 == l as usize % 3)
            .collect();
        let reached: Vec<usize> = (0..27).filter(|&i| counts[i] != 0).collect();
        let nz: Vec<u128> = counts.iter().copied().filter(|&v| v != 0).collect();
        let checks = [
            rec.to_vec() == counts,
            reached == cls,
            inv == rec,
            counts.iter().sum::<u128>() == u128::from(paths),
        ];
        ok &= checks.iter().all(|&c| c);
        rows_t += 1;
        let keys = table.len() as u64;
        out.push(format!(
            "TRANSPOSITION|L={l}|paths={paths}|enumeration_turns={turns}|distinct_positions={}|memo_keys={keys}|\
             memo_hits={}|reuse_avg={}|paths_min={}|paths_max={}|counts={}|table_sha256={table_sha}|\
             recurrence_updates={upd}|recurrence_equals_enumeration={}|invariant_class={}|\
             invariant_candidates={}|reachable_equals_invariant_class={}|invariant_updates={upd_inv}|\
             invariant_equals_recurrence={}|json=0",
            reached.len(),
            paths - keys,
            paths / keys,
            nz.iter().min().copied().unwrap_or(0),
            nz.iter().max().copied().unwrap_or(0),
            counts_text(&counts),
            mark(checks[0]),
            l % 3,
            cls.len(),
            mark(checks[1]),
            mark(checks[2])
        ));
    }
    let mut rows_p = 0;
    for l in [27u32, 54] {
        let (c, _) = recurrence(l, 0);
        let perm: Vec<u128> = PERM_012.iter().map(|&p| c[index(p)]).collect();
        let sum: u128 = c.iter().sum();
        let support = c.iter().filter(|&&v| v != 0).count();
        let checks = [
            sum == 3u128.pow(l),
            perm[0] == 3u128.pow(l - 2),
            perm.iter().all(|&v| v == perm[0]),
        ];
        ok &= checks.iter().all(|&v| v) && support == 9;
        rows_p += 1;
        out.push(format!(
            "PATHCOUNT|L={l}|route=recurrence_exact|sum={sum}|sum_equals_3^L={}|support={support}|count_000={}|\
             count_111={}|count_222={}|count_012={}|count_012_equals_3^(L-2)={}|perm_012_all_equal={}|\
             counts_sha256={}|json=0",
            mark(checks[0]),
            c[index([0, 0, 0])],
            c[index([1, 1, 1])],
            c[index([2, 2, 2])],
            perm[0],
            mark(checks[1]),
            mark(checks[2]),
            exact_sha(&c)
        ));
    }
    let ladder = cube_ladder(adjacency(), 40, &mut Ops::default());
    let mut e0 = vec![0u64; 27];
    e0[0] = 1;
    let mut rows_j = 0;
    for k in [3u32, 4, 40] {
        let l = 3u128.pow(k);
        let (mk, lad) = &ladder[k as usize];
        let mut ops = Ops { n: *lad };
        let c = matvec(mk, &e0, &mut ops);
        let eq = if k <= 4 {
            let (rec, _) = recurrence(l as u32, u128::from(P));
            let same = rec.iter().zip(c.iter()).all(|(a, b)| *a == u128::from(*b));
            ok &= same;
            mark(same)
        } else {
            "NOT_RUN"
        };
        let s = c.iter().sum::<u64>() % P;
        let want_sum = modpow(3, l % u128::from(P - 1), P);
        let want_012 = modpow(3, (l - 2) % u128::from(P - 1), P);
        let perm_ok = PERM_012.iter().all(|&p| c[index(p)] == want_012);
        ok &= s == want_sum && perm_ok;
        rows_j += 1;
        out.push(format!(
            "PATHCOUNT_JUMP|k={k}|L={l}|route=cube_mod_p|cubings={k}|cube_mulmods={}|equals_recurrence_mod_p={eq}|\
             support_mod_p={}|sum_mod_p={s}|sum_equals_3^L_mod_p={}|count_012_mod_p={}|\
             count_012_equals_3^(L-2)_mod_p={}|counts_sha256={}|json=0",
            ops.n,
            c.iter().filter(|&&v| v != 0).count(),
            mark(s == want_sum),
            c[index([0, 1, 2])],
            mark(perm_ok),
            state_sha(&c)
        ));
    }
    let hidden = [0usize, 2, 1, 1, 0, 2, 2, 1, 0, 1];
    let p = hidden.iter().fold(0, |p, &d| turn(p, d));
    let target = sha256(&digits(p));
    let mut tried_full = 0;
    let mut found = None;
    for i in 0..27 {
        tried_full += 1;
        if sha256(&digits(i)) == target {
            found = Some(i);
            break;
        }
    }
    let cls: Vec<usize> = (0..27)
        .filter(|&i| digit_sum(i) % 3 == hidden.len() % 3)
        .collect();
    let mut tried_inv = 0;
    let mut found_inv = None;
    for &i in &cls {
        tried_inv += 1;
        if sha256(&digits(i)) == target {
            found_inv = Some(i);
            break;
        }
    }
    let rec_ok = found == Some(p) && found_inv == Some(p);
    ok &= rec_ok;
    out.push(format!(
        "FINITE_RECOVER|moves_hidden={}|hidden_moves={}|position={}|bound_full=27|tried_full={tried_full}|\
         bound_invariant={}|tried_invariant={tried_inv}|recovered_equals_hidden={}|json=0",
        hidden.len(),
        hidden.iter().map(|d| d.to_string()).collect::<Vec<_>>().join(","),
        name(p),
        cls.len(),
        mark(rec_ok)
    ));
    out.push(format!(
        "UNSTRUCTURED|input_bits=64|candidates={}|structure=none|cost=one_sha256_per_candidate|json=0",
        1u128 << 64
    ));
    out.push(format!(
        "TRANSMEMO_VERIFY|status={}|transposition_rows={rows_t}|pathcount_rows={rows_p}|pathcount_jump_rows={rows_j}|json=0",
        if ok { "PASS" } else { "FAIL" }
    ));
    if timing {
        let n: u64 = 1 << 20;
        let t0 = Instant::now();
        for i in 0..n {
            sha256(&i.to_be_bytes());
        }
        let ns = t0.elapsed().as_nanos().max(1);
        let rate = u128::from(n) * 1_000_000_000 / ns;
        out.push(format!(
            "UNSTRUCTURED_TIMING|lane=rust|hashes={n}|ns={ns}|hashes_per_second={rate}|years_for_2^64={}|json=0",
            (1u128 << 64) / rate.max(1) / (365 * 24 * 3600)
        ));
    }
    print!("{}", rows_text(&out));
    ok
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if !args.iter().any(|a| a == "--self-test") {
        eprintln!("usage: transmemo --self-test [--timing]");
        std::process::exit(2);
    }
    if !self_test(args.iter().any(|a| a == "--timing")) {
        std::process::exit(1);
    }
}
