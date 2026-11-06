**VERDICT: PASS_WITH_HOLDS**

The Build lane's main claims reproduce on my own points with my own checks. There are 8 holds, none at HIGH: 2 MEDIUM, 3 LOW, 3 NOTE. The paused blind-run folder was never opened, listed or touched. I made no git push and no GitHub writes. I only read the Build folder. Every launch went through a guard that refuses reserved points: 87 driver jobs and 551 tamper launches, 0 refused, and 0 command lines containing 3/5, 2/27 or 6/10. float_used=0.

**What I did**
- Copied e4fc6ba into `C:/tmp/opus55/freeze-play-check-20261106/repo`. Its tree is `7019067a…`, the same as e4fc6ba's.
- Added the Build lane's 940bbfa files on top of it. Compared with the base, 940bbfa adds 11 files and changes 0. Their sha256 values match the report.
- In the first 1,497 lines, the copy differs from e4fc6ba's engine in only 3 places: the `forbid(unsafe)` line moved, a `use sha256` line was added, and one doc comment changed.
- I chose points A to D (all allowed) and wrote my own verifier and mutation battery in Python, built in WSL with Rust 1.81.

**Results**

*Gate and build (MEASURED)*
- fmt, clippy -D warnings on the binary and on the whole package, 3 of 3 unit tests, and the release build all exit 0.
- Built in a different target folder, both binaries are byte-identical to the Build lane's: `93e3ab7d…e489` and `d892971d…f270`.

*Point and stream stay byte-identical to e4fc6ba (MEASURED)*
- 20 of 20 stdout pairs from e4fc6ba's engine and the copy are byte-equal, array_bytes included. They cover:
  - point at A, B, C and D;
  - stream at A and B with the default settings, depth 3/tile 3 and depth 5/tile 7;
  - the 4 sabotages in both point and stream at A;
  - bessel and rows.
- RECEIPT-DYADIC reproduces with sha256 `66e13a71…`.

*Freeze and play (MEASURED)*
- Freeze's stdout equals e4fc6ba's point output at 5 of 5 runs (A, B, C, D, and A with levels 1,3).
- Leaf counts are 21, 14, 14, 14 and 64, as predicted from the window counts.
- A second freeze at D wrote leaves byte-equal to the first, 14 of 14.
- 23 plays produced 41 rows, and all 41 are byte-equal to the uninterrupted run. The plays started 1 to 183 windows before the read. They include a play two generations from t = 0 and reads at off-grid times (65/64, 5/2, 73/64, 127/64).
- My binary played 3 of the Build lane's leaves, and all 6 rows equal RECEIPT-DYADIC.
- 96 of 96 leaves rewritten by a play are byte-equal to the originals.

*Independent leaf verifier (MEASURED)*
- My verifier checks each leaf's chain and file-name/size/sidecar consistency. It also checks three equalities against the window count: eps_poly = window × rho, eps_round = window × per_window_units, and applies = run_degree × (2·window − 1). It recomputes each tree's level and parent for every leaf.
- It passes 223 of 223 of my leaves and 139 of 139 of the Build lane's.
- All 9 of my trees pass. All 9 of the Build lane's trees pass, and so do their SHA256SUMS files.

*Backward play*
- DERIVED, my own working: the backward window is G(x) = conj(F(conj x)), and conj(U_P) = U_P⁻¹ because H_P is real and symmetric. Each window F is within ρ + pwu of U_P. Adding up the steps gives ‖Gⁿ(Fⁿx) − x‖ ≤ 2n(ρ + pwu), which is the engine's forward increment plus its backward increment.
- MEASURED: the engine's bound equals my formula at 16 of 16 states, and l2 is inside it at 16 of 16. The bound is 35× to 4,660× the l2 distance.
- MEASURED: 11 of 11 back reads have C1, A1 and V2 inside the forward read's interval. The largest gaps are |ΔC1| ≤ 24e-18, |ΔA1| ≤ 25e-18 and |ΔV2| ≤ 32e-18. The epsilon carried back matches my arithmetic exactly, 11 of 11.

*Tamper and forgery (MEASURED)*
- 13 kinds of unsealed change at A w63 and D w126 were all refused: 26 of 26 exit 2 with 0 stdout bytes. Reasons: sha_mismatch, seal_mismatch, file_length, head_not_printable_ascii. The 4 unchanged controls play with rows equal.
- Changes with the chain resealed:

| Change | Result |
|---|---|
| eps_tail, no tree | plays, rows differ (2/2) |
| eps_tail, original tree | refused (2/2) |
| eps_tail, **tree rebuilt and resealed** | **plays, rows differ (2/2)** |
| state | plays (2/2) |
| eps_poly | plays (2/2) |
| window + 1 | plays at D, refused at A because the read fell on the wrong side |
| SH only | refused, sh_mismatch (2/2) |

*Planted faults (MEASURED)*
23 faults, each caught. The unmutated engine fails 0 of the 10 checks.

| Fault | Caught by |
|---|---|
| M01 leak_sum missing | play rows, verifier |
| M02 tail_max missing | c10 key-list check only |
| M03 play resets eps | play rows, leaf repro |
| M04 back adds nothing to eps | c8 eps arithmetic only |
| M05 stale window in HBP | tree= refuses, verifier, play rows |
| M06 stale window in back | c8 (2 reads instead of 3) |
| M07 SH check skipped | resealed-SH case only |
| M08 SHA check skipped | unit test, tamper |
| M09 seal check skipped | unit test, tamper |
| M10 seal omits HBP | unit test, verifier, tamper |
| M11 Im read from Re | unit test, play rows |
| M12 SHA skips word 0 | unit test, first_state_byte |
| M13 Im index one word early | unit test, play rows |
| M14 back one negation | c8: reads leave the interval |
| M15 back no conjugation | c8: reads leave the interval |
| M16 bound forward only | c8 derived bound only |
| M17 parent index | verifier tree check only |
| M18 −3/5 guard dropped | unit test |
| M19 tree= ignored | tamper with the original tree |
| M20 play drops Im | play rows |
| M21 check_leaf always ok | unit test |
| M22 freeze off by one | freeze exits 101 |
| M23 leak_max missing | c10 key-list check only |

M02, M07 and M23 got past my first nine checks. I added c10 (HBP key list) and the resealed-SH case after seeing that, and both are kept in the battery.

*Timing (MEASURED, but another workflow had the load at 5 to 7)*
- At D, a play from w126 reads t = 2 in 1.96 s, against 88.1 to 94.5 s from t = 0 (about 45×).
- Freeze took 85.3 s. The 14 leaves of 3.26 MB took 1,274 ms to write.
- From w54, windows ran at 626 ms each, against 688 to 738 ms each from t = 0. The Build lane's question about slower play windows is still open (CANNOT_SEE).

**HOLDs**

HOLD H1 | file=polariton-first-correction/freeze-play/README.md and src/bin/polariton_freeze_play.rs (from_leaf) | severity=MEDIUM
OLD<<<
A leaf with eps_tail changed and its chain resealed plays, and its rows differ from the uninterrupted ones; with `tree=` it is refused, and so is a tree with one changed byte.
>>>
NEW<<<
A leaf with eps_tail changed and its chain resealed plays, and its rows differ from the uninterrupted ones; with `tree=` naming the TREE.hbi the freeze wrote it is refused, and so is that tree with one changed byte. A TREE.hbi rebuilt around the forged leaf and resealed plays it too, since a tree's seal is computed from its own lines; the refusal holds when the tree's sha256 is held outside it: `tree_sha256=` checked against TREES.SHA256SUMS committed with the leaves.
>>>
plus in from_leaf: `let mut tree_sha = None;` and a `tree_sha256=` argument; after `check_tree(...)`: `if let Some(want) = &tree_sha { if std::fs::read(t).map(|b| sha256::hex_digest(&b)).unwrap_or_default() != *want { refuse("tree_sha256_mismatch"); } }`
REASON: MEASURED at 2 of 2 points (A w63, D w126): a forged eps_tail plus a resealed forged TREE.hbi gives exit 0 with rows that differ from the uninterrupted run.

HOLD H2 | file=polariton-first-correction/freeze-play/check.py (section 8) | severity=MEDIUM
OLD<<<
            l2, bound = int(s["l2_diff_s62"]), int(s["bound_s62"])
>>>
NEW<<<
            l2, bound = int(s["l2_diff_s62"]), int(s["bound_s62"])
            hf = from_leaf_hbp(name)  # FP_HBP fields of the from-leaf
            derived = 2 * (int(s["from_window"]) - int(s["window"])) * (int(hf["rho_s62"]) + int(hf["per_window_units_s62"]))
            bound_ok = int(bound == derived and l2 < derived)
>>>
For FP_BACK_READ, also require: eps_poly_e18 == e18_s60(((2b−a)·rho)>>2), the same for eps_round with per_window_units, eps_tail above the leaf's value, and applies == the leaf's applies + (b−a)·2·run_degree. BACK_STATE and BACK_INTERVAL should require these too. My reference implementation is `scripts/ck_check.py` sec_back.
REASON: M04 and M16 keep 3 of 3 back reads inside the forward intervals. check.py relies on the engine's own printed bound, so it would see them only through the 100·l2 ≤ bound margin (DERIVED from reading check.py; CANNOT_SEE as a run). My derived check catches both (MEASURED).

HOLD H3 | file=src/bin/polariton_freeze_play.rs (check_leaf) | severity=LOW
OLD<<<
fn check_leaf(leaf: &Leaf, s: &Setup, ids: &Ids) -> Result<(), String> {
    let mine = hbp_row(s, ids, &leaf.acc);
>>>
NEW<<<
fn check_leaf(leaf: &Leaf, s: &Setup, ids: &Ids) -> Result<(), String> {
    let w = leaf.acc.done;
    let applies = i128::try_from(s.run_degree).unwrap() * (2 * w - 1);
    if leaf.acc.eps_poly != w * s.rho
        || leaf.acc.eps_round != w * s.per_window_units
        || i128::from(leaf.acc.applies) != applies
    {
        return Err("accumulators_differ_from_window".to_string());
    }
    let mine = hbp_row(s, ids, &leaf.acc);
>>>
REASON: MEASURED: a resealed eps_poly (2/2) and a resealed window + 1 (at D) play with exit 0 and rows that differ. All 362 genuine leaves satisfy the three equalities. The applies equality holds for runs whose Im is non-zero after window 1, which was every leaf measured.

HOLD H4 | file=src/bin/polariton_freeze_play.rs (tests) | severity=LOW
OLD<<<
        assert!(check_leaf(&leaf, &s, &ids).is_ok());
>>>
NEW<<<
        assert!(check_leaf(&leaf, &s, &ids).is_ok());
        for key in ["window", "eps_tail_s62", "eps_poly_s62", "eps_round_s62", "leak_max_s62", "leak_sum_s62", "tail_max_s42", "applies"] {
            assert!(leaf.hbp.contains(&format!("|{key}=")), "{key} missing");
        }
        let head = String::from_utf8_lossy(&bytes[..bytes.len().min(4096)]).to_string();
        let at = head.find("FP_SH|of=sha256|sh=").unwrap() + 19;
        let mut forged = bytes.clone();
        forged[at..at + 16].copy_from_slice(b"0123456789abcdef");
        let seal_at = head.find("FP_HASH|over=HBI,HBP,SHA,SH|seal=").unwrap();
        let seal = sha256::hex_digest(&forged[..seal_at]);
        forged[seal_at + 33..seal_at + 97].copy_from_slice(seal.as_bytes());
        assert_eq!(verify_leaf(&forged).err().as_deref(), Some("sh_mismatch"));
>>>
REASON: M02 and M23 drop a running maximum from the leaf. At point A both maxima fall after the leaves, so play rows stay byte-equal. M07 passes every unsealed tamper because the seal covers the SH line (MEASURED).

HOLD H5 | file=src/bin/polariton_freeze_play.rs (Freezer::new) | severity=LOW
OLD<<<
        if let Err(e) = std::fs::create_dir_all(dir) {
            refuse(&format!("dir_not_made_{e}"));
        }
>>>
NEW<<<
        if let Err(e) = std::fs::create_dir_all(dir) {
            refuse(&format!("dir_not_made_{e}"));
        }
        if std::path::Path::new(dir).join("TREE.hbi").exists() {
            refuse("dir_holds_a_tree");
        }
>>>
REASON: MEASURED on a copy of A's folder: a play from w126 with dir= that folder replaced its TREE.hbi (21 leaves listed became 7). After that, the genuine w63 with tree= was refused (tree_has_no_leaf_at_this_window).

HOLD H6 | file=src/bin/polariton_freeze_play.rs (from_leaf, after times are sorted) | severity=NOTE
OLD<<<
    times.sort_unstable_by_key(|&(windows, _)| if back { -windows } else { windows });
>>>
NEW<<<
    times.sort_unstable_by_key(|&(windows, _)| if back { -windows } else { windows });
    let lowest = times.iter().map(|&(w, _)| w).min().unwrap_or(0);
    if states.iter().any(|l| l.acc.done < lowest || l.acc.done >= leaf.acc.done) {
        refuse("state_leaf_outside_the_windows_played_back");
    }
>>>
REASON: MEASURED: back from A w189 to t = 2 with state=w63 exits 0 and prints no FP_BACK_STATE row.

HOLD H7 | file=polariton-first-correction/freeze-play/README.md | severity=NOTE
OLD<<<
Scope: forward plays reproduce the uninterrupted integers, so their rows are the same bytes; a back read is a different integer computation, bounded by the same certificate.
>>>
NEW<<<
Scope: forward plays reproduce the uninterrupted integers, so their rows are the same bytes. A back read is a different integer computation; its epsilon is the leaf's plus the windows played back, so it is at least the epsilon at the leaf's later time (measured 7× to 6.3×10^11× the forward epsilon at the same t, and above 1 at N = 150, K = 60). A certified read at time t is a forward play from the nearest leaf at or before t.
>>>
REASON: DERIVED from MEASURED rows. The back tail terms are 43% to 62% of the back epsilon. They can be dropped, because playing the block evolution back undoes the forward evolution exactly: ‖ψ_back(a) − ψ_full(a)‖ ≤ eps_tail(b) + eps_poly(b) + eps_round(b) + (b − a)(ρ + pwu). That would make the back epsilon 1.76× to 2.64× smaller, but it would still be at least the leaf's.

HOLD H8 | file=polariton-first-correction/freeze-play/README.md | severity=NOTE
OLD<<<
the window, the engine's source sha256)
>>>
NEW<<<
the window, the engine's source sha256, which binds the arithmetic that ran, and the sha256 of polariton_first_correction.rs it was built beside)
>>>
REASON: DERIVED from the source. `integrator_source_sha256` hashes `include_bytes!("polariton_first_correction.rs")`, but the code that runs is the copy inside `polariton_freeze_play.rs`.

**Files (sha256), all under `C:/tmp/opus55/freeze-play-check-20261106/`**
- `CHECK-RESULTS.hbp` (466 rows) `2e2ee6d2699a6ae18087f55ec82db780db8f6e7a3527469341e72fd778c56e30`
- `CHECK-HEAD.hbp` `08f3353bb14e95cb6d55958e89bcbf697281a455c904d542e27cbfaa67e8e829`
- `SHA256SUMS` `e85f98b07906675a4aeb304853866dfd92e1d8df174f19a46e7038eff8a953e7`
- Scripts, in `scripts/`:
  - `ck.py` `87737150bf09e0f30ef8bc5207f4fbfa05551a7e57fcbe5ba76d6647f6306f79`
  - `ck_check.py` `312de2de20aeb7e8f79858130540c5d4c698569eaa34c03e8d787900eb90a297`
  - `ck_tamper.py` `269dd4641a7d9921418717478b467b6352e0833cc5287048c6e4d23d287af7a5`
  - `ck_mutants.py` `b3b6c01aa8083489d78666f848bf05fc9bd152b084c44eeae6e5c3178a25fc74`
  - `ck_mut.py` `acb4278fb60328c48bc037cefd96fc7530c16e39bc6018c573bab843fc3074b4`
  - `ck_theirs.py` `34edd50540c993362afa7d8a659eec6a0df1ed749ae3656193dfff9af7009e2f`
  - `ck_backtight.py` `457ff7b5bc6301c35e560ea5d433a5d1521fb666b613bb249e40d86fa666660b`
  - `gate.sh` `7666bf98db1fbf2555e5728c41c73b8e9eb4a5ee73b58d04cca0651c6bcaa26b`
  - `build_mutants.sh` `eacc9da3986e1560cdfa9b2fd6b258ae18c0056f8bc6dd062bb30f9df070227d`
- Mutation records, in `mutants/`:
  - `MUTANTS.hbp` `3c1d23c33b0d1f16b840c552429a51136460fbd94c7219311d1102d1ea3a9048`
  - `BUILD.hbp` `cf891a2522be9750d1dfbee551d2efed9d4d16d2109430ba1892ce1696a5b02d`
  - `BATTERY-ALL.hbp` `ae6bf024ec91b77cd0358450b00d42601488e436c7805452cbfc609d537bf891`
- `runs/tamper-main/TAMPER.hbp` `53e0a0db34417e7bf5ebd7fb49c10b76009fb75d3bdfdea17a0464a42796dde0`
- The `runs/JOBS-*.hbp` files and both binaries are listed with their sha256 in `SHA256SUMS`.
- Every script and result file is LF-only: 0 CR bytes.

**For the seat to decide**
1. **Forgery anchor (H1):** whether play and back should require `tree=` together with `tree_sha256=` by default. Without the sha256, a resealed tree lets a forged leaf play.
2. **Before landing:** whether H3 (the three equalities in check_leaf), H4 (the extra unit test) and H2 (the derived back checks in check.py) go in first. The copy's point, stream and rows modes run any point, exactly as e4fc6ba's engine does, so a CI landing needs the launch filter kept.
3. **Back certificate (H7):** keep the current valid but loose epsilon, or adopt the tighter one without the back tail terms. Either way, the README scope text needs updating.
4. **Disk:** my folder is 259 MB: 118 MB of leaves, 47 MB of mutants (23 binaries), plus the build targets. Keep it or remove it.
5. **Timing question 6** is still open. My D timings ran while another workflow had the load at 5 to 7, so they don't settle it.
