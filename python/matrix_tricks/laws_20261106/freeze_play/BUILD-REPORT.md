Freeze and play both work. I added them to a copy of the engine, and every check except one ran as predicted. The prediction was committed before any code. One timing band was missed in a run taken while another workflow was loading the machine; a rerun with nothing else of mine running met it. The paused blind run was left alone: I listed its folder once, its newest file time is 12:20 local (before this lane started), and nothing in it was read, written or started. No run touched a reserved point: 76 jobs were launched through the guard, 0 were refused, and 0 command lines contain 3/5 or 2/27. The repo's stage driver ran 60 more runs from plans filtered by the same guard, with 0 reserved rows.

**What I did**
- Copied e4fc6ba into my folder. Its tree hash `7019067a…` matches e4fc6ba's, so the copy is exact. Engine source sha256 is `ed75eca8…` as expected.
- Local commits only, nothing pushed: `5fd9e24` base, `f105c0d` prediction, `f70d96b` code, `940bbfa` results.
- New file `src/bin/polariton_freeze_play.rs`. The four old modes (point, stream, bessel, rows) are the original lines, moved into functions the new modes share. The new modes:
  - **freeze:** runs as point does and writes a leaf every 9 windows. The tree levels are j = 2,3,4,5 (9, 27, 81, 243 windows).
  - **play:** checks the whole chain, rebuilds the HBP row from the parameters the leaf names and requires it byte for byte. Any mismatch exits 2 with no row.
  - **back:** runs e^{+iHt} as conj(forward window(conj ψ)), which holds because H is real. The certificate accumulators carry on from the leaf's.
  - **size:** prints a run's constants and leaf size without running any window.
- Freeze, play, back and size refuse reserved points in code. The driver also checks every launch, including leaves named for state comparison.
- SHA-256 comes from the repo's own `gguf_qnn::sha256`.

**Results (all MEASURED unless tagged)**

*Gate and reproduction*
- **Gate:** fmt, clippy -D warnings (on the binary and the whole package), 3 unit tests (one flips every byte of a small leaf and each flip is refused), and build all exit 0.
- **Point mode is byte-identical to the receipts:**
  - RECEIPT-GRID whole file: sha256 `4a1a277885015a5e937a028d767c617878ae9a05e75bade9eee3e14f0e07e480`, 20 runs, 2,079 s at 2 parallel.
  - The 58 allowed rows of RECEIPT-BASELINE: `ff57d5e8f8dc5697709f10ed42693ac066e45d534aa112b87ebdefbc1a4a9de0`, 40 runs, 739 s.
  - RECEIPT-DYADIC: `66e13a71d948dbadfc486d36b22015432e4e8eb92310b8492bb265ad53d94cf9`.
- **Stream mode** matches the dyadic rows and 16 of 16 small-baseline rows (depth 8/tile 8 and depth 3/tile 3), with array_bytes set aside.

*Freeze and play*
- **Freeze's stdout** equals the receipt at 6 of 6 points: N=100, 243, 729, 2187, and N=6561 at K=1200 and 1600.
- **Leaf counts** match the prediction: 10, 35, 21, 14, 14, 14.
- **Play equals the uninterrupted run:** 14 of 14 plays give the same rows byte for byte. They started 1 to 263 windows before the read, and include N=6561 at K=1200 and 1600 at t = 119/64, 2 and 9/4.
- **Leaves reproduce:** 31 of 31 leaves written again by a play from an earlier leaf are byte-identical to the original leaves.

*Tamper tests*
- One changed byte in the state, the im section, the HBP row, the index, the seal, the SHA or the SH, and a truncated file: each exits 2 with 0 stdout bytes, 10 of 10 (P1 plus two at P6).
- The reasons given were sha_mismatch, seal_mismatch and file_length.
- An unchanged copy plays and its 3 rows equal the uninterrupted ones.

*File size per state (MEASURED)*

Each leaf is 16 × even_states plus a head of 1,273 to 1,336 bytes.

| Point | Leaf bytes |
|---|---:|
| K=1200 (N=6561) | 5,780,518 to 5,780,538 |
| K=1600 (N=6561) | 10,266,921 to 10,266,940 |

Freezing costs 1.2% (K=1200) and 1.6% (K=1600) over a plain point run.

*Time saved: read from the nearest leaf against a run from t = 0*

| Point | Read at | Leaf | First run (loaded machine) | Alone, one process at a time | Predicted band |
|---|---|---|---:|---:|---|
| K=1200 | 9/4 | 126 | 6.02× | 15.5 s vs 104.9 s, **6.75×** | 6.0–8.5 |
| K=1200 | 119/64 | 117 | 33.2× | 2.2 s vs 87.2 s, **39.2×** | 20–60 |
| K=1600 | 9/4 | 126 | **5.79× (MISS)** | 27.7 s vs 204.3 s, **7.37×** | 6.0–8.5 |
| K=1600 | 119/64 | 117 | 33.8× | 4.4 s vs 168.6 s, **38.4×** | 25–60 |

- The first plays ran at load 4 to 5 while another workflow was busy; the alone reruns ran at load 0.9 to 1.0. Both sets are kept in RESULTS.hbp.
- At K=1600 the 10 MB leaf's chain is verified 0.18 s after the process starts, and setup plus its check against the engine are done at 1.0 s.

*Backward play*
- **Reads:** 11 of 11 backward reads have C1, A1 and V2 inside the forward read's certified interval. The largest gaps are |ΔC1| ≤ 7e-18, |ΔA1| ≤ 19e-18 and |ΔV2| ≤ 175e-18, all at most 1/1000 of the forward error.
- **Overlap (DERIVED, MEASURED):** C1 and A1 intervals overlap 11 of 11; V2 overlaps 8 of 11. The 3 others are at N=100, K=40, where eps_back is about 5.6 and the engine prints V2 = (0,−1), meaning no V2 bound.
- **State vectors:** played back to an earlier leaf, the state is inside its derived bound at 6 of 6. The bound is 100.8, 1,273, 1,765, 1,302, 3,964 and 5,942 times the l2 distance.

*Projection to the full N=6561 space (DERIVED)*

`size` measured E=2714, degree 90 and 10,768,242 even states for (1/2, 1/324).
- Leaf payload: 172,291,872 bytes; 14 leaves to t=2 come to 2.41 GB.
- The work at full space is 27.46 times the work at K=1600.
- At the alone rate that is about 39.2 s per window: t=2 from t=0 in about 5,019 s (1.39 h), and any read from the nearest level-2 leaf in at most 8 windows, about 314 s, plus about 3 s to load and verify the leaf.
- The point layout needs 1.29 GB of arrays.
- The measured K=1600/K=1200 rate ratio was 2.03 against 1.994 from states × degree, so the projection may run about 2% low.

**Files (sha256)**

All under `C:/tmp/opus55/freeze-play-20261106/`.
- `repo/src/bin/polariton_freeze_play.rs` `4a903f22f434bb407df85eb12d0cfaa1a2acf2bbdf1063d095a69be5f0e700f3`
- `target/release/polariton_freeze_play` `93e3ab7dbbfed2a3ed093af26fb3735673fa4da043c95d921d9631eec76de489`
- `target/release/polariton_first_correction` (built from `ed75eca8…`, never run) `d892971d6c0023534a2cbfac2631a7d983eeeab6853ffb6524b589e74795f270`

In `repo/polariton-first-correction/freeze-play/`:
- `PREDICTION.hbp` `abc12c66254c68180e2d898786b638940adf66fbf1dfba78bebd79e8003e3bd2`
- `RESULTS.hbp` (142 rows, including an FP_VERDICT for each prediction) `50f2a00f3fee9ca7ff7218a72f66e26e6c67c775a9ddb74d9fedd6b9b5b4c6df`
- `README.md` `6d8ffbd22be6fa8c27c6cfc4f5601a2927e558d286ace47803750bc7faa3b598`
- `check.py` `ab84bf08e06b935af772363c407d9114cb4d401c6c4970ed8f1968a0fb32cde8`
- `fp_run.py` `65766b0fb0a90e967feb3cb034566eca206ecca8ac048bdf1413a936962aba85`
- `gate.sh` `a4010e11704c7e068d361ac6ab1434d4df311f50cc229c744547931bf7972bea`
- `smoke.sh` `7d5dc2d519d279203f56fe822e312b2f83f77adf6fd9244f4c189a3b25c75b2b`
- `SHA256SUMS`, `RESULTS.hbp.sha256`

In `leaves/`:
- 139 leaves, 342 MB, listed in `LEAVES.SHA256SUMS` `10b5093d0dfb2dcb85beb3830c078ccbbfa91ccb772928dfa8738b8fb50ff423`
- `TREES.SHA256SUMS` `8c903d6a981f260c4287e392a96d8a26862edf4b9046cdf7523987e9c42e06e4`
- The six main trees:
  - P1 `739de726…d696`
  - P2 `e8f9fb38…84f7`
  - P3 `fbb0ccd9…d3d8`
  - P4 `ec7e2383…7328`
  - P5 `8b544e4c…2615`
  - P6 `4b7ba3e3…ab38`

Run logs and reproduced receipts are in `runs/`, including `JOBS-*.hbp` and `GATE.hbp`.

**For the seat to decide**
1. **Landing:** whether `polariton_freeze_play.rs` and the `freeze-play/` folder go into the GitHub repo, with a CI step. Everything is local now.
2. **Forgery:** the seal proves the file is intact, not who wrote it. A leaf with eps_tail changed and its chain resealed plays when no tree is given (exit 0, rows differ from the uninterrupted run). The sealed TREE.hbi, or a published seal, is what refuses it. Should play require `tree=` by default?
3. **The missed band:** TIME at K=1600, 9/4 read 5.79× in the first run under outside load and 7.37× alone. My central guesses (7.6×/7.7× and 40×/45×) were above the alone results.
4. **Disk:** keep or remove the 342 MB of leaves in my folder.
5. **Full-space run:** whether to run N=6561, K=6561 with leaves. It projects to about 1.4 h to t=2 and 2.4 GB of leaves at level 2.
6. **Slower play windows (CANNOT_SEE):** in the alone pairs, windows played from a leaf took 13–18% longer than windows from t=0 (829 vs 728 ms at K=1200; 1,484–1,680 vs 1,419 ms at K=1600). These runs don't show why.
