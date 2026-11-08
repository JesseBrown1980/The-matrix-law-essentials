# Vienna nanoparticle toy model: integer-only replay repair (November 8, 2026)

## Why this branch exists
The merged PR #7 and later README update left the main repository's integer-only audit red: the original Vienna approximation script used Python Decimal operations and ordinary division. It was a valid high-precision *toy-model* study, but contradicted the repository's global integer-only contract. This is a source-compliance and independent-reproduction repair; it is not a revision to the underlying experimental interpretation.

## Historical evidence retained
- Original Python source from PR #7 retained byte-for-byte at `data/vienna_nanoparticle/evidence/ORIGINAL-DECIMAL-SOURCE-20261106.txt`, original Git blob `cd153dfaae045ff85ee2026da82268ce2299d353`.
- Original `data/vienna_nanoparticle/MEASURED-20261106.hbp` retained unchanged, original Git blob `ee5065794fe55a879b1fb89afdcae8b8636e501a`.
- Both predictions files and synthetic observation figures remain untouched.
- The historical high-precision solver may still be read from the nonexecuted text archive; it is **not imported** by any new executable tool.

## Replacement executable
`python/matrix_tricks/vienna_nanoparticle_slices_20261106.py --self-test` independently recreates the same finite Fourier grating toy using a 100-decimal-digit scaled **integer** representation. It computes the Talbot ratio, three grating Fourier modes, Bessel and trigonometric Taylor series, all 20 direct spatial visibilities for 9/27/81/243/729 grid sites, all four spectral visibilities, the 43,046,721-site *virtual* alias theorem and 80 triadic free-flight phase comparisons. It reads the historical results only as **comparison targets**, checks their Git blob identities, and checks numerical agreement to absolute error <10^-12.

This new arithmetic is not byte-equivalent to the old 75-digit Decimal series. Last-digit differences are not claimed as a contradiction: the scientifically useful observables are compared at a declared tolerance. The historical numeric receipt is not overwritten. There are no real experimental collision/interferometer datasets in either executable; the model remains a deliberately simplified phase-only toy and does not model dipole forces, nanoparticle visibility fits, or quantum gravity.

## Acceptance
Require the new fast hosted workflow and the existing main verification workflow:
1. strict integer Python audit reports zero violations on the replacement;
2. the full 20+4 observable cross-check and 80 triadic-jump checks pass;
3. all preserved historical Git blobs match;
4. the repository-level `./verify.sh` is green, including its integer-only audit, rustfmt, Rust 1.81/clippy and exact operator checks.

Draft MLE PRs #5 and #6 must be brought up to date after main passes; do not treat their stale branch failures as solved until they rerun at their own final heads.
