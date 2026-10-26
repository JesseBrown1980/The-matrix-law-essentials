#!/usr/bin/env python3
"""Certified E0, E1 and gap of a compressed operator (vacuum ball), integer-only.

1. `cluster` (Rust) resolves the whole low-energy cluster with a block
   Chebyshev-filtered subspace iteration and returns, for each level, an
   exact rational Rayleigh quotient and a rigorous residual radius:
   some eigenvalue lies in [theta_k - eps_k, theta_k + eps_k] (Krylov-Bogoliubov).
2. Ordering: thresholds t1 between intervals 0 and 1, t2 between intervals
   1 and 2; the interval-Sturm verifier counts eigenvalues below each. If
   count(t1) = 1 and count(t2) = 2, the interval below t1 holds E0 and the one
   between t1 and t2 holds E1, so
        E0 in I0,  E1 in I1,  gap in [lo1 - hi0, hi1 - lo0]
   holds rigorously for the compressed operator.
   Micro-unit thresholds use scaling invariance: JE = 1000, JB = 1000 r,
   threshold argument = x_micro (the verifier's threshold unit is 1/1000).
3. Without inertia (--no-inertia, for balls too large for the verifier) the
   intervals are reported with ordering=UNVERIFIED.

The compressed operator's E0, E1 are upper bounds on the full lattice's E0, E1
(Cauchy interlacing); its gap is NOT a bound on the full gap.

Band mode (--band): solve 12 levels and certify the separation of the
8-level vacuum cluster from the next band. With disjoint intervals around
levels 6..9, Sturm counts of exactly 7, 8 and 9 at thresholds between them
prove that the intervals hold E7 (top of the cluster) and E8 (bottom of the
band), so the band gap E8 - E7 and the cluster width E7 - E0 get rigorous
enclosures (E0 from a count of exactly 1 below a threshold above its interval).

usage: cluster_certify.py BALL_MATRIX R [--no-inertia] [--digits D] [--band]
Output json=0.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CLUSTER = os.path.join(ROOT, "target", "release", "cluster")
STURM = os.path.join(ROOT, "target", "release", "sturmrs")


def fields(line):
    return dict(kv.split("=", 1) for kv in line.strip().split("|")[1:] if "=" in kv)


def num(s):
    return int(s.split("/")[0])


def interval(s):
    a, b = s.split("]")[0].lstrip("[").split(",")
    return int(a), int(b)


def count_below(matrix, r, x_micro, digits):
    out = subprocess.run([STURM, matrix, "1000", str(1000 * r), str(x_micro), str(digits)],
                         capture_output=True, text=True).stdout
    f = fields(out)
    if f.get("status") != "COMPLETE" or f.get("ambiguous") != "0":
        return None
    return int(f["negative"])


def band(matrix, r, digits):
    out = subprocess.run([CLUSTER, matrix, "1", str(r), "12", "4"], capture_output=True, text=True).stdout
    levels, summary = [], None
    for line in out.splitlines():
        f = fields(line)
        if line.startswith("CLUSTER_LEVEL"):
            levels.append(interval(f["interval"]))
        elif line.startswith("CLUSTER|"):
            summary = f
    name = os.path.basename(matrix)
    if summary is None or summary.get("converged") != "1" or len(levels) < 10:
        print("BAND_CERT|matrix=%s|r=%d|status=NOT_CONVERGED|json=0" % (name, r))
        sys.exit(2)
    # need disjoint intervals: 0|1, 6|7, 7|8, 8|9
    for a, b in ((0, 1), (6, 7), (7, 8), (8, 9)):
        if not levels[a][1] < levels[b][0]:
            print("BAND_CERT|matrix=%s|r=%d|status=INTERVALS_OVERLAP|pair=%d,%d|json=0" % (name, r, a, b))
            sys.exit(2)
    th = {k: (levels[k - 1][1] + levels[k][0]) // 2 for k in (1, 7, 8, 9)}   # threshold below level k
    want = {1: 1, 7: 7, 8: 8, 9: 9}
    got = {k: count_below(matrix, r, th[k], digits) for k in th}
    (lo0, hi0), (lo7, hi7), (lo8, hi8) = levels[0], levels[7], levels[8]
    base = "BAND_CERT|matrix=%s|n=%s|r=%d|E0_in=[%d,%d]|E7_in=[%d,%d]|E8_in=[%d,%d]|cluster_width_in=[%d,%d]|band_gap_in=[%d,%d]|units=1/1000000|counts=%s" % (
        name, summary["n"], r, lo0, hi0, lo7, hi7, lo8, hi8, lo7 - hi0, hi7 - lo0, lo8 - hi7, hi8 - lo7,
        ",".join("t%d:%s" % (k, got[k]) for k in sorted(got)))
    if all(got[k] == want[k] for k in want):
        print(base + "|ordering=CERTIFIED|scope=compressed_operator_or_full_lattice_as_named|json=0")
    else:
        print(base + "|ordering=FAILED|json=0")
        sys.exit(2)


def main():
    matrix, r = sys.argv[1], int(sys.argv[2])
    if "--band" in sys.argv:
        digits = int(sys.argv[sys.argv.index("--digits") + 1]) if "--digits" in sys.argv else 40
        band(matrix, r, digits)
        return
    inertia = "--no-inertia" not in sys.argv
    digits = int(sys.argv[sys.argv.index("--digits") + 1]) if "--digits" in sys.argv else 40
    out = subprocess.run([CLUSTER, matrix, "1", str(r)], capture_output=True, text=True).stdout
    levels, summary = [], None
    for line in out.splitlines():
        f = fields(line)
        if line.startswith("CLUSTER_LEVEL"):
            levels.append((num(f["theta"]), num(f["eps_up"]), interval(f["interval"])))
        elif line.startswith("CLUSTER|"):
            summary = f
    name = os.path.basename(matrix)
    if summary is None or summary.get("converged") != "1" or len(levels) < 3:
        print("CLUSTER_CERT|matrix=%s|r=%d|status=NOT_CONVERGED|json=0" % (name, r))
        sys.exit(2)
    (t0, e0, (lo0, hi0)), (t1, e1, (lo1, hi1)), (_, _, (lo2, _)) = levels[0], levels[1], levels[2]
    if not (hi0 < lo1 and hi1 < lo2):
        print("CLUSTER_CERT|matrix=%s|r=%d|status=INTERVALS_OVERLAP|json=0" % (name, r))
        sys.exit(2)
    gap_lo, gap_hi = lo1 - hi0, hi1 - lo0
    base = "CLUSTER_CERT|matrix=%s|n=%s|r=%d|E0_in=[%d,%d]/1000000|E1_in=[%d,%d]/1000000|gap_in=[%d,%d]/1000000|outer=%s" % (
        name, summary["n"], r, lo0, hi0, lo1, hi1, gap_lo, gap_hi, summary["outer"])
    if not inertia:
        print(base + "|ordering=UNVERIFIED|scope=compressed_operator|json=0")
        return
    th1, th2 = (hi0 + lo1) // 2, (hi1 + lo2) // 2
    c1, c2 = count_below(matrix, r, th1, digits), count_below(matrix, r, th2, digits)
    if c1 == 1 and c2 == 2:
        print(base + "|t1=%d|count_t1=1|t2=%d|count_t2=2|ordering=CERTIFIED|scope=compressed_operator|json=0" % (th1, th2))
    else:
        print(base + "|t1=%d|count_t1=%s|t2=%d|count_t2=%s|ordering=FAILED|json=0" % (th1, c1, th2, c2))
        sys.exit(2)


if __name__ == "__main__":
    main()
