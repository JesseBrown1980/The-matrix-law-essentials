#!/usr/bin/env python3
"""Band-gap certificate for the full 3x3 zero-momentum operator from
shellschur inertia counts plus Krylov-Bogoliubov residual intervals.

Inputs: the sealed GitRAM sector receipt (CLUSTER_LEVEL lines per sector and
coupling) and shellschur logs run with JE=1000, JB=1000*R, so t_milli is the
threshold in micro-units of the physical energy.

Logic (every step an exact integer comparison):
  N(t) = sum over the four C2v sectors of neg(H_s - t). A sector count missing
  at t is taken as 0 only if the same sector has count 0 at some t' >= t
  (inertia is monotone in t).
  Require N(t1) = N(tB) = 8 with t1 < tB. Then E1..E8 < t1 and E9 >= tB.
  Each residual interval [lo, hi] contains an eigenvalue. Every interval with
  hi < t1 contains one of E1..E8, so E8 >= max lo over those intervals.
  Any interval with lo >= t1 contains an eigenvalue >= E9, so
  E9 <= min hi over those intervals.
  Hence E9 - E8 in [tB - t1, min_hi - max_lo] (micro-units).
usage: shell_gap_cert.py SEALED_HBP R SHELL_LOG...        json=0
"""
import re
import sys


def main():
    sealed, r, logs = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    ivs = []
    pat = re.compile(r"YM3X3_S(-?1)_(-?1)_r%d_CLUSTER_LEVEL\|k=(\d+)\|.*interval=\[(\d+),(\d+)\]" % r)
    for line in open(sealed):
        m = pat.match(line)
        if m:
            ivs.append((int(m.group(4)), int(m.group(5)), m.group(1) + "," + m.group(2), int(m.group(3))))
    counts = {}
    for path in logs:
        for line in open(path):
            if "SHELL_SCHUR|" not in line:
                continue
            f = dict(kv.split("=", 1) for kv in line.strip().split("|") if "=" in kv)
            if f.get("result") != "PASS" or f.get("JE") != "1000" or f.get("JB") != str(1000 * r):
                continue
            mt = re.search(r"s33_(-?1)_(-?1)\.matrix", f["matrix"])
            counts[(mt.group(1) + "," + mt.group(2), int(f["t_milli"]))] = int(f["negative"])
    sectors = ("1,1", "1,-1", "-1,1", "-1,-1")
    ts = sorted({t for (_, t) in counts})

    def n_at(t):
        tot = 0
        for s in sectors:
            if (s, t) in counts:
                tot += counts[(s, t)]
            elif any(counts.get((s, u)) == 0 for u in ts if u >= t):
                tot += 0
            else:
                return None
        return tot

    good = [t for t in ts if n_at(t) == 8]
    if len(good) < 2:
        print("SHELL_GAP_CERT|r=%d|result=FAIL|reason=need_two_thresholds_with_N=8|json=0" % r)
        sys.exit(2)
    t1, tb = good[0], good[-1]
    below = [iv for iv in ivs if iv[1] < t1]
    above = [iv for iv in ivs if iv[0] >= t1]
    straddle = [iv for iv in ivs if iv[0] < t1 <= iv[1]]
    if not below or not above or straddle:
        print("SHELL_GAP_CERT|r=%d|result=FAIL|reason=intervals_straddle_or_missing|json=0" % r)
        sys.exit(2)
    e8_lo = max(iv[0] for iv in below)
    e9_hi = min(iv[1] for iv in above)
    e1_lo_bound = min(iv[0] for iv in below)
    detail = ",".join("%s:%d@%d" % (s, counts.get((s, t), 0), t) for t in (t1, tb) for s in sectors)
    print("SHELL_GAP_CERT|L=3x3|r=%d|N(t1)=8|N(tB)=8|t1=%d/1000000|tB=%d/1000000|sector_counts=%s|E8_in=[%d,%d)/1000000|E9_in=[%d,%d]/1000000|gap_E9_minus_E8_in=[%d,%d]/1000000|cluster_intervals_below_t1=%d|cluster_lowest_interval_lo=%d/1000000|result=PASS|json=0"
          % (r, t1, tb, detail, e8_lo, t1, tb, e9_hi, tb - t1, e9_hi - e8_lo, len(below), e1_lo_bound))


if __name__ == "__main__":
    main()
