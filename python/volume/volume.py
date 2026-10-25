#!/usr/bin/env python3
"""Volume comparison through the HTTP-0 test ground: Chebyshev evolution on
the 2x3 torus (47,449 states) at the 2x2 rungs, step count doubled until
E0 and E1 each move by at most TOL micro-units between doublings.
Gap guess = the 2x2 matrix-jump gap at the same r. Integer-only, json=0."""
import os, sys, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))  # repository root
TOL = 5
GUESS = {0: 5999999, 1: 5608493, 2: 4002279, 3: 1473689, 5: 341237, 6: 257248, 8: 177326, 9: 154383,
         15: 88346, 21: 62323, 27: 48228, 45: 28794, 63: 20539, 81: 15966, 135: 9574, 189: 6838, 243: 5318}

def field(line, k):
    for kv in line.strip().split("|"):
        if kv.startswith(k + "="):
            return int(kv.split("=", 1)[1].split("/")[0])

def run(lat, r, steps):
    u = "http://127.0.0.1:5099/cheb/%s/1/%d/%d/%d" % (lat, r, steps, GUESS[r])
    line = urllib.request.urlopen(u, timeout=7200).read().decode()
    return field(line, "E0"), field(line, "E1"), field(line, "ms")

lat = sys.argv[1]
cap = int(sys.argv[2])
out = open(os.path.join(ROOT, "data", "volume", "volume_%s.hbp" % lat), "a")
for r in [int(x) for x in sys.argv[3:]]:
    steps, prev, total_ms = 500, None, 0
    while True:
        e0, e1, ms = run(lat, r, steps)
        total_ms += ms
        if prev and abs(e0 - prev[0]) <= TOL and abs(e1 - prev[1]) <= TOL:
            status = "settled"
            break
        if steps >= cap:
            status = "cap_reached"
            break
        prev, steps = (e0, e1), steps * 2
    out.write("VOLUME|lattice=%s|r=%d|steps=%d|E0=%d/1000000|E1=%d/1000000|gap=%d/1000000|status=%s|tol=%d|ms=%d|kind=MEASUREMENT_NOT_CERTIFICATE|json=0\n"
              % (lat, r, steps, e0, e1, e1 - e0, status, TOL, total_ms))
    out.flush()
    print("%s r=%-4d E0=%d.%06d E1=%d.%06d gap=%d.%06d (2x2 gap %d.%06d) %s steps=%d %ds"
          % (lat, r, e0 // 10**6, e0 % 10**6, e1 // 10**6, e1 % 10**6, (e1 - e0) // 10**6, (e1 - e0) % 10**6,
             GUESS[r] // 10**6, GUESS[r] % 10**6, status, steps, total_ms // 1000), flush=True)
