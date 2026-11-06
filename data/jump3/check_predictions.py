#!/usr/bin/env python3
"""Check the frozen PREDICT rows against the self-test outputs, one PREDICT_CHECK row each.

A PREDICT row names a row tag, a where-clause (field=value pairs joined by &),
a field, an operation and a value:
  op=eq      the field's text equals the value
  op=prefix  the field's text starts with the value
  op=le      the field is an integer at most the value
result=MET or MISSED; NO_ROW when no output row matches the where-clause and
AMBIGUOUS when more than one does. Every PREDICT row is reported, met or not.

usage: check_predictions.py PREDICTIONS.hbp OUTPUT.hbp [OUTPUT.hbp ...]
"""
import sys


def parse(line):
    parts = line.rstrip("\n").split("|")
    fields = {}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            fields[k] = v
    return parts[0], fields


def main():
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        return 2
    preds = []
    with open(sys.argv[1], encoding="ascii") as f:
        for line in f:
            tag, fields = parse(line)
            if tag == "PREDICT":
                preds.append(fields)
    rows = []
    for path in sys.argv[2:]:
        with open(path, encoding="ascii") as f:
            for line in f:
                if line.strip():
                    rows.append(parse(line))
    tally = {}
    out = []
    for p in preds:
        where = [w.split("=", 1) for w in p["where"].split("&")]
        hits = [f for tag, f in rows if tag == p["row"] and all(f.get(k) == v for k, v in where)]
        got = "none"
        if not hits:
            result = "NO_ROW"
        elif len(hits) > 1:
            result = "AMBIGUOUS"
        else:
            got = hits[0].get(p["field"], "none")
            if p["op"] == "eq":
                ok = got == p["value"]
            elif p["op"] == "prefix":
                ok = got.startswith(p["value"])
            elif p["op"] == "le":
                ok = got.isdigit() and int(got) <= int(p["value"])
            else:
                raise ValueError("unknown op " + p["op"])
            result = "MET" if ok else "MISSED"
        key = (p["basis"], result)
        tally[key] = tally.get(key, 0) + 1
        out.append("PREDICT_CHECK|id=%s|row=%s|where=%s|field=%s|op=%s|predicted=%s|got=%s|result=%s|basis=%s|json=0"
                   % (p["id"], p["row"], p["where"], p["field"], p["op"], p["value"], got, result, p["basis"]))
    results = ("MET", "MISSED", "NO_ROW", "AMBIGUOUS")
    total = {r: sum(v for (b, rr), v in tally.items() if rr == r) for r in results}
    by_basis = ",".join("%s:%s" % (b, "/".join("%d" % tally.get((b, r), 0) for r in results))
                        for b in ("KNOWN_FROM_SEAT", "DERIVED", "CONJECTURE"))
    out.append("PREDICT_SUMMARY|predictions=%d|met=%d|missed=%d|no_row=%d|ambiguous=%d|"
               "by_basis_met/missed/no_row/ambiguous=%s|json=0"
               % (len(preds), total["MET"], total["MISSED"], total["NO_ROW"], total["AMBIGUOUS"], by_basis))
    sys.stdout.write("".join(r + "\n" for r in out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
