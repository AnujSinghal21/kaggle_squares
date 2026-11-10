#!/usr/bin/env python3
"""Export packings n_min..n_max of a submission as gadget records ("G w h dm f nsq" + "x y angle" lines) with w=h=floor(sqrt n),
dm = n - w*w, f = side - w, positions normalised so the bounding box starts at the origin.  Used by tools gftmc_test / FDP_GADGETS."""
import csv, math, sys
src, n0, n1, dst = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
pk = {}
with open(src, newline="") as f:
    for r in csv.DictReader(f):
        n = int(r["n"])
        if n0 <= n <= n1:
            pk.setdefault(n, []).append((int(r["square_id"]), float(r["x"]), float(r["y"]), float(r["angle"])))
with open(dst, "w") as out:
    for n in sorted(pk):
        sq = [t[1:] for t in sorted(pk[n])]
        he = [0.5 * (abs(math.cos(a)) + abs(math.sin(a))) for _, _, a in sq]
        mx = min(x - h for (x, _, _), h in zip(sq, he)); my = min(y - h for (_, y, _), h in zip(sq, he))
        K = math.isqrt(n)
        xs = [x - mx + h for (x, _, _), h in zip(sq, he)]; ys = [y - my + h for (_, y, _), h in zip(sq, he)]
        L = max(max(xs), max(ys))
        out.write(f"G {K} {K} {n - K*K} {L - K!r} {n}\n")
        for (x, y, a) in sq:
            out.write(f"{x - mx!r} {y - my!r} {a!r}\n")
print("wrote", len(pk), "puzzle gadgets")
