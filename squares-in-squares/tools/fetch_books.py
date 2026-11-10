#!/usr/bin/env python3
"""Download the openly available packings ("books") and normalise each into books/<name>.csv
(submission format: id,n,square_id,x,y,angle). Raw downloads go to books/raw/<name>/.

Nothing downloaded is executed; files are only read as data. Each book has a matching
src/solvers/solve_<name>.cpp that cites the source. Run from the competition directory:

    python3 tools/fetch_books.py [name ...]      # default: all
"""
import csv, io, json, math, os, shutil, subprocess, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BOOKS = os.path.join(ROOT, "books")
RAW = os.path.join(BOOKS, "raw")
HEADER = "id,n,square_id,x,y,angle"

# name -> (kind, ref, file inside the download)
KAGGLE = {
    "alexchilton_starter_kit": ("dataset", "alexchilton/squares-in-squares-starter-kit", "best.csv"),
    "prvsiyan_auto_updater": ("dataset", "prvsiyan/squares-auto-updater-live", "submission.csv"),
    "lightningv08_atlas_merge": ("dataset", "lightningv08/squares-in-squares-layout-atlas", "alex_merge.csv"),
    "dariushafshar_certified": ("dataset", "dariushafshar/squares-in-squares-certified-packings", "index_packings_dilated.csv"),
    "prvsiyan_gpu_window": ("dataset", "prvsiyan/a-square-deal-frozen-public-inputs", "gpu_window_v2.csv"),
    "aqibraza_top_lb_blend": ("kernel", "aqibrazadev/squares-in-squares-top-lb-blending-visualization", "submission.csv"),
    "kmmmvs_dashboard": ("kernel", "kmmmvs/squares-dashboard-public-batch-20261025", "submission.csv"),
    "nilanjan_portfolio_merge": ("kernel", "nilanjansaha123/s-in-s-certified-portfolio-merge-proj-quench", "submission.csv"),
    "lightningv08_cuda_search": ("kernel", "lightningv08/squares-in-squares-cuda-packing-search", "submission.csv"),
    "kirilltym_projunion": ("kernel", "kirilltym/squares-in-squares-cpu-cuda-mps-eda-projunion", "submission.csv"),
}


def kaggle_fetch(name):
    kind, ref, fname = KAGGLE[name]
    d = os.path.join(RAW, name)
    os.makedirs(d, exist_ok=True)
    src = os.path.join(d, fname)
    if not os.path.exists(src):
        if kind == "dataset":
            cmd = ["kaggle", "datasets", "download", ref, "-f", fname, "-p", d, "--unzip"]
        else:
            cmd = ["kaggle", "kernels", "output", ref, "-p", d]
        subprocess.run(cmd, check=True)
    return src


def normalise_csv(src, dst):
    """Keep only the six standard columns, rows grouped by n, ids renumbered."""
    with open(src, newline="") as f:
        r = csv.DictReader(f)
        need = ["n", "square_id", "x", "y", "angle"]
        if not all(c in r.fieldnames for c in need):
            raise SystemExit(f"{src}: unexpected columns {r.fieldnames}")
        rows = sorted(((int(x["n"]), int(x["square_id"]), x["x"], x["y"], x["angle"]) for x in r))
    with open(dst, "w") as f:
        f.write(HEADER + "\n")
        for i, (n, s, x, y, a) in enumerate(rows):
            f.write(f"{i},{n},{s},{x},{y},{a}\n")
    return len(rows)


# ---- Couzo: exact rational certificates -> doubles --------------------------------------------
COUZO_RAW = "https://raw.githubusercontent.com/franciscouzo/square-packing/master/certificates/best"
DILATE = 1e-11  # open exact contacts into real gaps that survive the float64 round trip


def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


def couzo_fetch_convert():
    d = os.path.join(RAW, "couzo")
    os.makedirs(d, exist_ok=True)
    man = os.path.join(d, "manifest.json")
    if not os.path.exists(man):
        open(man, "wb").write(fetch(f"{COUZO_RAW}/manifest.json"))
    wit = json.load(open(man))["witnesses"]
    ns = sorted(int(k) for k in wit if int(k) <= 1000)

    def get(n):
        p = os.path.join(d, f"n{n}.cert")
        if not os.path.exists(p):
            open(p, "wb").write(fetch(f"{COUZO_RAW}/n{n}.cert"))
        return n

    with ThreadPoolExecutor(16) as ex:
        list(ex.map(get, ns))

    out = []
    for n in ns:
        lines = [l for l in open(os.path.join(d, f"n{n}.cert")).read().splitlines() if l and not l.startswith("#")]
        assert int(lines[0].split()[0]) == n, (n, lines[0][:40])
        pts = []
        for l in lines[1:1 + n]:
            x, y, t = (float(Fraction(tok)) for tok in l.split())
            pts.append([x, y, 2.0 * math.atan(t)])
        # dilate centres about the middle of the bounding box of centres
        cx = (min(p[0] for p in pts) + max(p[0] for p in pts)) / 2
        cy = (min(p[1] for p in pts) + max(p[1] for p in pts)) / 2
        for p in pts:
            p[0] = cx + (p[0] - cx) * (1 + DILATE)
            p[1] = cy + (p[1] - cy) * (1 + DILATE)
        out.append((n, pts))
    dst = os.path.join(BOOKS, "couzo.csv")
    with open(dst, "w") as f:
        f.write(HEADER + "\n")
        i = 0
        for n, pts in out:
            for s, (x, y, a) in enumerate(pts):
                f.write(f"{i},{n},{s},{x!r},{y!r},{a!r}\n")
                i += 1
    return len(out)


def main():
    os.makedirs(RAW, exist_ok=True)
    want = sys.argv[1:] or (["couzo"] + list(KAGGLE))
    for name in want:
        print(f"== {name}", flush=True)
        if name == "couzo":
            print(f"   {couzo_fetch_convert()} packings -> books/couzo.csv")
        else:
            src = kaggle_fetch(name)
            print(f"   {normalise_csv(src, os.path.join(BOOKS, name + '.csv'))} rows -> books/{name}.csv")


if __name__ == "__main__":
    main()
