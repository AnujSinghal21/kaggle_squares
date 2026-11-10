// Ribbon end-cap search: exact periodic interior (frozen) + hopped end caps.
//   Observation (2026-10-28, from the stored 5x15+5 / 5x17+7 gadgets): the tilted interior of a 5-wide ribbon gadget is an exactly periodic
//   sheared lattice.  Six columns of squares tilted by -theta; column j has centre x = (c+s)/2 + j*c, y = y0 - j*s + k/c (c=cos, s=sin),
//   so it fills the strip width W = 6c + s exactly (flush against both walls) at 6c squares per unit height (5.35 vs 5.0 upright).
//   theta is fixed by W = 5+f, so the interior is a closed-form one-parameter family; only the two end caps (~23 squares each) are unknown.
//   This tool freezes the interior rows with centre y in [ylo,yhi], keeps the rest as free squares (nb at the bottom, the remainder at
//   the top, random start) and basin-hops them (window.hpp) in the box [0,5+f]x[0,h+f].
//   Known-good benchmark: 5 17 7 0.8024 (23/46/23 split) should re-find the stored ends; then lower f to what the leaders have.
//   RIB_W=<w> selects other widths (w+1 columns, W = (w+1)c + s).
//   usage: ribbon_ends h dm f ylo yhi y0 nb seconds [seed]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "../window.hpp"
int main(int argc, char** argv) {
    if (argc < 9) { std::fprintf(stderr, "usage: ribbon_ends h dm f ylo yhi y0 nb seconds [seed]\n"); return 1; }
    const int w = getenv("RIB_W") ? atoi(getenv("RIB_W")) : 5, h = atoi(argv[1]), dm = atoi(argv[2]); const double f = atof(argv[3]), ylo = atof(argv[4]), yhi = atof(argv[5]), y0 = atof(argv[6]);
    const int nb = atoi(argv[7]); const double secs = atof(argv[8]); const uint64_t seed = argc > 9 ? atoll(argv[9]) : 1;
    const int n = w * h + dm; const double W = w + f, H = h + f;
    // theta from 6cos+sin = W (theta in (0.2,0.9))
    double lo = 0.2, hi = 0.9;
    for (int i = 0; i < 100; i++) { const double m = 0.5 * (lo + hi); ((w + 1) * std::cos(m) + std::sin(m) > W) ? lo = m : hi = m; }
    const double th = 0.5 * (lo + hi), c = std::cos(th), s = std::sin(th);
    window::Env env; env.S = W + 2 * squeeze::kMargin; env.yoff = h - w;
    for (int j = 0; j <= w; j++)
        for (int k = -40; k <= 40; k++) {
            const double x = (c + s) / 2 + j * c + squeeze::kMargin, y = y0 - j * s + k / c;
            if (y >= ylo && y <= yhi) env.frozen.push_back({x, y, -th});
        }
    env.index();
    const int nf = (int)env.frozen.size(), nfree = n - nf;
    std::printf("theta=%.5f (%.2f deg) frozen=%d free=%d (bottom %d top %d) W=%.5f H=%.5f\n", th, th * 180 / M_PI, nf, nfree, nb, nfree - nb, W, H);
    if (nfree < 0 || nb < 0 || nb > nfree) return 1;
    std::mt19937_64 rng(seed); std::uniform_real_distribution<double> U(0, 1);
    std::vector<double> v;
    for (int i = 0; i < nfree; i++) {
        const bool bot = i < nb;
        const double ya = bot ? 0.5 : yhi, yb = bot ? ylo : H - 0.5;
        v.push_back(0.5 + (W - 1) * U(rng)); v.push_back(ya + (yb - ya) * U(rng)); v.push_back(U(rng) * M_PI / 2);
    }
    window::HopOptions o; o.seconds = secs; o.threads = 8; o.seed = seed; o.max_iter = getenv("HOP_IT") ? atoi(getenv("HOP_IT")) : 300;
    o.temp = getenv("HOP_T") ? atof(getenv("HOP_T")) : 0.2; o.r_lo = getenv("HOP_RLO") ? atof(getenv("HOP_RLO")) : 1.2; o.r_span = getenv("HOP_RSPAN") ? atof(getenv("HOP_RSPAN")) : 3;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<double> out; double be = 0;
    const bool ok = window::hop(env, v, o, out, &be);
    if (ok) {
        std::FILE* fo = std::fopen("work/ribbon_found.txt", "a");
        std::fprintf(fo, "G %d %d %d %.12f %d\n", w, h, dm, f, n);
        for (const auto& q : env.frozen) std::fprintf(fo, "%.12f %.12f %.12f\n", q.x, q.y, q.t);
        for (int i = 0; i < nfree; i++) std::fprintf(fo, "%.12f %.12f %.12f\n", out[3 * i], out[3 * i + 1], out[3 * i + 2]);
        std::fclose(fo);
    }
    std::printf("%dx%d+%d f=%.4f: %s (best E %.3e) in %.1fs\n", w, h, dm, f, ok ? "FOUND" : "not found", be, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}
