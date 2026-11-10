// Seeded ribbon search: best tilted-column DP slab (slab::dp, scan tilt) + exact upright refill, missing squares dropped in near the
// ends, then basin hopping of ALL squares in the fixed box [0,w+f]x[0,h+f].
//   usage: ribbon_hop w h dm f seconds [seed] [th_lo th_hi]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "../fill.hpp"
#include "../slab.hpp"
#include "../window.hpp"
int main(int argc, char** argv) {
    const int w = atoi(argv[1]), h = atoi(argv[2]), dm = atoi(argv[3]); const double f = atof(argv[4]), secs = atof(argv[5]);
    const uint64_t seed = argc > 6 ? atoll(argv[6]) : 1; const double tlo = argc > 7 ? atof(argv[7]) : 18, thi = argc > 8 ? atof(argv[8]) : 40;
    const int n = w * h + dm; const double Sx = w + f, Sy = h + f, S = std::max(Sx, Sy);
    std::mt19937_64 rng(seed); std::uniform_real_distribution<double> U(0, 1);
    // 1. best seed over tilt
    int best = -1; double bth = 0; std::vector<Square> bT; std::vector<fill::P2> bU;
    for (double th = tlo; th <= thi; th += 0.25) {
        slab::Window win; win.x1 = Sx; win.y1 = Sy;
        auto F = slab::dp(S, th * M_PI / 180, {}, win, 0.02, 0.01, 12, true);
        std::vector<Square> T;
        for (auto& q : F.squares) { const double he = 0.5 * (std::fabs(std::cos(q.t)) + std::fabs(std::sin(q.t))); if (q.x - he > -0.01 && q.y - he > -0.01 && q.x + he < Sx + 0.01 && q.y + he < Sy + 0.01) T.push_back({q.x, q.y, q.t}); }
        auto R = fill::upright_fill(Sx, Sy, T, 1e-3, 30, rng);
        const int tot = (int)T.size() + (int)R.up.size();
        if (tot > best) { best = tot; bth = th; bT = T; bU = R.up; }
    }
    std::fprintf(stderr, "seed: theta %.2f, %zu tilted + %zu upright = %d (target %d)\n", bth, bT.size(), bU.size(), best, n);
    std::vector<double> v;
    for (auto& q : bT) { v.push_back(q.x); v.push_back(q.y); v.push_back(q.angle); }
    for (auto& c : bU) { v.push_back(c.x); v.push_back(c.y); v.push_back(0); }
    while ((int)v.size() / 3 > n) v.resize(v.size() - 3);
    while ((int)v.size() / 3 < n) { const bool top = U(rng) < 0.5; v.push_back(0.5 + (Sx - 1) * U(rng)); v.push_back(top ? Sy - 0.5 - 2 * U(rng) : 0.5 + 2 * U(rng)); v.push_back(U(rng) * M_PI / 2); }
    window::Env env; env.S = Sx + 2 * squeeze::kMargin; env.yoff = h - w; env.index();
    window::HopOptions o; o.seconds = secs; o.threads = 8; o.seed = seed; o.max_iter = 300; o.temp = 0.2; o.r_lo = 1.5; o.r_span = 4;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<double> out; double be = 0;
    const bool ok = window::hop(env, v, o, out, &be);
    std::printf("%dx%d+%d f=%.4f: %s (E %.2e) in %.0fs\n", w, h, dm, f, ok ? "FOUND" : "not found", be, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    if (ok) { std::FILE* fo = std::fopen("work/ribbon_found.txt", "a"); std::fprintf(fo, "G %d %d %d %.12f %d\n", w, h, dm, f, n); for (int i = 0; i < n; i++) std::fprintf(fo, "%.12f %.12f %.12f\n", out[3 * i], out[3 * i + 1], out[3 * i + 2]); std::fclose(fo); }
}
