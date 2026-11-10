// From-scratch gadget search: w*h+dm free squares, random start, basin hopping at a fixed rectangular box [0,w+f]x[0,h+f] (window.hpp engine).
//   usage: gadget_scratch w h dm f seconds [seed] [init=random|grid]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include "../window.hpp"
int main(int argc, char** argv) {
    const int w = atoi(argv[1]), h = atoi(argv[2]), dm = atoi(argv[3]); const double f = atof(argv[4]), secs = atof(argv[5]);
    const uint64_t seed = argc > 6 ? atoll(argv[6]) : 1; const bool grid = argc > 7 && !strcmp(argv[7], "grid");
    const int n = w * h + dm;
    window::Env env; env.S = w + f + 2 * squeeze::kMargin; env.yoff = h - w; env.index();
    std::mt19937_64 rng(seed); std::uniform_real_distribution<double> U(0, 1);
    std::vector<double> v;
    for (int i = 0; i < n; i++) {
        if (grid) { const int c = i % w, r = i / w; v.push_back(0.5 + c * (w + f - 1) / std::max(1, w - 1) * (w > 1)); v.push_back(0.5 + r * (h + f - 1) / std::max(1.0, (double)(n / w))); v.push_back(U(rng) * 0.1); }
        else { v.push_back(0.5 + (w + f - 1) * U(rng)); v.push_back(0.5 + (h + f - 1) * U(rng)); v.push_back(U(rng) * M_PI / 2); }
    }
    window::HopOptions o; if (getenv("HOP_T")) o.temp = atof(getenv("HOP_T")); if (getenv("HOP_PL")) o.p_local = atof(getenv("HOP_PL")); if (getenv("HOP_PA")) o.p_all = atof(getenv("HOP_PA")); if (getenv("HOP_RLO")) o.r_lo = atof(getenv("HOP_RLO")); if (getenv("HOP_RSPAN")) o.r_span = atof(getenv("HOP_RSPAN"));
    const auto t0 = std::chrono::steady_clock::now(); o.seconds = secs; o.threads = 8; o.seed = seed; o.max_iter = getenv("HOP_IT") ? atoi(getenv("HOP_IT")) : 300; o.verbose = false;
    std::vector<double> out; double be = 0;
    const bool ok = window::hop(env, v, o, out, &be);
    if (ok) { std::FILE* fo = std::fopen("work/scratch_found.txt", "a"); std::fprintf(fo, "G %d %d %d %.12f %d\n", w, h, dm, f, n); for (int i = 0; i < n; i++) std::fprintf(fo, "%.12f %.12f %.12f\n", out[3 * i], out[3 * i + 1], out[3 * i + 2]); std::fclose(fo); }
    std::printf("%dx%d+%d f=%.4f (%d squares): %s (best E %.3e) in %.1fs\n", w, h, dm, f, n, ok ? "FOUND" : "not found", be, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}
