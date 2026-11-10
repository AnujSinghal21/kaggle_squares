// usage: slab_probe S [g_lo g_hi] [th_step] [step] [tol]
// core grid g1 x g2 at the lower-left, right slab (tilted columns), top slab; prints the best totals.
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <thread>
#include "../slab.hpp"
int main(int argc, char** argv) {
    const double S = atof(argv[1]);
    const int glo = argc > 2 ? atoi(argv[2]) : (int)S - 12, ghi = argc > 3 ? atoi(argv[3]) : (int)S - 3;
    const double ths = argc > 4 ? atof(argv[4]) : 0.5, step = argc > 5 ? atof(argv[5]) : 0.02, tol = argc > 6 ? atof(argv[6]) : 0.01;
    struct Job { int g1, g2; };
    std::vector<Job> jobs;
    for (int g1 = glo; g1 <= ghi; g1++) for (int g2 = glo; g2 <= ghi; g2++) jobs.push_back({g1, g2});
    std::atomic<size_t> next{0}; std::mutex mu; int best = 0;
    auto worker = [&]() {
        for (;;) {
            size_t i = next++; if (i >= jobs.size()) return;
            const int g1 = jobs[i].g1, g2 = jobs[i].g2;
            std::vector<slab::Sq> core;
            for (int a = 0; a < g1; a++) for (int b = 0; b < g2; b++) core.push_back({a + .5, b + .5, 0});
            slab::Window wr; wr.x0 = g1 - 1e-9;
            int bR = -1; double thR = 0; slab::Fill FR;
            for (double th = 0; th < 90; th += ths) {
                auto F = slab::dp(S, th * M_PI / 180, core, wr, step, tol, 40);
                if (F.count > bR) { bR = F.count; thR = th; }
            }
            FR = slab::dp(S, thR * M_PI / 180, core, wr, step, tol, 40, true);
            std::vector<slab::Sq> obs = core; for (auto& q : FR.squares) obs.push_back(q);
            slab::Window wt; wt.y0 = g2 - 1e-9;
            int bT = -1; double thT = 0;
            for (double th = 0; th < 90; th += ths) {
                auto F = slab::dp(S, th * M_PI / 180, obs, wt, step, tol, 40);
                if (F.count > bT) { bT = F.count; thT = th; }
            }
            slab::Window wt2; wt2.y0 = g2 - 1e-9;
            auto FT = slab::dp(S, thT * M_PI / 180, obs, wt2, step, tol, 40, true);
            std::vector<slab::Sq> all = obs; for (auto& q : FT.squares) all.push_back(q);
            const int extra = slab::fill_corners(S, all, tol);
            const int total = g1 * g2 + bR + bT + extra;
            std::lock_guard<std::mutex> lk(mu);
            if (total > best) { best = total; std::fprintf(stderr, "total %d = core %d (%dx%d) + R %d (th %.2f) + T %d (th %.2f) + corners %d\n", total, g1 * g2, g1, g2, bR, thR, bT, thT, total - g1 * g2 - bR - bT); }
        }
    };
    std::vector<std::thread> pool; for (int t = 0; t < 8; t++) pool.emplace_back(worker); for (auto& t : pool) t.join();
    printf("S=%.5f best %d\n", S, best);
}
