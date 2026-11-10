// usage: recompose_probe file.csv n delta... : refill count at L-delta vs number of tilted squares removed
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <atomic>
#include <mutex>
#include "../book.hpp"
#include "../recompose.hpp"
static double STEP = 0.02, TOL = 0.01, WM = 0.8;
int main(int argc, char** argv) {
    if (getenv("STEP")) STEP = atof(getenv("STEP")); if (getenv("TOL")) TOL = atof(getenv("TOL")); if (getenv("WM")) WM = atof(getenv("WM"));
    const Packing p = book_packing(argv[1], atoi(argv[2]));
    auto sp = recompose::split(p);
    std::printf("n=%d L=%.6f axis %zu tilted %d\n", atoi(argv[2]), sp.L, sp.axis.size(), sp.tilted);
    for (int a = 3; a < argc; a++) {
        const double delta = atof(argv[a]), S = sp.L - delta;
        std::vector<slab::Sq> obs = recompose::shrink(sp, delta);
        auto cl = recompose::clusters(sp);
        int total = 0, need = 0; std::string desc;
        for (auto& c : cl) {
            slab::Window w; w.x0 = c.x0 - WM; w.x1 = c.x1 + WM; w.y0 = c.y0 - WM; w.y1 = c.y1 + WM;
            // window is for squares' hulls: centre bbox +- 0.8 ~ original extents
            std::atomic<int> next{0}; std::mutex mu; int best = -1; double bth = 0;
            auto worker = [&]() { for (;;) { int i = next++; if (i >= 360) return; const double th = i * 0.5;
                auto F = slab::dp(S, th * M_PI / 180, obs, w, STEP, TOL, 40);
                std::lock_guard<std::mutex> lk(mu); if (F.count > best) { best = F.count; bth = th; } } };
            std::vector<std::thread> pool; for (int t = 0; t < 8; t++) pool.emplace_back(worker); for (auto& t : pool) t.join();
            auto F = slab::dp(S, bth * M_PI / 180, obs, w, STEP, TOL, 40, true);
            for (auto& q : F.squares) obs.push_back(q);
            total += F.count; need += (int)c.idx.size();
            char b[96]; std::snprintf(b, sizeof b, " [cluster %zu: th %.1f -> %d]", c.idx.size(), bth, F.count); desc += b;
        }
        std::printf("delta %.4f S %.5f refill %d vs removed %d  %s%s\n", delta, S, total, sp.tilted, desc.c_str(), total >= sp.tilted ? "   <-- OK" : "");
    }
}
