// Ribbon generator: tilted columns (slab::dp, any tilt, columns slide independently) in the box [0,w+f]x[0,h+f] with flat walls,
// then exact upright refill around them (fill::upright_fill).  Reports the best total vs the target w*h+dm.
//   usage: ribbon_dp w h dm f [step] [tol]
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <atomic>
#include <mutex>
#include "../slab.hpp"
#include "../fill.hpp"
int main(int argc, char** argv) {
    const int w = atoi(argv[1]), h = atoi(argv[2]), dm = atoi(argv[3]); const double f = atof(argv[4]);
    const double step = argc > 5 ? atof(argv[5]) : 0.02, tol = argc > 6 ? atof(argv[6]) : 0.004;
    const double Sx = w + f, Sy = h + f, S = std::max(Sx, Sy);
    // corner blocks: nx x ny uprights anchored at each of the 4 corners (0..mx each way)
    const int mx = argc > 7 ? atoi(argv[7]) : 2;
    const int nb = (mx + 1) * (mx + 1);
    long total_cfg = 1; for (int c = 0; c < 4; c++) total_cfg *= nb;
    const double th_lo = argc > 8 ? atof(argv[8]) : 10, th_hi = argc > 9 ? atof(argv[9]) : 50;
    std::atomic<long> next{0}; std::mutex mu; int best = 0;
    auto worker = [&]() {
        for (;;) {
            const long ci = next++; if (ci >= total_cfg) return;
            std::vector<slab::Sq> obs; int nup = 0; long c = ci;
            for (int corner = 0; corner < 4; corner++) {
                const int b = (int)(c % nb); c /= nb; const int nx = b % (mx + 1), ny = b / (mx + 1);
                for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) {
                    const double x = (corner & 1) ? Sx - 0.5 - i : 0.5 + i, y = (corner & 2) ? Sy - 0.5 - j : 0.5 + j;
                    bool dup = false; for (auto& o : obs) if (std::fabs(o.x - x) < 0.99 && std::fabs(o.y - y) < 0.99) dup = true;
                    if (!dup) { obs.push_back({x, y, 0}); nup++; }
                }
            }
            for (double th = th_lo; th <= th_hi; th += 0.5) {
                slab::Window win; win.x1 = Sx; win.y1 = Sy;
                auto F = slab::dp(S, th * M_PI / 180, obs, win, step, tol, 12);
                const int tot = nup + F.count;
                std::lock_guard<std::mutex> lk(mu);
                if (tot > best) { best = tot; std::fprintf(stderr, "cfg %ld theta %.2f: uprights %d + columns %d = %d (target %d)\n", ci, th, nup, F.count, tot, w * h + dm); }
            }
        }
    };
    std::vector<std::thread> pool; for (int t = 0; t < 8; t++) pool.emplace_back(worker); for (auto& t : pool) t.join();
    std::printf("%dx%d+%d f=%.4f: best %d (target %d)\n", w, h, dm, f, best, w * h + dm);
}
