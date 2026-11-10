// Systematic "core insertion": take the tilted core of gadget A (w,h,dm), stretch it into the box of the next key (w,h+1,dm+1)
// with overhang f, try every insertion point/angle on a grid, relax the core, refill the uprights exactly, report the best total.
//   usage: insert_probe dump_A.txt f_new [dh=1] [ddm=1]
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <atomic>
#include <mutex>
#include <thread>
#include "../gftmc.hpp"
int main(int argc, char** argv) {
    std::ifstream in(argv[1]); std::stringstream all; all << in.rdbuf(); std::string txt = all.str();
    for (size_t k; (k = txt.find("\\n")) != std::string::npos;) txt.replace(k, 2, "\n");
    std::stringstream hs(txt);
    int w, h, dm; double f0; hs >> w >> h >> dm >> f0;
    std::vector<Square> sq; double x, y, t;
    while (hs >> x >> y >> t) sq.push_back({x, y, t});
    const double fn = atof(argv[2]); const int dh = argc > 3 ? atoi(argv[3]) : 1, ddm = argc > 4 ? atoi(argv[4]) : 1;
    const double ftol = argc > 5 ? atof(argv[5]) : 1e-9;
    const int h2 = h + dh, target = w * h2 + dm + ddm;
    const double Sx = w + fn, Sy = h2 + fn, yoff = h2 - w;
    std::vector<Square> T0; for (auto& q : sq) if (gftmc::is_tilted(q)) T0.push_back(q);
    double ang = 0; for (auto& q : T0) ang += q.angle; ang /= std::max<size_t>(1, T0.size());
    std::printf("core %zu tilted (mean angle %.2f deg), new box %.4f x %.4f, target %d\n", T0.size(), ang * 180 / M_PI, Sx, Sy, target);
    // stretch variants: scale y by Sy/(h+f0) and x by Sx/(w+f0)
    std::vector<Square> base = T0; for (auto& q : base) { q.x *= Sx / (w + f0); q.y *= Sy / (h + f0); }
    std::vector<std::pair<double, double>> pts;
    for (double px = 0.5; px < Sx - 0.4; px += 0.2) for (double py = 0.5; py < Sy - 0.4; py += 0.2) pts.push_back({px, py});
    std::atomic<size_t> next{0}; std::mutex mu; int best = 0;
    auto worker = [&](int tid) {
        std::mt19937_64 rng(11 + tid);
        for (;;) {
            const size_t i = next++; if (i >= pts.size()) return;
            for (int ai = 0; ai < 5; ai++) {
                std::vector<Square> T = base;
                T.push_back({pts[i].first, pts[i].second, ang + (ai - 2) * 0.1});
                if (!gftmc::relax_core(T, Sx, yoff)) continue;
                auto R = fill::upright_fill(Sx, Sy, T, ftol, 20, rng);
                const int tot = (int)T.size() + (int)R.up.size();
                std::lock_guard<std::mutex> lk(mu);
                if (tot >= target) {   // loosen-then-squeeze: full relax of core+uprights at growing overhang
                    std::vector<double> v; for (auto& q : T) { v.push_back(q.x); v.push_back(q.y); v.push_back(q.angle); }
                    for (auto& c : R.up) { v.push_back(c.x); v.push_back(c.y); v.push_back(0.0); }
                    while ((int)v.size() / 3 > target) v.resize(v.size() - 3);
                    for (double extra = 0; extra <= 0.12; extra += 0.01) {
                        std::vector<double> vv = v;
                        if (squeeze::relax(vv, Sx + extra + 2 * squeeze::kMargin, 3000, yoff) <= squeeze::kTol) {
                            std::lock_guard<std::mutex> lk(mu);
                            std::printf("  FEASIBLE with count %d at f=%.4f (+%.2f over request) insert (%.2f,%.2f)\n", tot, fn + extra, extra, pts[i].first, pts[i].second); std::fflush(stdout); break;
                        }
                    }
                }
                if (tot > best) { best = tot; std::printf("  insert at (%.2f,%.2f) ang+%.1f: %zu + %zu = %d / %d\n", pts[i].first, pts[i].second, (ai - 2) * 0.1, T.size(), R.up.size(), tot, target); std::fflush(stdout); }
            }
        }
    };
    std::vector<std::thread> pool; for (int t = 0; t < 8; t++) pool.emplace_back(worker, t); for (auto& t : pool) t.join();
    std::printf("best %d (target %d)\n", best, target);
}
