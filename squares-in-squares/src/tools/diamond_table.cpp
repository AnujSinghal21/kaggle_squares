// For every container side S = k + f (f scanned), the max count of {45-degree W x H block on the diagonal + 4 corner lattices}.
// Output: lines "n S_min" = smallest side found reaching at least n squares.
// usage: diamond_table k_min k_max f_step theta_deg tol
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <thread>
#include "../block.hpp"
int main(int argc, char** argv) {
    const int k0 = atoi(argv[1]), k1 = atoi(argv[2]);
    const double fs = atof(argv[3]), th = atof(argv[4]) * M_PI / 180, tol = atof(argv[5]);
    std::map<int, double> best;  // count -> min S
    std::mutex mu;
    for (int k = k0; k <= k1; k++) {
        std::vector<double> Ss; for (double f = 0.0; f < 1.0; f += fs) Ss.push_back(k + f);
        std::atomic<size_t> next{0};
        auto worker = [&]() {
            for (;;) {
                size_t i = next++; if (i >= Ss.size()) return;
                const double S = Ss[i]; int bc = 0;
                const int wmax = (int)(S * 0.75);
                for (int W = 3; W <= wmax; W++) for (int H = W; H <= wmax && H <= W + 3; H++)
                    for (double c = S * 0.3; c < S * 0.7; c += 0.01) {
                        block::Params P; P.S = S; P.theta = th; P.cx = c; P.cy = c; P.W = W; P.H = H; P.tol = tol;
                        auto B = block::build(P);
                        if (B.n_block < 0) continue;
                        bc = std::max(bc, B.total());
                    }
                std::lock_guard<std::mutex> lk(mu);
                for (int c = bc; c > 0; c--) { auto it = best.find(c); if (it == best.end() || S < it->second) best[c] = S; else break; }
            }
        };
        std::vector<std::thread> pool; for (int t = 0; t < 8; t++) pool.emplace_back(worker); for (auto& t : pool) t.join();
        std::fprintf(stderr, "k=%d done\n", k);
    }
    for (auto& kv : best) printf("%d %.6f\n", kv.first, kv.second);
}
