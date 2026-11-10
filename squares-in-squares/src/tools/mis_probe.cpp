// Probe: max independent set count of {4 axis lattices + tilted lattice(s)} at container side S.
//   usage: mis_probe S seconds [tol] [thmin thmax] [nlat]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <thread>
#include "../mis.hpp"
int main(int argc, char** argv) {
    const double S = std::atof(argv[1]), secs = std::atof(argv[2]);
    const double tol = argc > 3 ? std::atof(argv[3]) : 0.02;
    const double thmin = argc > 4 ? std::atof(argv[4]) : 15, thmax = argc > 5 ? std::atof(argv[5]) : 35;
    const int nlat = argc > 6 ? std::atoi(argv[6]) : 1;
    const int threads = 8;
    std::mutex mu; int best = 0; long evals = 0; double bth = 0;
    const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(secs);
    auto worker = [&](int tid) {
        std::mt19937_64 rng(99 + tid); std::uniform_real_distribution<double> U(0, 1);
        while (std::chrono::steady_clock::now() < t_end) {
            std::vector<mis::Cand> c; mis::add_axis_lattices(c, S);
            double th0 = 0;
            for (int l = 0; l < nlat; l++) {
                const double th = (thmin + (thmax - thmin) * U(rng)) * M_PI / 180; if (l == 0) th0 = th;
                mis::add_lattice(c, S, th, U(rng) * 2 + S / 2 - 1, U(rng) * 2 + S / 2 - 1);
            }
            auto g = mis::build(std::move(c), tol);
            auto sol = mis::solve(g, 150, rng);
            std::lock_guard<std::mutex> lk(mu); evals++;
            if ((int)sol.size() > best) { best = sol.size(); bth = th0 * 180 / M_PI; std::fprintf(stderr, "[t%d] count %d theta %.2f (cands %zu)\n", tid, best, bth, g.c.size()); }
        }
    };
    std::vector<std::thread> pool; for (int t = 0; t < threads; t++) pool.emplace_back(worker, t);
    for (auto& t : pool) t.join();
    std::printf("S=%.5f best %d evals %ld\n", S, best, evals);
}
