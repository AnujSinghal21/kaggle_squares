// Approach: diagonal-strip generator (src/strip.hpp) + L-BFGS relax + smooth polish.
// RESULT (2026-10-27): dead end. Random lattice-band generator reaches counts far below the known packing (291 vs 299 at n=299); superseded by solve_fault.cpp / the motif lab.
//   Idea:    best-known packings for n > ~300 are two corner-anchored grids with a band of tilted squares between
//            them.  Sample the band parameters (tilt, phase, offset, width) at many container sides S below the
//            base side, keep configurations whose square count reaches n, and let the penetration-energy
//            L-BFGS close the small overlaps.  Anything that relaxes to zero energy is a valid packing at S.
//   Run:     STRIP_TARGETS=299,335 STRIP_SECONDS=60 make run SOLVER=strip
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <sstream>
#include <thread>

#include "../solve.hpp"
#include "../squeeze.hpp"
#include "../strip.hpp"

namespace {
bool wanted(int n) {
    const char* e = std::getenv("STRIP_TARGETS");
    if (!e) return false;
    std::stringstream ss(e);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        const size_t dash = tok.find('-');
        if (dash == std::string::npos) { if (std::atoi(tok.c_str()) == n) return true; }
        else if (n >= std::atoi(tok.substr(0, dash).c_str()) && n <= std::atoi(tok.substr(dash + 1).c_str())) return true;
    }
    return false;
}
double env_d(const char* k, double d) { const char* e = std::getenv(k); return e ? std::atof(e) : d; }
}  // namespace

Packing solve(int n, const Packing& base) {
    if (base.empty() || !wanted(n)) return base;
    const double seconds = env_d("STRIP_SECONDS", 60), depth = env_d("STRIP_DEPTH", 0.04);
    const int threads = (int)env_d("STRIP_THREADS", std::max(1u, std::thread::hardware_concurrency()));
    const double L0 = container_side(base);
    std::fprintf(stderr, "n=%d: strip search below %.9f for %.0fs on %d threads\n", n, L0, seconds, threads);

    std::mutex mu;
    double best_S = L0;
    std::vector<double> best_v;
    long tried = 0, counted = 0, relaxed = 0, ok = 0;
    int near_miss = -999;
    const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);

    auto worker = [&](int tid) {
        std::mt19937_64 rng(777ULL * n + tid * 104729ULL);
        std::uniform_real_distribution<double> U(0, 1);
        while (std::chrono::steady_clock::now() < t_end) {
            strip::Config c;
            { std::lock_guard<std::mutex> lk(mu); c.S = best_S - 1e-4 - U(rng) * depth; }
            c.theta = (8 + 37 * U(rng)) * M_PI / 180;
            c.px = U(rng); c.py = U(rng);
            c.s = (U(rng) - 0.5) * 4;
            c.h = 0.8 + 1.6 * U(rng);
            c.tol = 0.01 + 0.08 * U(rng);
            strip::Generated g = strip::generate(c);
            const int cnt = (int)g.v.size() / 3;
            long cn = 0;
            {
                std::lock_guard<std::mutex> lk(mu);
                tried++;
                if (cnt - n > near_miss && cnt - n <= 0) near_miss = cnt - n;
                if (cnt >= n) cn = ++counted;
            }
            if (cnt < n || cnt > n + 12) continue;
            std::vector<double> v = g.v;
            strip::trim_to(v, g.n_strip, n, rng);
            const double S = c.S + 2 * squeeze::kMargin + 1e-9;
            const double E = squeeze::relax(v, S, 2500);
            { std::lock_guard<std::mutex> lk(mu); relaxed++; }
            if (E > squeeze::kTol) continue;
            std::lock_guard<std::mutex> lk(mu);
            ok++;
            if (S < best_S) {
                best_S = S; best_v = v;
                std::fprintf(stderr, "  [t%d] n=%d feasible at S=%.9f (theta=%.2f s=%.2f h=%.2f cnt=%d)\n", tid, n, S, c.theta * 180 / M_PI, c.s, c.h, cnt);
            }
            (void)cn;
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    std::fprintf(stderr, "  tried %ld, count>=n %ld, relaxed %ld, feasible %ld, best count-n below target %d\n", tried, counted, relaxed, ok, near_miss);
    if (best_v.empty()) return base;
    double S = best_S;
    squeeze::polish(best_v, S);
    squeeze::State st; st.v = best_v; st.S = S;
    Packing p = squeeze::to_packing(st);
    std::fprintf(stderr, "n=%d: -> side %.9f (was %.9f)\n", n, container_side(p), L0);
    return p;
}
