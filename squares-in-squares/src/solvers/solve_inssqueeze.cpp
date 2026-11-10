// Approach: insert one square into the (n-1)-packing, enlarge the box until it relaxes, then squeeze back down.
// RESULT (2026-10-27): no gain. Random insertion into the (n-1) packing does not relax even with +0.12 box growth (n=299 from the 45-degree diamond 298).
//   Idea:    L(n) is often only ~0.02 above L(n-1) at the leaders (e.g. 299 vs our 45-degree-family 298 at f=.778), but the
//            (n-1) packing has no room for a new square without re-nesting everything.  Random insertion points, a box
//            enlarged by a few hundredths, penalty relax, then the shrinking-box hop search (single-thread per trial).
//   Run:     INS_TARGETS=299 INS_SRC=submissions/x.csv INS_SECONDS=30 make run SOLVER=inssqueeze
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <random>
#include <sstream>
#include <thread>

#include "../book.hpp"
#include "../solve.hpp"
#include "../squeeze.hpp"

namespace {
bool wanted(int n) {
    const char* e = std::getenv("INS_TARGETS");
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
    const char* src = std::getenv("INS_SRC");
    if (!src) { std::fprintf(stderr, "set INS_SRC to a submission csv\n"); return base; }
    const Packing down = book_packing(src, n - 1);
    if (down.empty()) return base;
    const double seconds = env_d("INS_SECONDS", 30), trials = env_d("INS_TRIALS", 64);
    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    const double L0 = container_side(base), Ld = container_side(down);
    std::fprintf(stderr, "n=%d: insert-and-squeeze from the (n-1) packing (side %.9f), base %.9f\n", n, Ld, L0);
    std::atomic<int> next{0};
    std::mutex mu;
    Packing best = base;
    double best_L = L0;
    auto worker = [&](int tid) {
        std::mt19937_64 rng(31 * tid + n);
        std::uniform_real_distribution<double> U(0, 1);
        for (;;) {
            const int tr = next++;
            if (tr >= trials) return;
            const double slack = 0.03 + 0.12 * U(rng);          // how much bigger than the base side we start from
            const double S = std::min(L0 + 0.02, Ld + slack);
            squeeze::State st = squeeze::from_packing(down);
            st.S = S + 2 * squeeze::kMargin;
            // new square at a random place, random tilt
            st.v.push_back(0.5 + (S - 1) * U(rng)); st.v.push_back(0.5 + (S - 1) * U(rng)); st.v.push_back(U(rng) * M_PI / 2);
            const double E = squeeze::relax(st.v, st.S, 6000);
            if (E > squeeze::kTol) { std::lock_guard<std::mutex> lk(mu); std::fprintf(stderr, "  trial %d: S=%.4f did not relax (E=%.2e)\n", tr, S, E); continue; }
            Packing p = squeeze::to_packing(st);
            squeeze::Options o;
            o.seconds = seconds; o.threads = 1; o.seed = 1000 + tr;
            Packing q = squeeze::run(p, o);
            const double L = container_side(q);
            std::lock_guard<std::mutex> lk(mu);
            std::fprintf(stderr, "  trial %d: start S=%.4f -> side %.9f\n", tr, S, L);
            if (L < best_L) { best_L = L; best = q; }
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker, t);
    for (auto& t : pool) t.join();
    std::fprintf(stderr, "n=%d: -> side %.9f (was %.9f)\n", n, best_L, L0);
    return best;
}
