// Approach: delete one square from the (n+1)-packing, then squeeze the hole shut.
// RESULT (2026-10-27): no gain. Deleting each tilted square of the (n+1) packing and squeezing never lowered L (n=299 from 300, all 43 tilted squares).
//   Idea:    the leaders' L(n) often sits ~0.02 below L(n+1) (e.g. n=299 vs 300, 335 vs 336): the n+1 packing has no slack,
//            but once a square is removed the strip can re-nest.  One thread per deleted square (single-thread squeeze).
//   Run:     DELSQ_TARGETS=299 DELSQ_SECONDS=20 DELSQ_SRC=submissions/x.csv make run SOLVER=delsqueeze
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
    const char* e = std::getenv("DELSQ_TARGETS");
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
    const char* src = std::getenv("DELSQ_SRC");
    if (!src) { std::fprintf(stderr, "set DELSQ_SRC to a submission csv\n"); return base; }
    const Packing up = book_packing(src, n + 1);
    if (up.empty()) return base;
    const double seconds = env_d("DELSQ_SECONDS", 20);
    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    const double L0 = container_side(base);
    std::fprintf(stderr, "n=%d: delete-and-squeeze from the (n+1) packing (side %.9f), base %.9f, %d candidates\n", n, container_side(up), L0, n + 1);
    std::atomic<int> next{0};
    std::mutex mu;
    Packing best = base;
    double best_L = L0;
    auto worker = [&]() {
        for (;;) {
            const int d = next++;
            if (d >= n + 1) return;
            if (std::getenv("DELSQ_TILTED")) { const double a = std::fmod(std::fabs(up[d].angle), M_PI / 2); if (std::min(a, M_PI / 2 - a) < 0.02) continue; }
            Packing p = up;
            p.erase(p.begin() + d);
            squeeze::Options o;
            o.seconds = seconds; o.threads = 1; o.seed = 7 + d;
            Packing q = squeeze::run(p, o);
            const double L = container_side(q);
            std::lock_guard<std::mutex> lk(mu);
            std::fprintf(stderr, "  delete %3d (angle %.1f deg): side %.9f\n", d, up[d].angle * 180 / M_PI, L);
            if (L < best_L) { best_L = L; best = q; }
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    std::fprintf(stderr, "n=%d: -> side %.9f (was %.9f)\n", n, best_L, L0);
    return best;
}
