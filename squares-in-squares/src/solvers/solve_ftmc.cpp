// Approach: free tilted Monte Carlo (FTMC) on the tilted core, exact upright refill (src/fill.hpp), then full relax + polish.
// RESULT (2026-10-27): whole-puzzle FTMC fails (n=299 stays at 292/299; from a DP seed n=300 reaches only 293). Works only at gadget scale with neighbour-core seeds: see solve_gquench.cpp / gftmc.hpp.
//   Idea:    after Chris Deotte's "FTMC + fill_rect" briefing: perturbing a whole packing jams on the upright lattice's zero modes, so
//            keep only the tilted squares T as variables.  For a smaller box, mutate T (jitter / slide / re-angle / insert / delete),
//            relax T alone, and recompute the best upright fill (max independent set of wall-lattice + flush candidates).  A state
//            with |T| + |fill| >= n is a packing of n squares in the smaller box.
//   Run:     FTMC_TARGETS=299 FTMC_FRAC=1.0 FTMC_LEADER=analysis/data/per_puzzle_leaderboard.csv FTMC_SECONDS=300 make run SOLVER=ftmc
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <random>
#include <sstream>
#include <thread>

#include "../fault.hpp"
#include "../fill.hpp"
#include "../recompose.hpp"
#include "../solve.hpp"
#include "../squeeze.hpp"
#include "../window.hpp"

namespace {
bool wanted(int n) {
    const char* e = std::getenv("FTMC_TARGETS");
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
double leader_L(int n) {
    const char* path = std::getenv("FTMC_LEADER");
    if (!path) return -1;
    std::FILE* f = std::fopen(path, "r");
    if (!f) return -1;
    char line[512]; double out = -1;
    while (std::fgets(line, sizeof line, f)) { int k = 0; double L = 0; if (std::sscanf(line, "%d,%lf", &k, &L) == 2 && k == n) { out = L; break; } }
    std::fclose(f);
    return out;
}

using Core = std::vector<Square>;

// relax the tilted core alone inside [0,S]^2; returns true if (nearly) overlap-free
bool relax_core(Core& T, double S) {
    window::Env env; env.S = S; env.index();
    std::vector<double> v;
    for (const Square& q : T) { v.push_back(q.x); v.push_back(q.y); v.push_back(q.angle); }
    const double E = window::relax(env, v, 400);
    for (size_t i = 0; i < T.size(); i++) { T[i].x = v[3 * i]; T[i].y = v[3 * i + 1]; T[i].angle = v[3 * i + 2]; }
    return E <= 1e-14;
}

int total_of(const Core& T, double S, int iters, std::mt19937_64& rng, fill::Result* out = nullptr) {
    fill::Result R = fill::upright_fill(S, S, T, 1e-9, iters, rng);
    const int tot = (int)T.size() + (int)R.up.size();
    if (out) *out = std::move(R);
    return tot;
}
}  // namespace

Packing solve(int n, const Packing& base) {
    if (base.empty() || !wanted(n)) return base;
    const double L0 = container_side(base), seconds = env_d("FTMC_SECONDS", 120), lead = leader_L(n);
    double S = L0 - env_d("FTMC_DROP", 0.01);
    if (lead > 0) S = L0 - env_d("FTMC_FRAC", 1.0) * (L0 - lead);
    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    const recompose::Split sp = recompose::split(base);
    Core T0;
    for (const auto& t : sp.tilt) T0.push_back({t.x * S / L0, t.y * S / L0, t.t});
    if (std::getenv("FTMC_SEED_DP")) {   // benchmark mode: ignore the base core, start from the best column-DP strip at this side
        T0.clear(); int bestc = -1; fault::Result bR; fault::Params bP;
        for (double th = 14; th <= 42; th += 0.5) for (int a = 0; a <= 3; a++) for (int r = 1; r <= 3; r++) {
            fault::Params P; P.S = S; P.theta = th * M_PI / 180; P.a = a; P.b = a - 1 - r; P.tol = 0.02; P.max_h = 4;
            auto R = fault::solve(P, true); if (R.total() > bestc) { bestc = R.total(); bR = R; bP = P; }
        }
        const double ct = std::cos(bP.theta), st = std::sin(bP.theta);
        for (auto& c : bR.cols) for (int q = 0; q < c.h; q++) { const double w = c.w + (q - (c.h - 1) * 0.5); T0.push_back({c.u * ct - w * st, c.u * st + w * ct, bP.theta}); }
        std::fprintf(stderr, "  DP seed: theta %.1f a=%d, %zu tilted, DP total %d\n", bP.theta * 180 / M_PI, bP.a, T0.size(), bestc);
    }
    std::fprintf(stderr, "n=%d: FTMC target side %.6f (base %.6f), core %zu tilted squares, %d threads\n", n, S, L0, T0.size(), threads);
    {   // starting state
        Core T = T0; std::mt19937_64 rng(1);
        const bool ok = relax_core(T, S);
        std::fprintf(stderr, "  start: core relaxes %s, total %d (need %d)\n", ok ? "yes" : "NO", total_of(T, S, 40, rng), n);
    }
    std::mutex mu;
    std::atomic<bool> done{false};
    Core winner; int best_total = 0;
    const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    auto worker = [&](int tid) {
        std::mt19937_64 rng(777 + 131 * tid + n);
        std::uniform_real_distribution<double> U(0, 1);
        std::normal_distribution<double> N01(0, 1);
        Core cur = T0;
        relax_core(cur, S);
        int cur_total = total_of(cur, S, 20, rng);
        long it = 0, accepted = 0;
        while (!done && std::chrono::steady_clock::now() < t_end) {
            Core T = cur;
            const int mode = (int)(U(rng) * 6);
            const int m = (int)T.size();
            if (m == 0) break;
            if (mode == 0) {            // jitter a neighbourhood
                const int c = (int)(U(rng) * m);
                const double r = 1.0 + 3.0 * U(rng), sp_ = U(rng) < 0.5 ? 0.02 : 0.1, sa = U(rng) < 0.5 ? 0.02 : 0.15;
                for (auto& q : T) { const double dx = q.x - cur[c].x, dy = q.y - cur[c].y; if (dx * dx + dy * dy < r * r) { q.x += sp_ * N01(rng); q.y += sp_ * N01(rng); q.angle += sa * N01(rng); } }
            } else if (mode == 1) {     // rotate all tilt angles a little
                const double d = 0.01 * N01(rng);
                for (auto& q : T) q.angle += d;
            } else if (mode == 2) {     // slide a random half of the core along the core axis
                double cx = 0, cy = 0; for (auto& q : T) { cx += q.x; cy += q.y; } cx /= m; cy /= m;
                double sxx = 0, sxy = 0, syy = 0; for (auto& q : T) { sxx += (q.x - cx) * (q.x - cx); sxy += (q.x - cx) * (q.y - cy); syy += (q.y - cy) * (q.y - cy); }
                const double a = 0.5 * std::atan2(2 * sxy, sxx - syy), ux = std::cos(a), uy = std::sin(a);
                const double cut = (U(rng) - 0.5) * 2.0, ds = 0.3 * N01(rng);
                for (auto& q : T) if ((q.x - cx) * ux + (q.y - cy) * uy > cut * 4) { q.x += ds * ux; q.y += ds * uy; }
            } else if (mode == 3) {     // insert a tilted square near a random core square
                const int c = (int)(U(rng) * m);
                T.push_back({T[c].x + (U(rng) - 0.5) * 2, T[c].y + (U(rng) - 0.5) * 2, T[c].angle + 0.05 * N01(rng)});
            } else if (mode == 4) {     // delete one
                T.erase(T.begin() + (int)(U(rng) * m));
            } else {                    // move one square to the neighbourhood of another
                const int a = (int)(U(rng) * m), b = (int)(U(rng) * m);
                T[a].x = T[b].x + (U(rng) - 0.5) * 2; T[a].y = T[b].y + (U(rng) - 0.5) * 2; T[a].angle = T[b].angle + 0.05 * N01(rng);
            }
            it++;
            if (!relax_core(T, S)) continue;
            const int tot = total_of(T, S, 10, rng);
            if (tot >= cur_total || (tot == cur_total - 1 && U(rng) < 0.01)) {
                cur = T; cur_total = tot; accepted++;
                std::lock_guard<std::mutex> lk(mu);
                if (tot > best_total) { best_total = tot; std::fprintf(stderr, "  [t%d] it %ld: total %d / %d (core %zu)\n", tid, it, tot, n, T.size()); }
                if (tot >= n) { done = true; winner = T; }
            }
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker, t);
    for (auto& t : pool) t.join();
    if (!done) { std::fprintf(stderr, "n=%d: FTMC did not reach %d (best %d)\n", n, n, best_total); return base; }
    std::mt19937_64 rng(5);
    fill::Result R;
    total_of(winner, S, 300, rng, &R);
    std::vector<double> v;
    for (const Square& q : winner) { v.push_back(q.x); v.push_back(q.y); v.push_back(q.angle); }
    for (const auto& c : R.up) { v.push_back(c.x); v.push_back(c.y); v.push_back(0.0); }
    while ((int)v.size() / 3 > n) v.resize(v.size() - 3);
    double Sa = S + 2 * squeeze::kMargin + 1e-9;
    bool good = false;
    for (double slack = 0; slack <= 4e-3 && !good; slack = slack == 0 ? 1e-9 : slack * 4) {
        std::vector<double> t2 = v;
        if (squeeze::relax(t2, Sa + slack, 6000) <= squeeze::kTol) { v = t2; Sa += slack; good = true; }
    }
    if (!good) { std::fprintf(stderr, "n=%d: found a count-feasible state but the full relax failed\n", n); return base; }
    squeeze::polish(v, Sa);
    squeeze::State st; st.v = v; st.S = Sa;
    Packing out = squeeze::to_packing(st);
    std::fprintf(stderr, "n=%d: FTMC success, side %.9f (was %.9f)\n", n, container_side(out), L0);
    return out;
}
