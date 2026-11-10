// Approach: window re-solve. Free the tilted cluster (+ the grid layer next to it), freeze the rest, shrink the box by
// RESULT (2026-10-27): window re-solve (free cluster, frozen shifted grid) is rigid: even drop 2e-5 is infeasible for strips; delete/extra modes failed at 299 (E stays ~1e-3). Useful only as an engine (window.hpp) for periodic cells and end-cap probes.
//   WIN_DROP and basin-hop only the free squares until the penalty energy reaches zero.
//   Idea:    the leaders' packings differ from ours by ~0.02 in L mostly inside the tilted strips/slabs; the grid halves just
//            slide with the walls.  Searching only the cluster is ~10x cheaper than the full-packing engine.
//   Run:     WIN_TARGETS=299 WIN_DROP=0.0201 WIN_SECONDS=300 WIN_RADIUS=2.2 make run SOLVER=window
#include <atomic>
#include <chrono>
#include <cstdio>
#include <mutex>
#include <random>
#include <cstdlib>
#include <sstream>
#include <thread>

#include "../book.hpp"
#include "../recompose.hpp"
#include "../slab.hpp"
#include "../solve.hpp"
#include "../window.hpp"

namespace {
bool wanted(int n) {
    const char* e = std::getenv("WIN_TARGETS");
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

static Packing solve_one(int n, const Packing& base, double drop);

static double leader_L(int n) {
    const char* path = std::getenv("WIN_LEADER");
    if (!path) return -1;
    std::FILE* f = std::fopen(path, "r");
    if (!f) return -1;
    char line[512];
    double out = -1;
    while (std::fgets(line, sizeof line, f)) {
        int k = 0; double L = 0;
        if (std::sscanf(line, "%d,%lf", &k, &L) == 2 && k == n) { out = L; break; }
    }
    std::fclose(f);
    return out;
}

Packing solve(int n, const Packing& base) {
    if (base.empty() || !wanted(n)) return base;
    const double L0 = container_side(base), lead = leader_L(n);
    if (lead > 0) {   // leader-driven: try growing fractions of the known gap, keep the last success
        const double gap = L0 - lead;
        if (gap < 1e-6) return base;
        Packing cur = base;
        std::stringstream ss(std::getenv("WIN_SCALES") ? std::getenv("WIN_SCALES") : "0.05,0.2,0.5,1.0");
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            const double drop = std::atof(tok.c_str()) * gap - (L0 - container_side(cur));
            if (drop <= 1e-9) continue;
            Packing nx = solve_one(n, cur, drop);
            if (container_side(nx) < container_side(cur) - 1e-9) { cur = nx; std::fprintf(stderr, "n=%d: scale %s ok -> %.9f (leader %.9f)\n", n, tok.c_str(), container_side(cur), lead); }
            else { std::fprintf(stderr, "n=%d: scale %s failed\n", n, tok.c_str()); break; }
        }
        return cur;
    }
    return solve_one(n, base, env_d("WIN_DROP", 0.01));
}

static Packing solve_one(int n, const Packing& base, double drop) {
    const double seconds = env_d("WIN_SECONDS", 120), radius = env_d("WIN_RADIUS", 2.2);
    const int extra = (int)env_d("WIN_EXTRA", 0);
    const bool del = std::getenv("WIN_DELETE") != nullptr;   // start from the (n+1) packing and delete one free square
    Packing start = base;
    if (del) {
        const char* src = std::getenv("WIN_SRC");
        if (!src) { std::fprintf(stderr, "WIN_DELETE needs WIN_SRC\n"); return base; }
        start = book_packing(src, n + 1);
        if (start.empty()) return base;
    }
    if (extra > 0) {   // start from the (n - extra) packing of WIN_SRC; the box is grown by -WIN_DROP
        const char* src = std::getenv("WIN_SRC");
        if (!src) { std::fprintf(stderr, "WIN_EXTRA needs WIN_SRC\n"); return base; }
        start = book_packing(src, n - extra);
        if (start.empty()) return base;
    }
    const recompose::Split sp = recompose::split(start);
    const double Snew = sp.L - drop;
    const std::vector<slab::Sq> moved = recompose::shrink(sp, drop);
    // which kept axis squares become free: within `radius` of any tilted square
    window::Env env;
    env.S = Snew;
    std::vector<double> free_v;
    int n_tilted = 0;
    auto shift_of = [&](double x, double y, double* sx, double* sy) {   // shift of the nearest kept axis square
        double best = 1e18; *sx = *sy = 0;
        for (size_t i = 0; i < sp.axis.size(); i++) {
            const double dx = sp.axis[i].x - x, dy = sp.axis[i].y - y, d = dx * dx + dy * dy;
            if (d < best) { best = d; *sx = moved[i].x - sp.axis[i].x; *sy = moved[i].y - sp.axis[i].y; }
        }
    };
    for (size_t i = 0; i < sp.axis.size(); i++) {
        bool near = false;
        for (const slab::Sq& t : sp.tilt) {
            const double dx = t.x - sp.axis[i].x, dy = t.y - sp.axis[i].y;
            if (dx * dx + dy * dy < radius * radius) { near = true; break; }
        }
        if (near) { free_v.push_back(moved[i].x); free_v.push_back(moved[i].y); free_v.push_back(moved[i].t); }
        else env.frozen.push_back({moved[i].x, moved[i].y, moved[i].t});
    }
    const int n_free_axis = (int)free_v.size() / 3;
    for (const slab::Sq& t : sp.tilt) {
        double sx, sy; shift_of(t.x, t.y, &sx, &sy);
        free_v.push_back(t.x + sx); free_v.push_back(t.y + sy); free_v.push_back(t.t);
        n_tilted++;
    }
    {   // extra squares go near random free squares
        std::mt19937_64 rng(4321 + n);
        std::uniform_real_distribution<double> U(0, 1);
        const int base_free = (int)free_v.size() / 3;
        for (int e = 0; e < extra && base_free > 0; e++) {
            const int q = (int)(rng() % base_free);
            free_v.push_back(free_v[3 * q] + (U(rng) - 0.5) * 3); free_v.push_back(free_v[3 * q + 1] + (U(rng) - 0.5) * 3); free_v.push_back(U(rng) * M_PI / 2);
        }
    }
    env.index();
    std::fprintf(stderr, "n=%d: window re-solve at S=%.6f (drop %.5f): %zu frozen, %d free grid + %d tilted\n", n, Snew, drop, env.frozen.size(), n_free_axis, n_tilted);
    if (del) {
        const int nf = (int)free_v.size() / 3;
        double best_all = 1e300;
        for (int d = 0; d < nf; d++) {
            std::vector<double> v = free_v;
            v.erase(v.begin() + 3 * d, v.begin() + 3 * d + 3);
            window::relax(env, v, 3000);
            window::HopOptions o;
            o.seconds = seconds; o.threads = (int)std::max(1u, std::thread::hardware_concurrency()); o.seed = 91 + d; o.max_iter = 250;
            std::vector<double> sol; double be = 0;
            const bool ok = window::hop(env, v, o, sol, &be);
            if (be < best_all) best_all = be;
            std::fprintf(stderr, "  delete free square %d/%d (angle %.1f): %s best E %.3e (overall best %.3e)\n", d, nf, std::fmod(free_v[3 * d + 2] * 180 / M_PI, 90.0), ok ? "FEASIBLE" : "no", be, best_all);
            if (!ok) continue;
            Packing p;
            for (const window::Fixed& f : env.frozen) p.push_back({f.x, f.y, f.t});
            for (size_t i = 0; i < sol.size() / 3; i++) p.push_back({sol[3 * i], sol[3 * i + 1], sol[3 * i + 2]});
            std::vector<double> allv;
            for (const Square& q : p) { allv.push_back(q.x); allv.push_back(q.y); allv.push_back(q.angle); }
            double S = Snew + 2 * squeeze::kMargin + 1e-9;
            squeeze::polish(allv, S);
            squeeze::State st; st.v = allv; st.S = S;
            Packing out = squeeze::to_packing(st);
            std::fprintf(stderr, "n=%d: FOUND via deletion, polished side %.9f (was %.9f)\n", n, container_side(out), sp.L);
            return out;
        }
        return base;
    }
    if (std::getenv("WIN_DPONLY")) {
        std::vector<slab::Sq> obs;
        for (const window::Fixed& f : env.frozen) obs.push_back({f.x, f.y, f.t});
        int best = 0; double bth = 0;
        for (double th = 5; th < 175; th += 0.5) {
            const int c = slab::dp(Snew, th * M_PI / 180, obs, slab::Window(), 0.02, 0.02, 14).count;
            if (c > best) { best = c; bth = th; }
        }
        std::fprintf(stderr, "  DP-only: best %d of %d needed at theta %.1f (radius %.1f)\n", best, (int)free_v.size() / 3, bth, radius);
        return base;
    }
    if (std::getenv("WIN_SEEDS")) {
        // random structured seeds: column-DP fill of the widened corridor (tilt scan), remaining squares dropped in at random
        const double th_lo = env_d("WIN_TH_LO", 18), th_hi = env_d("WIN_TH_HI", 34), seed_sec = env_d("WIN_SEED_SECONDS", 20);
        const int total = (int)free_v.size() / 3;
        std::atomic<bool> found{false};
        std::mutex mu;
        std::vector<double> sol;
        double gbest = 1e300;
        std::vector<slab::Sq> obs;
        for (const window::Fixed& f : env.frozen) obs.push_back({f.x, f.y, f.t});
        const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
        auto worker = [&](int tid) {
            std::mt19937_64 rng(1234 + 77 * tid + n);
            std::uniform_real_distribution<double> U(0, 1);
            long seeds = 0;
            while (!found && std::chrono::steady_clock::now() < t_end) {
                const double th = (th_lo + (th_hi - th_lo) * U(rng)) * M_PI / 180;
                slab::Fill F = slab::dp(Snew, th, obs, slab::Window(), 0.02, 0.02, 12, true);
                std::vector<double> v;
                for (const slab::Sq& q : F.squares) { v.push_back(q.x); v.push_back(q.y); v.push_back(q.t); }
                // drop DP squares at random if it found more than needed
                while ((int)v.size() / 3 > total) { const size_t k = rng() % (v.size() / 3); v.erase(v.begin() + 3 * k, v.begin() + 3 * k + 3); }
                const int have = (int)v.size() / 3;
                for (int i = have; i < total; i++) {
                    const slab::Sq& a = free_v.size() ? slab::Sq{free_v[3 * (rng() % total)], free_v[3 * (rng() % total) + 1], 0} : slab::Sq{0, 0, 0};
                    v.push_back(a.x + (U(rng) - 0.5) * 1.0); v.push_back(a.y + (U(rng) - 0.5) * 1.0); v.push_back(U(rng) * M_PI / 2);
                }
                window::HopOptions o;
                o.seconds = seed_sec; o.threads = 1; o.seed = rng(); o.max_iter = 250;
                std::vector<double> out; double be = 0;
                const bool ok = window::hop(env, v, o, out, &be);
                std::lock_guard<std::mutex> lk(mu);
                seeds++;
                if (be < gbest) { gbest = be; std::fprintf(stderr, "  [t%d] seed %ld theta %.1f dp=%d/%d best E %.3e\n", tid, seeds, th * 180 / M_PI, have, total, be); }
                if (ok) { found = true; sol = out; }
            }
        };
        std::vector<std::thread> pool;
        for (int t = 0; t < (int)std::max(1u, std::thread::hardware_concurrency()); t++) pool.emplace_back(worker, t);
        for (auto& t : pool) t.join();
        if (!found) { std::fprintf(stderr, "n=%d: seeds found no zero-energy state (best E %.3e)\n", n, gbest); return base; }
        Packing p;
        for (const window::Fixed& f : env.frozen) p.push_back({f.x, f.y, f.t});
        for (size_t i = 0; i < sol.size() / 3; i++) p.push_back({sol[3 * i], sol[3 * i + 1], sol[3 * i + 2]});
        std::vector<double> allv;
        for (const Square& q : p) { allv.push_back(q.x); allv.push_back(q.y); allv.push_back(q.angle); }
        double S = Snew + 2 * squeeze::kMargin + 1e-9;
        squeeze::polish(allv, S);
        squeeze::State st; st.v = allv; st.S = S;
        Packing out = squeeze::to_packing(st);
        std::fprintf(stderr, "n=%d: FOUND, polished side %.9f (was %.9f)\n", n, container_side(out), sp.L);
        return out;
    }
    std::vector<double> v = free_v;
    const double E0 = window::relax(env, v, 3000);
    std::fprintf(stderr, "  initial relaxed energy %.3e\n", E0);
    window::HopOptions o;
    o.seconds = seconds; o.threads = (int)std::max(1u, std::thread::hardware_concurrency()); o.seed = 17 + n; o.verbose = true;
    std::vector<double> sol;
    double bestE = 0;
    if (!window::hop(env, v, o, sol, &bestE)) {
        std::fprintf(stderr, "n=%d: no zero-energy state found (best E=%.3e)\n", n, bestE);
        return base;
    }
    Packing p;
    for (const window::Fixed& f : env.frozen) p.push_back({f.x, f.y, f.t});
    for (size_t i = 0; i < sol.size() / 3; i++) p.push_back({sol[3 * i], sol[3 * i + 1], sol[3 * i + 2]});
    std::fprintf(stderr, "n=%d: FOUND packing at side %.9f (was %.9f)\n", n, container_side(p), sp.L);
    std::vector<double> allv;
    for (const Square& q : p) { allv.push_back(q.x); allv.push_back(q.y); allv.push_back(q.angle); }
    double S = Snew + 2 * squeeze::kMargin + 1e-9;
    squeeze::polish(allv, S);
    squeeze::State st; st.v = allv; st.S = S;
    Packing out = squeeze::to_packing(st);
    std::fprintf(stderr, "n=%d: polished side %.9f\n", n, container_side(out));
    return out;
}
