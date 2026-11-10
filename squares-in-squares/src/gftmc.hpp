// Gadget-level FTMC (free tilted Monte Carlo) + exact upright refill + rectangular squeeze.
//   Own implementation of the procedure described in Chris Deotte's "FTMC + fill_rect" notebook: for a gadget (w, h, dm) with
//   overhang f, mutate only the tilted core, relax it alone in the smaller box, refill the uprights exactly (src/fill.hpp),
//   and when the count reaches w*h+dm lock the result in with a squeeze of the rectangular box [0,w+f] x [0,h+f].
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <vector>

#include "fill.hpp"
#include "gadget.hpp"
#include "squeeze.hpp"
#include "window.hpp"

namespace gftmc {
using gadget::Donor;

inline void tips_of(const std::vector<Square>& sq, int w, int h, double* ybr, double* xtl) {
    *ybr = 9; *xtl = 9;
    for (const Square& q : sq) {
        const double he = 0.5 * (std::fabs(std::cos(q.angle)) + std::fabs(std::sin(q.angle)));
        if (q.x + he > w + 1e-9) *ybr = std::min(*ybr, q.y - he);
        if (q.y + he > h + 1e-9) *xtl = std::min(*xtl, q.x - he);
    }
}

// ---- gadget-level FTMC: mutate the tilted core, relax it alone, refill uprights exactly, lock in with a rectangular squeeze ----
inline bool is_tilted(const Square& q) { double a = std::fmod(std::fabs(q.angle), M_PI / 2); a = std::min(a, M_PI / 2 - a); return a > 0.02; }

inline bool relax_core(std::vector<Square>& T, double Sx, double yoff) {
    window::Env env; env.S = Sx + 2 * squeeze::kMargin; env.yoff = yoff; env.index();
    std::vector<double> v;
    for (const Square& q : T) { v.push_back(q.x); v.push_back(q.y); v.push_back(q.angle); }
    const double E = window::relax(env, v, 400);
    for (size_t i = 0; i < T.size(); i++) { T[i].x = v[3 * i]; T[i].y = v[3 * i + 1]; T[i].angle = v[3 * i + 2]; }
    return E <= 1e-14;
}

inline void mutate(std::vector<Square>& T, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> U(0, 1);
    std::normal_distribution<double> N01(0, 1);
    const int m = (int)T.size();
    const int mode = (int)(U(rng) * 8);
    if (m == 0) return;
    if (mode == 0 || mode == 1) {
        const int c = (int)(U(rng) * m);
        const double r = 1.0 + 3.0 * U(rng), sp = U(rng) < 0.5 ? 0.02 : 0.1, sa = U(rng) < 0.5 ? 0.02 : 0.15;
        const Square ref = T[c];
        for (auto& q : T) { const double dx = q.x - ref.x, dy = q.y - ref.y; if (dx * dx + dy * dy < r * r) { q.x += sp * N01(rng); q.y += sp * N01(rng); q.angle += sa * N01(rng); } }
    } else if (mode == 2) {
        const double d = 0.01 * N01(rng);
        for (auto& q : T) q.angle += d;
    } else if (mode == 3) {
        double cx = 0, cy = 0; for (auto& q : T) { cx += q.x; cy += q.y; } cx /= m; cy /= m;
        double sxx = 0, sxy = 0, syy = 0; for (auto& q : T) { sxx += (q.x - cx) * (q.x - cx); sxy += (q.x - cx) * (q.y - cy); syy += (q.y - cy) * (q.y - cy); }
        const double a = 0.5 * std::atan2(2 * sxy, sxx - syy), ux = std::cos(a), uy = std::sin(a);
        const double cut = (U(rng) - 0.5) * 6.0, ds = 0.3 * N01(rng);
        for (auto& q : T) if ((q.x - cx) * ux + (q.y - cy) * uy > cut) { q.x += ds * ux; q.y += ds * uy; }
    } else if (mode == 4) {
        const int c = (int)(U(rng) * m);
        T.push_back({T[c].x + (U(rng) - 0.5) * 2, T[c].y + (U(rng) - 0.5) * 2, T[c].angle + 0.05 * N01(rng)});
    } else if (mode == 5) {
        T.erase(T.begin() + (int)(U(rng) * m));
    } else if (mode == 6) {
        const int a = (int)(U(rng) * m), b = (int)(U(rng) * m);
        T[a].x = T[b].x + (U(rng) - 0.5) * 2; T[a].y = T[b].y + (U(rng) - 0.5) * 2; T[a].angle = T[b].angle + 0.05 * N01(rng);
    } else {   // re-angle the whole core to a rational slope (angle of the lattice vector (q, p))
        static const int pq[][2] = {{1, 2}, {2, 3}, {3, 4}, {1, 3}, {2, 5}, {1, 1}, {3, 5}, {1, 4}};
        const auto& r = pq[rng() % 8];
        const double th = std::atan2((double)r[0], (double)r[1]);
        const double old = T[0].angle;
        for (auto& q : T) q.angle = th + (q.angle - old) * 0.3 + 0.01 * N01(rng);
    }
}

struct CoreSeed { std::vector<Square> T; double W, H; };   // a tilted core with its box [0,W] x [0,H]

inline bool ftmc_gadget(const gadget::Gadget& g, const Donor& d, double seconds, uint64_t seed, gadget::Gadget& out, const std::vector<CoreSeed>* alts = nullptr, double goal_f = -1) {
    std::vector<Square> sq0;
    double f0 = g.f;
    if (!g.sq.empty()) sq0 = g.sq; else gadget::crop(d, g.sym, g.w, g.h, g.dm, sq0, &f0);
    const int target = g.w * g.h + g.dm;
    const double yoff = g.h - g.w;
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0, 1);
    std::vector<Square> best_sq = sq0;
    double best_f = f0;
    bool improved = false;
    const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    auto tilted_of = [&](const std::vector<Square>& v) { std::vector<Square> T; for (const Square& q : v) if (is_tilted(q)) T.push_back(q); return T; };
    std::vector<Square> T0 = tilted_of(best_sq);
    if (T0.empty()) return false;
    while (std::chrono::steady_clock::now() < t_end) {
        static const double deltas[] = {0.003, 0.008, 0.015, 0.03, 0.06, 0.1, 0.15};
        double f_try = best_f - deltas[rng() % 7];
        if (goal_f > 0 && goal_f < best_f - 1e-4 && U(rng) < 0.5) f_try = goal_f;   // aim straight at the f that would complete a chain
        if (f_try < 0.2) continue;
        const double Sx = g.w + f_try, Sy = g.h + f_try;
        std::vector<Square> cur = T0;
        double bw = g.w + best_f, bh = g.h + best_f;
        if (alts && !alts->empty() && U(rng) < 0.6) { const CoreSeed& cs = (*alts)[rng() % alts->size()]; cur = cs.T; bw = cs.W; bh = cs.H; }
        for (auto& q : cur) { q.x *= Sx / bw; q.y *= Sy / bh; }
        if (!relax_core(cur, Sx, yoff)) continue;
        auto total = [&](const std::vector<Square>& T, fill::Result* R) {
            fill::Result r = fill::upright_fill(Sx, Sy, T, 1e-9, 10, rng);
            const int tot = (int)T.size() + (int)r.up.size();
            if (R) *R = std::move(r);
            return tot;
        };
        int cur_total = total(cur, nullptr);
        int ep_best = cur_total; const int ep_start = cur_total;
        for (int step = 0; step < 250 && std::chrono::steady_clock::now() < t_end; step++) {
            std::vector<Square> T = cur;
            mutate(T, rng);
            if (!relax_core(T, Sx, yoff)) continue;
            fill::Result R;
            const int tot = total(T, &R);
            if (tot >= cur_total || (tot == cur_total - 1 && U(rng) < 0.02)) { cur = T; cur_total = tot; }
            ep_best = std::max(ep_best, tot);
            if (tot < target) continue;
            // candidate: verify with a full relax in the rectangular box
            std::vector<double> v;
            for (const Square& q : T) { v.push_back(q.x); v.push_back(q.y); v.push_back(q.angle); }
            for (const auto& c : R.up) { v.push_back(c.x); v.push_back(c.y); v.push_back(0.0); }
            while ((int)v.size() / 3 > target) v.resize(v.size() - 3);
            squeeze::State st; st.yoff = yoff;
            bool ok = false;
            for (double slack = 1e-9; slack <= 3e-4 && !ok; slack *= 8) {
                squeeze::State t = st; t.v = v; t.S = Sx + 2 * squeeze::kMargin + slack;
                if (squeeze::relax(t.v, t.S, 5000, yoff) <= squeeze::kTol) { st = t; ok = true; }
            }
            if (!ok) continue;
            squeeze::Options o; o.seconds = 3; o.threads = 1; o.seed = seed + step; o.max_iter = 250;
            squeeze::State r2 = squeeze::run_state(st, o);
            const double f_new = r2.S - g.w;
            if (std::getenv("GFTMC_DEBUG")) std::fprintf(stderr, "    candidate at f_try %.4f -> squeezed f %.5f (best %.5f)\n", f_try, f_new, best_f);
            if (f_new < best_f - 1e-6) {
                best_f = f_new; best_sq.clear();
                for (int i = 0; i < squeeze::count(r2); i++) best_sq.push_back({r2.v[3 * i], r2.v[3 * i + 1], r2.v[3 * i + 2]});
                T0 = tilted_of(best_sq);
                improved = true;
                break;
            }
        }
        if (std::getenv("GFTMC_DEBUG")) std::fprintf(stderr, "  episode f_try %.4f core %zu: start total %d, best %d, target %d\n", f_try, cur.size(), ep_start, ep_best, target);
    }
    if (!improved) return false;
    out = g; out.f = best_f; out.sq = best_sq;
    tips_of(out.sq, g.w, g.h, &out.ybr, &out.xtl);
    out.qb = gadget::bin_of(out.ybr); out.qt = gadget::bin_of(out.xtl);
    return true;
}


}  // namespace gftmc
