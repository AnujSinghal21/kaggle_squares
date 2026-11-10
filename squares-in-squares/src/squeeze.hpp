// Squeeze engine: shrink the box around a packing of n unit squares.
//
//   energy  E(v, S) = sum over close pairs of P_ij^2  +  sum over squares of wall_violation^2
//   P_ij    = max(0, min over the 4 SAT axes of overlap)   (penetration depth, = Kaggle's SAT overlap)
//   solve   L-BFGS on E at a fixed box side S; E ~ 0 means a feasible packing with a small real gap.
//   search  shrink S, and when stuck kick a few squares (position + angle) and relax again
//           (basin hopping, Metropolis on E), from several threads that share the best result.
//
// For unit squares P_ij has a closed form (no corner loops): with d = c_j - c_i and phi = theta_j - theta_i,
//   H = 0.5 + 0.5(|cos phi| + |sin phi|),   overlap on axis k = H - |d . u_k|,
// where u_k are the four edge normals (two from each square).
#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#include "common.hpp"

namespace squeeze {

constexpr double kMargin = 3e-11;  // every pair/wall is inflated by this, so a converged state has real gaps
constexpr double kTol = 1e-24;     // E below this = feasible

struct State {
    std::vector<double> v;  // x0,y0,t0, x1,y1,t1, ...
    double S = 0;           // box width (and height - yoff)
    double yoff = 0;        // box is [0,S] x [0,S+yoff]; non-zero only for gadget boxes
};

inline int count(const State& s) { return (int)s.v.size() / 3; }

// E and gradient (g resized to v.size()). Pairs closer than sqrt(2) can overlap.
inline double energy(const std::vector<double>& v, double S, std::vector<double>& g, double* dES = nullptr, double yoff = 0) {
    const int n = (int)v.size() / 3;
    g.assign(v.size(), 0.0);
    double E = 0;
    if (dES) *dES = 0;
    if (!(S > 0.5) || !(S < 1e4) || !(std::fabs(yoff) < 1e4)) return 1e30;   // runaway box side during polish: reject the step
    // walls: AABB half extent h(theta) must fit inside [m, S-m]
    for (int i = 0; i < n; i++) {
        const double x = v[3 * i], y = v[3 * i + 1], t = v[3 * i + 2];
        const double c = std::cos(t), s = std::sin(t);
        const double h = 0.5 * (std::fabs(c) + std::fabs(s)) + kMargin;
        const double dh = 0.5 * (-s * (c >= 0 ? 1 : -1) + c * (s >= 0 ? 1 : -1));
        const double w[4] = {h - x, x + h - S, h - y, y + h - S - yoff};
        const double dx[4] = {-1, 1, 0, 0}, dy[4] = {0, 0, -1, 1};
        for (int k = 0; k < 4; k++)
            if (w[k] > 0) {
                E += w[k] * w[k];
                g[3 * i] += 2 * w[k] * dx[k];
                g[3 * i + 1] += 2 * w[k] * dy[k];
                g[3 * i + 2] += 2 * w[k] * dh;
                if (dES && (k == 1 || k == 3)) *dES -= 2 * w[k];
            }
    }
    // pairs via a uniform grid of cell size 1.5
    const double cell = 1.5;
    const int gx = (int)std::ceil((S + std::max(0.0, yoff) + 2) / cell) + 2;
    std::vector<int> head(gx * gx, -1), next(n, -1);
    auto cid = [&](double x, double y) {
        int cx = std::min(gx - 1, std::max(0, (int)((x + 1) / cell))), cy = std::min(gx - 1, std::max(0, (int)((y + 1) / cell)));
        return cy * gx + cx;
    };
    for (int i = 0; i < n; i++) { const int c = cid(v[3 * i], v[3 * i + 1]); next[i] = head[c]; head[c] = i; }
    const double reach2 = 2.0 + 4 * kMargin + 1e-6;
    for (int i = 0; i < n; i++) {
        const double xi = v[3 * i], yi = v[3 * i + 1], ti = v[3 * i + 2];
        const int cx = std::min(gx - 1, std::max(0, (int)((xi + 1) / cell))), cy = std::min(gx - 1, std::max(0, (int)((yi + 1) / cell)));
        for (int oy = -1; oy <= 1; oy++)
            for (int ox = -1; ox <= 1; ox++) {
                const int ax = cx + ox, ay = cy + oy;
                if (ax < 0 || ay < 0 || ax >= gx || ay >= gx) continue;
                for (int j = head[ay * gx + ax]; j >= 0; j = next[j]) {
                    if (j <= i) continue;
                    const double dx = v[3 * j] - xi, dy = v[3 * j + 1] - yi;
                    if (dx * dx + dy * dy > reach2) continue;
                    const double tj = v[3 * j + 2];
                    const double ci = std::cos(ti), si = std::sin(ti), cj = std::cos(tj), sj = std::sin(tj);
                    const double phi = tj - ti, cp = std::cos(phi), sp = std::sin(phi);
                    const double H = 0.5 + 0.5 * (std::fabs(cp) + std::fabs(sp)) + kMargin;
                    const double dH = 0.5 * (-sp * (cp >= 0 ? 1 : -1) + cp * (sp >= 0 ? 1 : -1));  // dH/dtj (= -dH/dti)
                    const double d[4] = {dx * ci + dy * si, -dx * si + dy * ci, dx * cj + dy * sj, -dx * sj + dy * cj};
                    int k = 0;
                    for (int q = 1; q < 4; q++) if (std::fabs(d[q]) > std::fabs(d[k])) k = q;
                    const double P = H - std::fabs(d[k]);
                    if (P <= 0) continue;
                    const double sg = d[k] >= 0 ? 1.0 : -1.0;
                    // dd[k]/d(dx,dy,ti,tj)
                    double ddx, ddy, ddti = 0, ddtj = 0;
                    if (k == 0)      { ddx = ci;  ddy = si;  ddti = d[1]; }
                    else if (k == 1) { ddx = -si; ddy = ci;  ddti = -d[0]; }
                    else if (k == 2) { ddx = cj;  ddy = sj;  ddtj = d[3]; }
                    else             { ddx = -sj; ddy = cj;  ddtj = -d[2]; }
                    E += P * P;
                    const double f = 2 * P;
                    const double gdx = f * (-sg * ddx), gdy = f * (-sg * ddy);
                    g[3 * j] += gdx;       g[3 * i] -= gdx;
                    g[3 * j + 1] += gdy;   g[3 * i + 1] -= gdy;
                    g[3 * i + 2] += f * (-dH - sg * ddti);
                    g[3 * j + 2] += f * (dH - sg * ddtj);
                }
            }
    }
    return E;
}

// Generic L-BFGS with backtracking line search on f(x, g) -> value. Stops when stop(value) is true,
// after max_iter iterations, or when no progress is possible. Returns the final value.
template <class F, class Stop>
double lbfgs(std::vector<double>& x, F f, int max_iter, Stop stop) {
    const int N = (int)x.size(), M = 8;
    std::vector<double> g, g_new, x_new(N), d(N);
    std::vector<std::vector<double>> sv, yv;
    std::vector<double> rho;
    double E = f(x, g);
    for (int it = 0; it < max_iter && !stop(E); it++) {
        d = g;
        const int h = (int)sv.size();
        std::vector<double> alpha(h);
        for (int i = h - 1; i >= 0; i--) {
            double a = 0; for (int k = 0; k < N; k++) a += sv[i][k] * d[k];
            alpha[i] = rho[i] * a;
            for (int k = 0; k < N; k++) d[k] -= alpha[i] * yv[i][k];
        }
        if (h > 0) {
            double sy = 0, yy = 0;
            for (int k = 0; k < N; k++) { sy += sv[h - 1][k] * yv[h - 1][k]; yy += yv[h - 1][k] * yv[h - 1][k]; }
            const double gam = sy / yy;
            for (int k = 0; k < N; k++) d[k] *= gam;
        } else {
            for (int k = 0; k < N; k++) d[k] *= 0.25;
        }
        for (int i = 0; i < h; i++) {
            double b = 0; for (int k = 0; k < N; k++) b += yv[i][k] * d[k];
            const double c = alpha[i] - rho[i] * b;
            for (int k = 0; k < N; k++) d[k] += c * sv[i][k];
        }
        for (int k = 0; k < N; k++) d[k] = -d[k];
        double gd = 0; for (int k = 0; k < N; k++) gd += g[k] * d[k];
        if (gd >= 0) { sv.clear(); yv.clear(); rho.clear(); for (int k = 0; k < N; k++) d[k] = -g[k] * 0.25; gd = 0; for (int k = 0; k < N; k++) gd += g[k] * d[k]; }
        double step = 1.0, E_new = E;
        bool ok = false;
        for (int ls = 0; ls < 40; ls++) {
            for (int k = 0; k < N; k++) x_new[k] = x[k] + step * d[k];
            E_new = f(x_new, g_new);
            if (E_new <= E + 1e-4 * step * gd) { ok = true; break; }
            step *= 0.5;
        }
        if (!ok) break;
        std::vector<double> s(N), y(N);
        double sy = 0;
        for (int k = 0; k < N; k++) { s[k] = x_new[k] - x[k]; y[k] = g_new[k] - g[k]; sy += s[k] * y[k]; }
        if (sy > 1e-30) {
            if ((int)sv.size() == M) { sv.erase(sv.begin()); yv.erase(yv.begin()); rho.erase(rho.begin()); }
            sv.push_back(s); yv.push_back(y); rho.push_back(1.0 / sy);
        }
        x = x_new; g = g_new;
        const double dE = E - E_new;
        E = E_new;
        if (std::fabs(dE) <= 1e-13 * std::fabs(E) && it > 30) break;
    }
    return E;
}

// Relax at a fixed box side; returns E (<= kTol means feasible).
inline double relax(std::vector<double>& v, double S, int max_iter, double yoff = 0) {
    return lbfgs(v, [&](const std::vector<double>& x, std::vector<double>& g) { return energy(x, S, g, nullptr, yoff); }, max_iter,
                 [](double E) { return E <= kTol; });
}

// Smooth end-game: minimise S + mu/2 * E(v, S) over (v, S) with growing mu, then find the smallest
// feasible side just above the result. Returns true and updates (v, S) if it ends feasible and smaller.
inline bool polish(std::vector<double>& v, double& S, double yoff = 0) {
    std::vector<double> x = v;
    x.push_back(S);
    const int N = (int)v.size();
    for (double mu : {1e3, 1e5, 1e7, 1e9, 1e11}) {
        lbfgs(x, [&](const std::vector<double>& y, std::vector<double>& g) {
            std::vector<double> vv(y.begin(), y.begin() + N), gv;
            double dES = 0;
            const double E = energy(vv, y[N], gv, &dES, yoff);
            g.assign(N + 1, 0.0);
            for (int k = 0; k < N; k++) g[k] = 0.5 * mu * gv[k];
            g[N] = 1.0 + 0.5 * mu * dES;
            return y[N] + 0.5 * mu * E;
        }, 1500, [](double) { return false; });
    }
    std::vector<double> w(x.begin(), x.begin() + N);
    for (double delta = 1e-10; delta < 1e-3; delta *= 2) {
        std::vector<double> t = w;
        if (relax(t, x[N] + delta, 3000, yoff) <= kTol) {
            if (x[N] + delta < S) { v = t; S = x[N] + delta; return true; }
            return false;
        }
    }
    return false;
}

inline State from_packing(const Packing& p) {
    State s;
    s.v.reserve(p.size() * 3);
    double lo_x = INFINITY, lo_y = INFINITY;
    for (const Square& q : p) {
        const double h = 0.5 * (std::fabs(std::cos(q.angle)) + std::fabs(std::sin(q.angle)));
        lo_x = std::min(lo_x, q.x - h); lo_y = std::min(lo_y, q.y - h);
    }
    for (const Square& q : p) { s.v.push_back(q.x - lo_x + kMargin); s.v.push_back(q.y - lo_y + kMargin); s.v.push_back(q.angle); }
    s.S = container_side(p) + 2 * kMargin + 1e-9;
    return s;
}

inline Packing to_packing(const State& s) {
    Packing p;
    for (int i = 0; i < count(s); i++) p.push_back({s.v[3 * i], s.v[3 * i + 1], s.v[3 * i + 2]});
    return p;
}

struct Options {
    double seconds = 30;
    int threads = 4;
    uint64_t seed = 1;
    int max_iter = 400;
    bool verbose = false;
    double fixed_drop = 0;  // > 0: skip gradual shrinking, hop at the fixed target side (start - fixed_drop)
};

// Shrinks the box of a feasible state (S and yoff fixed relation) for opt.seconds; returns the best state found.
inline State run_state(State best, const Options& opt) {
    const int n = count(best);
    std::mutex mu;
    const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(opt.seconds);
    auto worker = [&](int tid) {
        std::mt19937_64 rng(opt.seed * 1000003ULL + tid * 7919ULL);
        std::uniform_real_distribution<double> U(0, 1);
        std::normal_distribution<double> N01(0, 1);
        State mine;
        { std::lock_guard<std::mutex> lk(mu); mine = best; }
        double step = 2e-3;
        State cur;  // state at the current target side
        double cur_E = INFINITY, target = 0;
        int stuck = 0;
        while (std::chrono::steady_clock::now() < t_end) {
            {   // adopt a better global best
                std::lock_guard<std::mutex> lk(mu);
                if (best.S < mine.S - 1e-12) { mine = best; cur_E = INFINITY; }
            }
            if (cur_E == INFINITY) {   // new target: scale the best state down
                target = opt.fixed_drop > 0 ? mine.S - opt.fixed_drop : mine.S * (1 - step);
                cur.S = target; cur.yoff = mine.yoff;
                cur.v = mine.v;
                const double f = target / mine.S;
                for (int i = 0; i < n; i++) { cur.v[3 * i] *= f; cur.v[3 * i + 1] *= f; }
                cur_E = relax(cur.v, target, opt.max_iter, cur.yoff);
                stuck = 0;
            }
            if (cur_E <= kTol) {  // feasible at target
                mine = cur; cur_E = INFINITY;
                std::lock_guard<std::mutex> lk(mu);
                if (mine.S < best.S) { best = mine; if (opt.verbose) std::fprintf(stderr, "  [t%d] S=%.9f\n", tid, best.S); }
                step = std::min(step * 1.3, 5e-3);
                continue;
            }
            // kick a neighbourhood of a random square
            State trial = cur;
            const int c = (int)(U(rng) * n);
            const double r = 0.8 + 2.2 * U(rng), sig_p = (U(rng) < 0.5 ? 0.03 : 0.15), sig_a = (U(rng) < 0.5 ? 0.05 : 0.4);
            const bool whole = U(rng) < 0.02;
            for (int i = 0; i < n; i++) {
                const double dx = trial.v[3 * i] - cur.v[3 * c], dy = trial.v[3 * i + 1] - cur.v[3 * c + 1];
                if (!whole && dx * dx + dy * dy > r * r) continue;
                trial.v[3 * i] += sig_p * N01(rng);
                trial.v[3 * i + 1] += sig_p * N01(rng);
                trial.v[3 * i + 2] += sig_a * N01(rng);
            }
            const double E = relax(trial.v, target, opt.max_iter, cur.yoff);
            const double T = 0.02 * cur_E;  // Metropolis temperature relative to the current energy
            if (E < cur_E || U(rng) < std::exp(-(E - cur_E) / std::max(T, 1e-18))) { cur = trial; cur_E = E; stuck = 0; }
            else if (opt.fixed_drop > 0) {
                if (opt.verbose && tid == 0 && (++stuck % 200) == 0) std::fprintf(stderr, "  [t0] kicks without success: %d  cur_E=%.3e\n", stuck, cur_E);
            }
            else if (++stuck > 60) {  // give up this target: polish the best state smoothly, then try a smaller step
                step = std::max(step * 0.5, 1e-7);
                cur_E = INFINITY;
                if (step <= 1e-6) {
                    State q = mine;
                    if (polish(q.v, q.S, q.yoff)) { mine = q; std::lock_guard<std::mutex> lk(mu); if (mine.S < best.S) best = mine; }
                }
            }
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < opt.threads; t++) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    polish(best.v, best.S, best.yoff);
    return best;
}

// Shrinks the square box around `start` (must be feasible at its own side) for `seconds`. Returns the best packing found.
inline Packing run(const Packing& start, const Options& opt) {
    State best = from_packing(start);
    const double start_side = container_side(start);
    {   // make the start strictly feasible under our margins: the base only has ~1e-11 gaps, we want 3e-11
        bool feasible = false;
        for (double slack = 1e-9; slack <= 2e-5 && !feasible; slack *= 10) {
            State t = from_packing(start);
            t.S += slack;
            feasible = relax(t.v, t.S, 3000) <= kTol;
            if (feasible) best = t;
        }
        if (!feasible) {
            if (opt.verbose) std::fprintf(stderr, "  could not make the start feasible under margins; returning it unchanged\n");
            return start;
        }
    }
    Packing out = to_packing(run_state(best, opt));
    return container_side(out) < start_side - 1e-10 ? out : start;  // never hand back something no better than the start
}

}  // namespace squeeze
