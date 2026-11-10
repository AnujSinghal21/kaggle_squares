// Window re-solve: a few "free" squares move inside a fixed environment of "frozen" squares and the container walls.
// Same penetration-depth energy as squeeze.hpp, but only free squares are variables, and free-vs-frozen pairs come from a
// static spatial hash, so one relax costs ~10x less than the full-packing engine.  Used to re-pack a tilted cluster
// (strip, slab, diamond edge) at a smaller container side while the grid pieces around it keep their lattice positions.
#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#include "common.hpp"
#include "squeeze.hpp"

namespace window {

constexpr double kWinTol = 1e-16;   // 'feasible enough'; squeeze::polish later makes it exact with real gaps

struct Fixed { double x, y, t; };

struct Env {
    double S = 0;
    bool walls = true;          // false for periodic cells
    double yoff = 0;            // box is [0,S] x [0,S+yoff]
    bool wall_y = true;         // false: only the left/right walls (vertical periodic strips)
    double band_lo = -1e300, band_hi = 1e300;   // free squares must keep y - x of their centre in [band_lo, band_hi] (periodic cells)
    double wrap_a0 = -1e300;    // periodic cells: keep (x+y)/2 of every free square in [wrap_a0, wrap_a0 + px)
    double px = 0, py = 0;      // period vector (0,0 = not periodic): free squares also see their images shifted by +-(px,py)
    std::vector<Fixed> frozen;
    // static hash of frozen squares, cell 1.5
    int gx = 0;
    std::vector<std::vector<int>> cell;
    void index() {
        gx = (int)std::ceil((S + std::max(yoff, 0.0) + 2) / 1.5) + 2;
        cell.assign((size_t)gx * gx, {});
        for (int i = 0; i < (int)frozen.size(); i++) {
            const int cx = std::min(gx - 1, std::max(0, (int)((frozen[i].x + 1) / 1.5))), cy = std::min(gx - 1, std::max(0, (int)((frozen[i].y + 1) / 1.5)));
            cell[(size_t)cy * gx + cx].push_back(i);
        }
    }
};

// Energy of free squares v (x,y,t triples) in env; g gets the gradient wrt v.
inline double energy(const Env& env, const std::vector<double>& v, std::vector<double>& g) {
    using squeeze::kMargin;
    const int n = (int)v.size() / 3;
    g.assign(v.size(), 0.0);
    double E = 0;
    const double S = env.S;
    for (int i = 0; i < n && env.walls; i++) {
        const bool wy = env.wall_y;
        const double x = v[3 * i], y = v[3 * i + 1], t = v[3 * i + 2];
        const double c = std::cos(t), s = std::sin(t);
        const double h = 0.5 * (std::fabs(c) + std::fabs(s)) + kMargin;
        const double dh = 0.5 * (-s * (c >= 0 ? 1 : -1) + c * (s >= 0 ? 1 : -1));
        const double w[4] = {h - x, x + h - S, h - y, y + h - S - env.yoff};
        const double dx[4] = {-1, 1, 0, 0}, dy[4] = {0, 0, -1, 1};
        for (int k = 0; k < 4; k++)
            if (w[k] > 0 && (wy || k < 2)) { E += w[k] * w[k]; g[3 * i] += 2 * w[k] * dx[k]; g[3 * i + 1] += 2 * w[k] * dy[k]; g[3 * i + 2] += 2 * w[k] * dh; }
    }
    if (env.band_hi < 1e200)
        for (int i = 0; i < n; i++) {
            const double dd = v[3 * i + 1] - v[3 * i];
            const double viol = dd > env.band_hi ? dd - env.band_hi : (dd < env.band_lo ? dd - env.band_lo : 0.0);
            if (viol != 0) { E += viol * viol; g[3 * i + 1] += 2 * viol; g[3 * i] -= 2 * viol; }
        }
    // pair kernel: penetration between square (xi,yi,ti) and (xj,yj,tj); adds to gi (+opt gj)
    auto pair = [&](double xi, double yi, double ti, double xj, double yj, double tj, double* gi, double* gj) {
        const double dx = xj - xi, dy = yj - yi;
        if (dx * dx + dy * dy > 2.0 + 4 * kMargin + 1e-6) return;
        const double ci = std::cos(ti), si = std::sin(ti), cj = std::cos(tj), sj = std::sin(tj);
        const double phi = tj - ti, cp = std::cos(phi), sp = std::sin(phi);
        const double H = 0.5 + 0.5 * (std::fabs(cp) + std::fabs(sp)) + kMargin;
        const double dH = 0.5 * (-sp * (cp >= 0 ? 1 : -1) + cp * (sp >= 0 ? 1 : -1));
        const double d[4] = {dx * ci + dy * si, -dx * si + dy * ci, dx * cj + dy * sj, -dx * sj + dy * cj};
        int k = 0;
        for (int q = 1; q < 4; q++) if (std::fabs(d[q]) > std::fabs(d[k])) k = q;
        const double P = H - std::fabs(d[k]);
        if (P <= 0) return;
        const double sg = d[k] >= 0 ? 1.0 : -1.0;
        double ddx, ddy, ddti = 0, ddtj = 0;
        if (k == 0)      { ddx = ci;  ddy = si;  ddti = d[1]; }
        else if (k == 1) { ddx = -si; ddy = ci;  ddti = -d[0]; }
        else if (k == 2) { ddx = cj;  ddy = sj;  ddtj = d[3]; }
        else             { ddx = -sj; ddy = cj;  ddtj = -d[2]; }
        E += P * P;
        const double f = 2 * P, gdx = f * (-sg * ddx), gdy = f * (-sg * ddy);
        gi[0] -= gdx; gi[1] -= gdy; gi[2] += f * (-dH - sg * ddti);
        if (gj) { gj[0] += gdx; gj[1] += gdy; gj[2] += f * (dH - sg * ddtj); }
    };
    // free-free: O(n^2) with distance reject (n is ~100-300)
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            const double dx = v[3 * j] - v[3 * i], dy = v[3 * j + 1] - v[3 * i + 1];
            if (dx * dx + dy * dy > 2.0 + 1e-5) continue;
            pair(v[3 * i], v[3 * i + 1], v[3 * i + 2], v[3 * j], v[3 * j + 1], v[3 * j + 2], &g[3 * i], &g[3 * j]);
        }
    if (env.px != 0 || env.py != 0)   // periodic images (shift +-1 period) of free-free pairs, including a square with its own image
        for (int sgn = 1; sgn <= 1; sgn++)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    const double xj = v[3 * j] + env.px, yj = v[3 * j + 1] + env.py;
                    const double dx = xj - v[3 * i], dy = yj - v[3 * i + 1];
                    if (dx * dx + dy * dy > 2.0 + 1e-5) continue;
                    pair(v[3 * i], v[3 * i + 1], v[3 * i + 2], xj, yj, v[3 * j + 2], &g[3 * i], &g[3 * j]);
                }
    // free-frozen
    for (int i = 0; i < n; i++) {
        const int cx = std::min(env.gx - 1, std::max(0, (int)((v[3 * i] + 1) / 1.5))), cy = std::min(env.gx - 1, std::max(0, (int)((v[3 * i + 1] + 1) / 1.5)));
        for (int oy = -1; oy <= 1; oy++)
            for (int ox = -1; ox <= 1; ox++) {
                const int ax = cx + ox, ay = cy + oy;
                if (ax < 0 || ay < 0 || ax >= env.gx || ay >= env.gx) continue;
                for (int q : env.cell[(size_t)ay * env.gx + ax]) {
                    const Fixed& f = env.frozen[q];
                    pair(v[3 * i], v[3 * i + 1], v[3 * i + 2], f.x, f.y, f.t, &g[3 * i], nullptr);
                }
            }
    }
    return E;
}

inline void wrap(const Env& env, std::vector<double>& v) {
    if (env.wrap_a0 < -1e200) return;
    if (env.px == 0 && env.py != 0) {   // vertical period: wrap y
        for (size_t i = 0; i + 2 < v.size(); i += 3) v[i + 1] -= std::floor((v[i + 1] - env.wrap_a0) / env.py) * env.py;
        return;
    }
    if (env.px == 0) return;
    for (size_t i = 0; i + 2 < v.size(); i += 3) {
        const double a = 0.5 * (v[i] + v[i + 1]);
        const double sh = std::floor((a - env.wrap_a0) / env.px) * env.px;
        if (sh != 0) { v[i] -= sh; v[i + 1] -= sh; }
    }
}

inline double relax(const Env& env, std::vector<double>& v, int max_iter) {
    wrap(env, v);
    const double E = squeeze::lbfgs(v, [&](const std::vector<double>& x, std::vector<double>& g) { return energy(env, x, g); }, max_iter,
                          [](double E) { return E <= kWinTol; });
    wrap(env, v);
    return E;
}

struct HopOptions {
    double seconds = 30;
    int threads = 8;
    uint64_t seed = 1;
    int max_iter = 300;
    bool verbose = false;
    double shift_dx = 0, shift_dy = 0;   // if non-zero: extra move 'translate every free square along (dx,dy) by a random amount'
    double shift_sigma = 0.25;
    double temp = 0.03;          // Metropolis temperature relative to the current energy
    double p_local = 1.0 / 3, p_all = 1.0 / 3;   // move mix: local kick / whole-set kick / relocate (rest)
    double r_lo = 0.8, r_span = 2.5;             // kick radius range
};

// Basin hopping on the free squares at fixed container side.  Returns true if some thread reached E <= kTol; v_out holds it.
inline bool hop(const Env& env, const std::vector<double>& start, const HopOptions& opt, std::vector<double>& v_out, double* best_E_out = nullptr) {
    const int n = (int)start.size() / 3;
    std::atomic<bool> done{false};
    std::mutex mu;
    double best_E = 1e300;
    const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(opt.seconds);
    auto worker = [&](int tid) {
        std::mt19937_64 rng(opt.seed * 1000003ULL + tid * 7919ULL);
        std::uniform_real_distribution<double> U(0, 1);
        std::normal_distribution<double> N01(0, 1);
        std::vector<double> cur = start;
        double cur_E = relax(env, cur, opt.max_iter);
        long it = 0;
        while (!done && std::chrono::steady_clock::now() < t_end) {
            if (cur_E <= kWinTol) { std::lock_guard<std::mutex> lk(mu); if (!done) { done = true; v_out = cur; best_E = cur_E; } return; }
            std::vector<double> trial = cur;
            int mode; { const double u = U(rng); mode = u < opt.p_local ? 0 : (u < opt.p_local + opt.p_all ? 1 : 2); }
            if ((opt.shift_dx != 0 || opt.shift_dy != 0) && U(rng) < 0.25) mode = 3;
            if (mode == 3) {   // global slide along the strip direction
                const double a = opt.shift_sigma * N01(rng);
                for (int i = 0; i < n; i++) { trial[3 * i] += a * opt.shift_dx; trial[3 * i + 1] += a * opt.shift_dy; }
            } else if (mode == 2) {  // relocate: pick a random free square, put it at a random free-ish spot (near another random square)
                const int a = (int)(U(rng) * n), b = (int)(U(rng) * n);
                trial[3 * a] = cur[3 * b] + (U(rng) - 0.5) * 1.5; trial[3 * a + 1] = cur[3 * b + 1] + (U(rng) - 0.5) * 1.5; trial[3 * a + 2] = U(rng) * M_PI / 2;
            } else {
                const int c = (int)(U(rng) * n);
                const double r = opt.r_lo + opt.r_span * U(rng), sp = (U(rng) < 0.5 ? 0.03 : 0.2), sa = (U(rng) < 0.5 ? 0.05 : 0.5);
                for (int i = 0; i < n; i++) {
                    const double dx = cur[3 * i] - cur[3 * c], dy = cur[3 * i + 1] - cur[3 * c + 1];
                    if (mode == 0 && dx * dx + dy * dy > r * r) continue;
                    trial[3 * i] += sp * N01(rng); trial[3 * i + 1] += sp * N01(rng); trial[3 * i + 2] += sa * N01(rng);
                }
            }
            const double E = relax(env, trial, opt.max_iter);
            const double T = opt.temp * std::max(cur_E, 1e-6);
            if (E < cur_E || U(rng) < std::exp(-(E - cur_E) / T)) { cur = trial; cur_E = E; }
            if (++it % 500 == 0) { std::lock_guard<std::mutex> lk(mu); if (cur_E < best_E) best_E = cur_E; if (opt.verbose && tid == 0) std::fprintf(stderr, "  [t0] it %ld E=%.3e best %.3e\n", it, cur_E, best_E); }
        }
        std::lock_guard<std::mutex> lk(mu);
        if (cur_E < best_E) best_E = cur_E;
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < opt.threads; t++) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    if (best_E_out) *best_E_out = best_E;
    return done;
}

}  // namespace window
