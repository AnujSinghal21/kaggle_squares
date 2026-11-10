// Diagonal-strip generator: two axis-aligned grids anchored in opposite corners, separated by a band of
// tilted squares (a rotated lattice) along the main diagonal.  This is the structure of every best-known
// packing with n > ~300 that sits at L = k + 0.8 (see analysis: strips of 2-wide squares tilted ~24 deg).
//
//   UL grid: centres (i+.5, S-j-.5)    LR grid: centres (S-i-.5, j+.5)   (offset by (f,-f) from each other, f = S-floor(S))
//   strip  : lattice  o + a*u + b*w  (u = (cos th, sin th), w = (-sin th, cos th)), kept if its centre lies in the
//            band |(x-y)/sqrt2 - s| <= h and the square fits in [0,S]^2.
// Grid squares that penetrate a strip square (or an already kept grid square) by more than `tol` are dropped.
// The result is only a *start*: squeeze::relax closes the small residual overlaps.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "common.hpp"

namespace strip {

struct Config {
    double S = 0;            // container side
    double theta = 0.42;     // lattice tilt (rad)
    double px = 0, py = 0;   // lattice phase: reference point offset from the container centre (in lattice axes)
    double s = 0;            // band centre, signed distance from the main diagonal (x - y)/sqrt2
    double h = 1.3;          // band half width (centres)
    double tol = 0.04;       // max penetration kept at generation time
};

// Penetration depth (>0 = overlapping) between unit squares, closed form (same as squeeze.hpp).
inline double pen(double xi, double yi, double ti, double xj, double yj, double tj) {
    const double dx = xj - xi, dy = yj - yi;
    if (dx * dx + dy * dy >= 2.0) return -1;
    const double ci = std::cos(ti), si = std::sin(ti), cj = std::cos(tj), sj = std::sin(tj);
    const double phi = tj - ti;
    const double H = 0.5 + 0.5 * (std::fabs(std::cos(phi)) + std::fabs(std::sin(phi)));
    const double d[4] = {dx * ci + dy * si, -dx * si + dy * ci, dx * cj + dy * sj, -dx * sj + dy * cj};
    double m = 0;
    for (int q = 0; q < 4; q++) m = std::max(m, std::fabs(d[q]));
    return H - m;
}

struct Generated {
    std::vector<double> v;  // x,y,t triples
    int n_strip = 0;        // the first n_strip squares are the tilted ones
};

// Candidate set (strip + both grids). Squares are flagged in `kind`: 0 strip, 1 UL, 2 LR.
inline Generated generate(const Config& c) {
    Generated g;
    const double S = c.S, r2 = std::sqrt(2.0);
    const double ct = std::cos(c.theta), st = std::sin(c.theta);
    const double hh = 0.5 * (std::fabs(ct) + std::fabs(st));
    // strip lattice
    const double ox = S / 2 + c.px * ct - c.py * st, oy = S / 2 + c.px * st + c.py * ct;
    const int R = (int)std::ceil(S * 1.2) + 2;
    for (int a = -R; a <= R; a++)
        for (int b = -R; b <= R; b++) {
            const double x = ox + a * ct - b * st, y = oy + a * st + b * ct;
            if (x < hh || y < hh || x > S - hh || y > S - hh) continue;
            if (std::fabs((x - y) / r2 - c.s) > c.h) continue;
            g.v.push_back(x); g.v.push_back(y); g.v.push_back(c.theta);
        }
    g.n_strip = (int)g.v.size() / 3;
    const int ns = g.n_strip;
    auto hits_strip = [&](double x, double y) {
        for (int i = 0; i < ns; i++) {
            const double dx = g.v[3 * i] - x, dy = g.v[3 * i + 1] - y;
            if (dx * dx + dy * dy > 2.0) continue;
            if (pen(x, y, 0, g.v[3 * i], g.v[3 * i + 1], g.v[3 * i + 2]) > c.tol) return true;
        }
        return false;
    };
    const int K = (int)std::floor(S);
    std::vector<double> ul, lr;
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++) {
            const double x = i + 0.5, y = S - j - 0.5;
            if ((x - y) / r2 >= c.s || hits_strip(x, y)) continue;
            ul.push_back(x); ul.push_back(y);
        }
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++) {
            const double x = S - i - 0.5, y = j + 0.5;
            if ((x - y) / r2 < c.s || hits_strip(x, y)) continue;
            // against kept UL squares (axis-aligned both): overlap if both |dx|,|dy| < 1 - small
            bool bad = false;
            for (size_t q = 0; q < ul.size() && !bad; q += 2)
                if (std::fabs(ul[q] - x) < 1 - c.tol && std::fabs(ul[q + 1] - y) < 1 - c.tol) bad = true;
            if (bad) continue;
            lr.push_back(x); lr.push_back(y);
        }
    for (size_t q = 0; q < ul.size(); q += 2) { g.v.push_back(ul[q]); g.v.push_back(ul[q + 1]); g.v.push_back(0); }
    for (size_t q = 0; q < lr.size(); q += 2) { g.v.push_back(lr[q]); g.v.push_back(lr[q + 1]); g.v.push_back(0); }
    return g;
}

// Remove squares until exactly n remain: worst-penetrating first, then random grid squares touching the strip.
inline bool trim_to(std::vector<double>& v, int n_strip, int n, std::mt19937_64& rng) {
    int cnt = (int)v.size() / 3;
    if (cnt < n) return false;
    while (cnt > n) {
        std::vector<double> sum(cnt, 0.0);
        for (int i = 0; i < cnt; i++)
            for (int j = i + 1; j < cnt; j++) {
                const double p = pen(v[3 * i], v[3 * i + 1], v[3 * i + 2], v[3 * j], v[3 * j + 1], v[3 * j + 2]);
                if (p > 0) { sum[i] += p; sum[j] += p; }
            }
        int pick = -1;
        double bestv = 0;
        for (int i = 0; i < cnt; i++) if (sum[i] > bestv) { bestv = sum[i]; pick = i; }
        if (pick < 0) {   // no overlaps: drop a random grid square near the strip (most slack where it matters)
            std::vector<int> near;
            for (int i = n_strip; i < cnt; i++)
                for (int j = 0; j < n_strip; j++) {
                    const double dx = v[3 * i] - v[3 * j], dy = v[3 * i + 1] - v[3 * j + 1];
                    if (dx * dx + dy * dy < 6.0) { near.push_back(i); break; }
                }
            if (near.empty()) pick = n_strip + (int)(rng() % (uint64_t)std::max(1, cnt - n_strip));
            else pick = near[rng() % near.size()];
        }
        v.erase(v.begin() + 3 * pick, v.begin() + 3 * pick + 3);
        if (pick < n_strip) n_strip--;
        cnt--;
    }
    return true;
}

}  // namespace strip
