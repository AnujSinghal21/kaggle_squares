// Tilted-block composer: one W x H block of unit squares tilted by theta (centre (cx,cy)) inside four corner-anchored
// axis lattices.  Lattice squares overlapping the block deeper than `tol` are dropped; the four lattices (UL, UR, LL, LR)
// are assigned by quadrant around the block centre and resolved greedily where they meet.  This is the structure of the
// 45-degree "diamond" packings (e.g. n=298: an 11x11 block, f = frac(11*sqrt2)/2) and of tilted-island packings.
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

#include "strip.hpp"

namespace block {

struct Params {
    double S, theta, cx, cy;
    int W, H;
    double tol = 0.02;
};

struct Built {
    std::vector<double> v;   // block squares first
    int n_block = 0;
    int total() const { return (int)v.size() / 3; }
};

inline Built build(const Params& P) {
    Built B;
    const double S = P.S, ct = std::cos(P.theta), st = std::sin(P.theta);
    const double hh = 0.5 * (std::fabs(ct) + std::fabs(st));
    for (int i = 0; i < P.W; i++)
        for (int j = 0; j < P.H; j++) {
            const double a = i + 0.5 - P.W * 0.5, b = j + 0.5 - P.H * 0.5;
            const double x = P.cx + a * ct - b * st, y = P.cy + a * st + b * ct;
            if (x < hh - P.tol || y < hh - P.tol || x > S - hh + P.tol || y > S - hh + P.tol) { B.v.clear(); B.n_block = -1; return B; }
            B.v.push_back(x); B.v.push_back(y); B.v.push_back(P.theta);
        }
    B.n_block = P.W * P.H;
    const int K = (int)std::floor(S);
    // bucket the block squares for quick rejection
    const int nb = B.n_block;
    std::vector<double> kx, ky;   // kept axis squares
    auto conflicts = [&](double x, double y) {
        for (int q = 0; q < nb; q++) {
            const double dx = B.v[3 * q] - x, dy = B.v[3 * q + 1] - y;
            if (dx * dx + dy * dy < 2.0 && strip::pen(x, y, 0.0, B.v[3 * q], B.v[3 * q + 1], B.v[3 * q + 2]) > P.tol) return true;
        }
        for (size_t q = 0; q < kx.size(); q++)
            if (std::fabs(kx[q] - x) < 1 - P.tol && std::fabs(ky[q] - y) < 1 - P.tol) return true;
        return false;
    };
    // quadrant lattices: (right?, top?)  x = right ? S-i-.5 : i+.5 ;  y = top ? S-j-.5 : j+.5
    const int order[4][2] = {{0, 1}, {1, 1}, {0, 0}, {1, 0}};  // UL, UR, LL, LR
    for (const auto& o : order) {
        const bool right = o[0], top = o[1];
        for (int i = 0; i < K; i++)
            for (int j = 0; j < K; j++) {
                const double x = right ? S - i - .5 : i + .5, y = top ? S - j - .5 : j + .5;
                if ((x > P.cx) != right || (y > P.cy) != top) continue;
                if (!conflicts(x, y)) { kx.push_back(x); ky.push_back(y); }
            }
    }
    for (size_t q = 0; q < kx.size(); q++) { B.v.push_back(kx[q]); B.v.push_back(ky[q]); B.v.push_back(0); }
    return B;
}

}  // namespace block
