// Column DP in an arbitrary free region ("slab fill").
//
// Given obstacles (axis-aligned or tilted unit squares), a container [0,S]^2 and a tilt theta, fill the free area with
// tilted columns: a column is a stack of h unit squares along w = (-sin, cos); columns are spaced >= 1 apart along
// u = (cos, sin), so they are independent along w and each slides to nest against the obstacles/walls.  The best
// total count is a 1-D DP over u (same idea as fault.hpp, but the "stairs" are any obstacle list and the region can be
// limited by a window).  Used for the tilted slabs between an integer grid core and the container walls.
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

#include "strip.hpp"

namespace slab {

struct Sq { double x, y, t; };
struct Column { double u, w; int h; };

struct Window { double x0 = 0, x1 = 1e9, y0 = 0, y1 = 1e9; };   // square centres' hull must lie inside (up to tol)

struct Fill {
    int count = 0;
    std::vector<Column> cols;
    std::vector<Sq> squares;   // filled only when want_squares
};

inline Fill dp(double S, double theta, const std::vector<Sq>& obs, const Window& win, double step, double tol,
               int max_h, bool want_squares = false) {
    Fill F;
    const double ct = std::cos(theta), st = std::sin(theta);
    const double hh = 0.5 * (std::fabs(ct) + std::fabs(st));
    // spatial hash of obstacles
    const int gxn = (int)(S / 1.5) + 3;
    std::vector<std::vector<int>> cell((size_t)gxn * gxn);
    for (int i = 0; i < (int)obs.size(); i++) {
        const int cx = std::min(gxn - 1, std::max(0, (int)(obs[i].x / 1.5))), cy = std::min(gxn - 1, std::max(0, (int)(obs[i].y / 1.5)));
        cell[(size_t)cy * gxn + cx].push_back(i);
    }
    const double wx0 = std::max(win.x0, 0.0) + hh - tol, wx1 = std::min(win.x1, S) - hh + tol;
    const double wy0 = std::max(win.y0, 0.0) + hh - tol, wy1 = std::min(win.y1, S) - hh + tol;
    if (wx1 < wx0 || wy1 < wy0) return F;
    // bounding box of the window (centre positions) in (u, w)
    double umin = 1e18, umax = -1e18, wmin = 1e18, wmax = -1e18;
    for (double x : {wx0, wx1})
        for (double y : {wy0, wy1}) {
            const double u = x * ct + y * st, w = -x * st + y * ct;
            umin = std::min(umin, u); umax = std::max(umax, u); wmin = std::min(wmin, w); wmax = std::max(wmax, w);
        }
    const int N = (int)std::lround(1.0 / step);
    const int nu = (int)((umax - umin) / step) + 1, nw = (int)((wmax - wmin) / step) + 1;
    std::vector<unsigned char> ok((size_t)nu * nw, 0);
    for (int iu = 0; iu < nu; iu++) {
        const double u = umin + iu * step;
        for (int iw = 0; iw < nw; iw++) {
            const double w = wmin + iw * step;
            const double x = u * ct - w * st, y = u * st + w * ct;
            if (x < wx0 || x > wx1 || y < wy0 || y > wy1) continue;
            bool good = true;
            const int cx = (int)(x / 1.5), cy = (int)(y / 1.5);
            for (int oy = -1; oy <= 1 && good; oy++)
                for (int ox = -1; ox <= 1 && good; ox++) {
                    const int ax = cx + ox, ay = cy + oy;
                    if (ax < 0 || ay < 0 || ax >= gxn || ay >= gxn) continue;
                    for (int q : cell[(size_t)ay * gxn + ax])
                        if (strip::pen(obs[q].x, obs[q].y, obs[q].t, x, y, theta) > tol) { good = false; break; }
                }
            if (good) ok[(size_t)iu * nw + iw] = 1;
        }
    }
    std::vector<int> H(nu, 0), Hw(nu, 0);
    for (int iu = 0; iu < nu; iu++) {
        int best = 0, bw = 0;
        const unsigned char* row = &ok[(size_t)iu * nw];
        for (int iw = 0; iw < nw; iw++) {
            if (!row[iw]) continue;
            int h = 1;
            while (h < max_h && iw + h * N < nw && row[iw + h * N]) h++;
            if (h > best) { best = h; bw = iw; }
        }
        H[iu] = best; Hw[iu] = bw;
    }
    std::vector<int> V(nu, 0), pre(nu, -1), pm(nu, 0), pa(nu, -1);
    for (int iu = 0; iu < nu; iu++) {
        int base = 0, from = -1;
        if (iu - N >= 0) { base = pm[iu - N]; from = pa[iu - N]; }
        V[iu] = H[iu] > 0 ? H[iu] + base : 0;
        pre[iu] = from;
        if (iu == 0) { pm[0] = V[0]; pa[0] = V[0] > 0 ? 0 : -1; }
        else if (V[iu] > pm[iu - 1]) { pm[iu] = V[iu]; pa[iu] = iu; }
        else { pm[iu] = pm[iu - 1]; pa[iu] = pa[iu - 1]; }
    }
    F.count = nu ? pm[nu - 1] : 0;
    if (want_squares) {
        int cur = nu ? pa[nu - 1] : -1;
        while (cur >= 0 && H[cur] > 0) {
            const int h = H[cur];
            const double u = umin + cur * step, w0 = wmin + Hw[cur] * step;
            F.cols.push_back({u, w0, h});
            for (int q = 0; q < h; q++) {
                const double w = w0 + q;
                F.squares.push_back({u * ct - w * st, u * st + w * ct, theta});
            }
            cur = pre[cur];
        }
    }
    return F;
}

// Greedily add axis-aligned squares from the four corner-anchored lattices that fit in the leftover free space.
inline int fill_corners(double S, std::vector<Sq>& placed, double tol, std::vector<Sq>* added = nullptr) {
    const int K = (int)std::floor(S);
    int cnt = 0;
    const int gxn = (int)(S / 1.5) + 3;
    std::vector<std::vector<int>> cell((size_t)gxn * gxn);
    auto cid = [&](double x, double y) { return std::min(gxn - 1, std::max(0, (int)(y / 1.5))) * gxn + std::min(gxn - 1, std::max(0, (int)(x / 1.5))); };
    for (int i = 0; i < (int)placed.size(); i++) cell[cid(placed[i].x, placed[i].y)].push_back(i);
    const int order[4][2] = {{0, 1}, {1, 0}, {0, 0}, {1, 1}};
    for (const auto& o : order)
        for (int i = 0; i < K; i++)
            for (int j = 0; j < K; j++) {
                const double x = o[0] ? S - i - .5 : i + .5, y = o[1] ? S - j - .5 : j + .5;
                bool good = true;
                const int cx = std::min(gxn - 1, (int)(x / 1.5)), cy = std::min(gxn - 1, (int)(y / 1.5));
                for (int oy = -1; oy <= 1 && good; oy++)
                    for (int ox = -1; ox <= 1 && good; ox++) {
                        const int ax = cx + ox, ay = cy + oy;
                        if (ax < 0 || ay < 0 || ax >= gxn || ay >= gxn) continue;
                        for (int q : cell[(size_t)ay * gxn + ax])
                            if (strip::pen(placed[q].x, placed[q].y, placed[q].t, x, y, 0.0) > tol) { good = false; break; }
                    }
                if (!good) continue;
                placed.push_back({x, y, 0.0});
                cell[cid(x, y)].push_back((int)placed.size() - 1);
                if (added) added->push_back({x, y, 0.0});
                cnt++;
            }
    return cnt;
}

}  // namespace slab
