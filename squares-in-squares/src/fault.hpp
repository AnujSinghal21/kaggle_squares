// Column-DP strip generator ("fault line" composition).
//
// Container side S = k + f.  Grids: UL squares (i+.5, f+j+.5) kept when j-i >= a, LR squares (f+i+.5, j+.5) kept when
// j-i <= b (b < a; r = a-b-1 diagonals are left empty).  The empty corridor is filled with tilted *columns*: a column is a
// 1 x h rectangle (h unit squares stacked) tilted by theta; columns are spaced >= 1 apart along u = (cos, sin) so they are
// independent along w = (-sin, cos) and each may slide to nest in the grid stairs.
//   valid1(u, w)  : a unit square at tilt theta centred at u*U + w*W lies in the container and avoids all kept grid squares
//                   (up to `tol` penetration)
//   column(h) at (u, w) is valid iff the h unit squares at w + (q - (h-1)/2), q = 0..h-1 are all valid
//   best strip    = max over u_0 < u_1 < ... (gaps >= 1) of sum h(u_i)  -> 1-D DP on a u grid
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

#include "strip.hpp"

namespace fault {

struct Column { double u, w; int h; };

struct Result {
    int grid = 0;      // kept grid squares
    int strip = 0;     // squares in columns
    std::vector<Column> cols;
    int total() const { return grid + strip; }
};

struct Params {
    double S, theta;
    int a, b;
    double step = 0.02;      // grid resolution along u and w (1/step must be an integer)
    double tol = 0.015;      // allowed penetration of a column square into the grid
    int max_h = 4;
};

inline Result solve(const Params& P, bool want_cols = false) {
    Result R;
    const double S = P.S;
    const int k = (int)std::floor(S);
    const double f = S - k;
    // kept grid squares (store centres)
    std::vector<double> gx, gy;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++) {
            if (j - i >= P.a) { gx.push_back(i + .5); gy.push_back(f + j + .5); }
            else if (j - i <= P.b) { /* LR uses its own lattice, added below */ }
        }
    const int nul = (int)gx.size();
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            if (j - i <= P.b) { gx.push_back(f + i + .5); gy.push_back(j + .5); }
    (void)nul;
    R.grid = (int)gx.size();
    // spatial hash of kept squares
    const int gxn = (int)(S / 1.5) + 3;
    std::vector<std::vector<int>> cell(gxn * gxn);
    for (int i = 0; i < R.grid; i++) cell[(int)(gy[i] / 1.5) * gxn + (int)(gx[i] / 1.5)].push_back(i);

    const double ct = std::cos(P.theta), st = std::sin(P.theta);
    const double hh = 0.5 * (std::fabs(ct) + std::fabs(st));
    const double r2 = std::sqrt(2.0);
    // corridor in (x - y)/sqrt2 : between the two stairs
    // UL stair squares have centre diff (x-y) = -(f + j - i) <= -(f + a); LR: (f + i - j) >= f - b ... centre line
    const double cmid = 0.5 * ((-(f + P.a)) + (f - P.b)) / r2;
    const double half = (P.a - P.b + 2 * f + 3.0) / r2 * 0.5 + 1.0;
    const int N = (int)std::lround(1.0 / P.step);
    // u range covering the container
    const double umin = -0.2, umax = S * (ct + st) + 0.2;
    const int nu = (int)((umax - umin) / P.step) + 1;
    // w range
    const double wmin = -S * st - 0.2, wmax = S * ct + 0.2;
    const int nw = (int)((wmax - wmin) / P.step) + 1;
    std::vector<unsigned char> ok((size_t)nu * nw, 0);
    for (int iu = 0; iu < nu; iu++) {
        const double u = umin + iu * P.step;
        for (int iw = 0; iw < nw; iw++) {
            const double w = wmin + iw * P.step;
            const double x = u * ct - w * st, y = u * st + w * ct;
            if (x < hh - 1e-9 || y < hh - 1e-9 || x > S - hh + 1e-9 || y > S - hh + 1e-9) continue;
            if (std::fabs((x - y) / r2 - cmid) > half) continue;
            bool good = true;
            const int cx = (int)(x / 1.5), cy = (int)(y / 1.5);
            for (int oy = -1; oy <= 1 && good; oy++)
                for (int ox = -1; ox <= 1 && good; ox++) {
                    const int ax = cx + ox, ay = cy + oy;
                    if (ax < 0 || ay < 0 || ax >= gxn || ay >= gxn) continue;
                    for (int q : cell[ay * gxn + ax])
                        if (strip::pen(gx[q], gy[q], 0.0, x, y, P.theta) > P.tol) { good = false; break; }
                }
            if (good) ok[(size_t)iu * nw + iw] = 1;
        }
    }
    // column heights per u: longest chain with stride N in w
    std::vector<int> H(nu, 0), Hw(nu, 0);
    for (int iu = 0; iu < nu; iu++) {
        int best = 0, bw = 0;
        for (int iw = 0; iw < nw; iw++) {
            int h = 0;
            while (h < P.max_h && iw + h * N < nw && ok[(size_t)iu * nw + iw + h * N]) h++;
            if (h > best) { best = h; bw = iw; }
        }
        H[iu] = best; Hw[iu] = bw;
    }
    // DP: V[iu] = H[iu] + max_{iu' <= iu - N} V[iu']
    std::vector<int> V(nu, 0), pre(nu, -1), pm(nu, 0);  // pm = prefix argmax
    std::vector<int> prefmax(nu, 0), prefarg(nu, -1);
    for (int iu = 0; iu < nu; iu++) {
        int base = 0, from = -1;
        if (iu - N >= 0) { base = prefmax[iu - N]; from = prefarg[iu - N]; }
        V[iu] = H[iu] > 0 ? H[iu] + base : 0;
        pre[iu] = from;
        if (iu == 0) { prefmax[0] = V[0]; prefarg[0] = V[0] > 0 ? 0 : -1; }
        else if (V[iu] > prefmax[iu - 1]) { prefmax[iu] = V[iu]; prefarg[iu] = iu; }
        else { prefmax[iu] = prefmax[iu - 1]; prefarg[iu] = prefarg[iu - 1]; }
    }
    R.strip = prefmax[nu - 1];
    if (want_cols) {
        int cur = prefarg[nu - 1];
        while (cur >= 0 && H[cur] > 0) {
            const int h = H[cur];
            // centre of the column in w: first square at index Hw, last at Hw+(h-1)N
            const double w0 = wmin + Hw[cur] * P.step + (h - 1) * 0.5;
            R.cols.push_back({umin + cur * P.step, w0, h});
            cur = pre[cur];
        }
    }
    return R;
}

}  // namespace fault
