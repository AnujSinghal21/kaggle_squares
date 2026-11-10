// Compose a full packing from a periodic strip motif (found by tools/periodic_hop).
//   Cell convention (periodic_hop): UL squares (i+.5, f+j+.5) kept when j-i >= 0, LR squares (f+i+.5, j+.5) kept when
//   j-i <= -1-r; motif squares repeat with period (p, p).  A container of side S = k + f uses the same lattices;
//   shifting the motif by (0, a) moves the corridor to the cut j-i >= a / j-i <= a-1-r.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "strip.hpp"

namespace motif {

struct Motif { int r = 2, p = 4; double f = 0; std::vector<double> sq; };   // x,y,theta(rad) triples

inline bool load(const std::string& path, Motif& m, int r, int p) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    m.r = r; m.p = p; m.sq.clear();
    while (std::getline(in, line)) {
        if (line.rfind("r=", 0) == 0) { const size_t k = line.find("f="); if (k != std::string::npos) m.f = std::atof(line.c_str() + k + 2); continue; }
        std::istringstream ss(line);
        double x, y, d;
        if (ss >> x >> y >> d) { m.sq.push_back(x); m.sq.push_back(y); m.sq.push_back(d * M_PI / 180); }
    }
    return !m.sq.empty();
}

struct Built { std::vector<double> v; int n_strip = 0, grid = 0; int total() const { return (int)v.size() / 3; } };

// Container of side S; corridor offset a; diagonal phase s (integer); strict = drop strip squares that touch walls.
inline Built build(const Motif& m, double S, int a, int s, double tol = 2e-4) {
    Built B;
    const int k = (int)std::floor(S);
    const double f = S - k;
    // strip: tile
    const int T = k / m.p + 3;
    for (int t = -T; t <= T; t++)
        for (size_t q = 0; q < m.sq.size(); q += 3) {
            const double x = m.sq[q] + t * m.p + s, y = m.sq[q + 1] + t * m.p + s + a, th = m.sq[q + 2];
            const double hh = 0.5 * (std::fabs(std::cos(th)) + std::fabs(std::sin(th)));
            if (x < hh - tol || y < hh - tol || x > S - hh + tol || y > S - hh + tol) continue;
            B.v.push_back(x); B.v.push_back(y); B.v.push_back(th);
        }
    B.n_strip = (int)B.v.size() / 3;
    const int ns = B.n_strip;
    auto hits = [&](double x, double y) {
        for (int i = 0; i < ns; i++) {
            const double dx = B.v[3 * i] - x, dy = B.v[3 * i + 1] - y;
            if (dx * dx + dy * dy < 2.0 && strip::pen(x, y, 0.0, B.v[3 * i], B.v[3 * i + 1], B.v[3 * i + 2]) > tol) return true;
        }
        return false;
    };
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++) {
            if (j - i >= a) { const double x = i + .5, y = f + j + .5; if (!hits(x, y)) { B.v.push_back(x); B.v.push_back(y); B.v.push_back(0); B.grid++; } }
            else if (j - i <= a - 1 - m.r) { const double x = f + i + .5, y = j + .5; if (!hits(x, y)) { B.v.push_back(x); B.v.push_back(y); B.v.push_back(0); B.grid++; } }
        }
    return B;
}

}  // namespace motif
