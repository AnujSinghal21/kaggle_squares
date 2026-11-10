// Crop-and-recompose: strip every tilted square out of a packing, shrink the container by delta (grid pieces that rest
// against the right / top wall slide with it), and refill the hole with tilted column slabs found by slab::dp.
#pragma once
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include "common.hpp"
#include "slab.hpp"

namespace recompose {

struct Split {
    std::vector<slab::Sq> axis;   // kept (near axis-aligned) squares, positions in the original frame
    std::vector<slab::Sq> tilt;   // the removed tilted squares (same frame)
    int tilted = 0;
    double L = 0;
};

inline double axis_dev(double a) { a = std::fmod(std::fabs(a), M_PI / 2); return std::min(a, M_PI / 2 - a); }

// Normalise so the container is [0,L]^2, split into kept axis squares and the tilted count.
inline Split split(const Packing& p, double tilt_thresh = 0.03) {
    Split s;
    s.L = container_side(p);
    double lo_x = 1e18, lo_y = 1e18;
    for (const Square& q : p) {
        const double h = 0.5 * (std::fabs(std::cos(q.angle)) + std::fabs(std::sin(q.angle)));
        lo_x = std::min(lo_x, q.x - h); lo_y = std::min(lo_y, q.y - h);
    }
    for (const Square& q : p) {
        if (axis_dev(q.angle) < tilt_thresh) s.axis.push_back({q.x - lo_x, q.y - lo_y, q.angle});
        else { s.tilted++; s.tilt.push_back({q.x - lo_x, q.y - lo_y, q.angle}); }
    }
    return s;
}

// Shift for each kept square when the container shrinks by delta: squares connected (through side contacts) to the
// right wall move left, to the top wall move down; others stay.  A shift that would leave the container is dropped.
inline std::vector<slab::Sq> shrink(const Split& s, double delta, double contact = 1e-3) {
    const int n = (int)s.axis.size();
    const double L = s.L;
    std::vector<char> right(n, 0), top(n, 0);
    std::vector<int> st;
    // build adjacency lazily by scanning pairs near each other (n up to ~1000: O(n^2) is fine here)
    std::vector<std::vector<int>> hadj(n), vadj(n);   // horizontal / vertical neighbours
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            const double dx = s.axis[j].x - s.axis[i].x, dy = s.axis[j].y - s.axis[i].y;
            if (std::fabs(std::fabs(dx) - 1) < contact && std::fabs(dy) < 1 - 1e-3 + 0.0) { hadj[i].push_back(j); hadj[j].push_back(i); }
            if (std::fabs(std::fabs(dy) - 1) < contact && std::fabs(dx) < 1 - 1e-3 + 0.0) { vadj[i].push_back(j); vadj[j].push_back(i); }
        }
    for (int i = 0; i < n; i++) {
        if (s.axis[i].x + 0.5 > L - contact) { right[i] = 1; st.push_back(i); }
    }
    while (!st.empty()) { const int i = st.back(); st.pop_back(); for (int j : hadj[i]) if (!right[j]) { right[j] = 1; st.push_back(j); } }
    for (int i = 0; i < n; i++) if (s.axis[i].y + 0.5 > L - contact) { top[i] = 1; st.push_back(i); }
    while (!st.empty()) { const int i = st.back(); st.pop_back(); for (int j : vadj[i]) if (!top[j]) { top[j] = 1; st.push_back(j); } }
    std::vector<slab::Sq> out;
    for (int i = 0; i < n; i++) {
        slab::Sq q = s.axis[i];
        if (right[i] && q.x - 0.5 - delta > -1e-9) q.x -= delta;
        if (top[i] && q.y - 0.5 - delta > -1e-9) q.y -= delta;
        out.push_back(q);
    }
    return out;
}

// Connected clusters of tilted squares (centre distance < 1.6) with their centre bounding boxes.
struct Cluster { std::vector<int> idx; double x0, x1, y0, y1; };
inline std::vector<Cluster> clusters(const Split& s, double link = 1.6) {
    const int m = (int)s.tilt.size();
    std::vector<int> par(m);
    for (int i = 0; i < m; i++) par[i] = i;
    auto find = [&](int a) { while (par[a] != a) a = par[a] = par[par[a]]; return a; };
    for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++) {
            const double dx = s.tilt[i].x - s.tilt[j].x, dy = s.tilt[i].y - s.tilt[j].y;
            if (dx * dx + dy * dy < link * link) par[find(i)] = find(j);
        }
    std::map<int, Cluster> g;
    for (int i = 0; i < m; i++) {
        Cluster& c = g[find(i)];
        if (c.idx.empty()) { c.x0 = c.x1 = s.tilt[i].x; c.y0 = c.y1 = s.tilt[i].y; }
        c.idx.push_back(i);
        c.x0 = std::min(c.x0, s.tilt[i].x); c.x1 = std::max(c.x1, s.tilt[i].x);
        c.y0 = std::min(c.y0, s.tilt[i].y); c.y1 = std::max(c.y1, s.tilt[i].y);
    }
    std::vector<Cluster> out;
    for (auto& kv : g) out.push_back(kv.second);
    std::sort(out.begin(), out.end(), [](const Cluster& a, const Cluster& b) { return a.idx.size() > b.idx.size(); });
    return out;
}

}  // namespace recompose
