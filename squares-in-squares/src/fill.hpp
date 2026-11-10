// Exact-ish upright refill around a fixed tilted core.
//   Given a box [0,Sx] x [0,Sy] and tilted squares T, pick the largest set of axis-aligned unit squares that fits in the box without
//   overlapping T or each other.  Candidates: the four wall-anchored integer lattices, plus squares flush against the extreme
//   vertices of every tilted square (up to `depth` steps), each paired with the lattice coordinates of the other axis.
//   The conflict graph among uprights is a "unit-square intersection" graph, so a simple iterated local search finds the optimum
//   in practice; a 1-2 swap local search plus random restarts is used (cf. the ARW heuristic).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <set>
#include <vector>

#include "common.hpp"
#include "strip.hpp"

namespace fill {

struct P2 { double x, y; };

struct Result {
    std::vector<P2> up;     // chosen upright centres
    int candidates = 0;
};

inline double snap(double v) { return std::round(v * 1e9) / 1e9; }

// bit (ax + 2*ay) set if the centre lies on the wall-anchored lattice (ax: 0 = from the left wall, 1 = from the right; ay likewise)
inline int lattice_mask(double Sx, double Sy, const P2& c) {
    auto on = [](double v) { return std::fabs(v - std::round(v - 0.5) - 0.5) < 1e-7; };
    int m = 0;
    for (int ax = 0; ax < 2; ax++) for (int ay = 0; ay < 2; ay++) {
        const double fx = ax ? (Sx - c.x) : c.x, fy = ay ? (Sy - c.y) : c.y;
        if (on(fx) && on(fy)) m |= 1 << (ax + 2 * ay);
    }
    return m;
}

// candidate upright centres inside the box that do not overlap any tilted square by more than tol
inline std::vector<P2> candidates(double Sx, double Sy, const std::vector<Square>& T, double tol, int depth = 4) {
    std::set<std::pair<long long, long long>> seen;
    std::vector<P2> out;
    auto add = [&](double x, double y) {
        if (x < 0.5 - 1e-9 || y < 0.5 - 1e-9 || x > Sx - 0.5 + 1e-9 || y > Sy - 0.5 + 1e-9) return;
        x = std::min(std::max(x, 0.5), Sx - 0.5); y = std::min(std::max(y, 0.5), Sy - 0.5);
        const auto key = std::make_pair((long long)std::llround(x * 1e8), (long long)std::llround(y * 1e8));
        if (!seen.insert(key).second) return;
        for (const Square& t : T) {
            const double dx = t.x - x, dy = t.y - y;
            if (dx * dx + dy * dy < 2.0 && strip::pen(x, y, 0.0, t.x, t.y, t.angle) > tol) return;
        }
        out.push_back({x, y});
    };
    const int KX = (int)std::ceil(Sx), KY = (int)std::ceil(Sy);
    // four wall-anchored lattices
    std::vector<double> xs, ys;
    for (int i = 0; i < KX; i++) { xs.push_back(i + 0.5); xs.push_back(Sx - 0.5 - i); }
    for (int j = 0; j < KY; j++) { ys.push_back(j + 0.5); ys.push_back(Sy - 0.5 - j); }
    std::sort(xs.begin(), xs.end()); xs.erase(std::unique(xs.begin(), xs.end(), [](double a, double b) { return std::fabs(a - b) < 1e-9; }), xs.end());
    std::sort(ys.begin(), ys.end()); ys.erase(std::unique(ys.begin(), ys.end(), [](double a, double b) { return std::fabs(a - b) < 1e-9; }), ys.end());
    // pure lattices: x from the anchor family, y from the anchor family: all four combinations are the cross products of {i+.5} / {S-.5-i}
    for (int ax = 0; ax < 2; ax++) for (int ay = 0; ay < 2; ay++)
        for (int i = 0; i < KX; i++) for (int j = 0; j < KY; j++) add(ax ? Sx - 0.5 - i : i + 0.5, ay ? Sy - 0.5 - j : j + 0.5);
    // flush candidates near each tilted square
    for (const Square& t : T) {
        const double c = std::cos(t.angle), s = std::sin(t.angle);
        const double hx = 0.5 * (std::fabs(c) + std::fabs(s));
        const double xr = t.x + hx + 0.5, xl = t.x - hx - 0.5, yt = t.y + hx + 0.5, yb = t.y - hx - 0.5;
        for (int k = 0; k < depth; k++) {
            for (double xf : {xr + k, xl - k})
                for (double y : ys) if (std::fabs(y - t.y) < 3.6) add(xf, y);
            for (double yf : {yt + k, yb - k})
                for (double x : xs) if (std::fabs(x - t.x) < 3.6) add(x, yf);
        }
        // doubly flush corners
        for (double xf : {xr, xl}) for (double yf : {yt, yb}) add(xf, yf);
    }
    return out;
}

// ILS for a maximum independent set among axis-aligned unit squares (conflict iff |dx|<1-tol and |dy|<1-tol).
inline std::vector<int> mis_axis(const std::vector<P2>& C, double tol, int iters, std::mt19937_64& rng, const std::vector<std::vector<int>>* starts = nullptr) {
    const int N = (int)C.size();
    if (N == 0) return {};
    std::vector<std::vector<int>> adj(N);
    // bucket grid, cell 1
    double mx = 0, my = 0;
    for (auto& c : C) { mx = std::max(mx, c.x); my = std::max(my, c.y); }
    const int gx = (int)mx + 2, gy = (int)my + 2;
    std::vector<std::vector<int>> cell((size_t)gx * gy);
    for (int i = 0; i < N; i++) cell[(size_t)(int)C[i].y * gx + (int)C[i].x].push_back(i);
    for (int i = 0; i < N; i++) {
        const int cx = (int)C[i].x, cy = (int)C[i].y;
        for (int oy = -1; oy <= 1; oy++) for (int ox = -1; ox <= 1; ox++) {
            const int ax = cx + ox, ay = cy + oy;
            if (ax < 0 || ay < 0 || ax >= gx || ay >= gy) continue;
            for (int j : cell[(size_t)ay * gx + ax]) if (j > i && std::fabs(C[i].x - C[j].x) < 1 - tol && std::fabs(C[i].y - C[j].y) < 1 - tol) { adj[i].push_back(j); adj[j].push_back(i); }
        }
    }
    std::vector<char> in(N, 0);
    std::vector<int> tight(N, 0);
    int size = 0;
    auto insert = [&](int v) { in[v] = 1; size++; for (int u : adj[v]) tight[u]++; };
    auto remove = [&](int v) { in[v] = 0; size--; for (int u : adj[v]) tight[u]--; };
    std::vector<int> order(N);
    for (int i = 0; i < N; i++) order[i] = i;
    std::vector<int> best_set;
    int best = -1;
    auto record = [&]() { if (size > best) { best = size; best_set.clear(); for (int i = 0; i < N; i++) if (in[i]) best_set.push_back(i); } };
    auto reset = [&]() { std::fill(in.begin(), in.end(), 0); std::fill(tight.begin(), tight.end(), 0); size = 0; };
    auto greedy_fill = [&]() {
        std::shuffle(order.begin(), order.end(), rng);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return adj[a].size() < adj[b].size(); });
        for (int v : order) if (!in[v] && tight[v] == 0) insert(v);
    };
    auto local = [&]() {
        bool improved = true;
        while (improved) {
            improved = false;
            for (int v = 0; v < N; v++) if (!in[v] && tight[v] == 0) { insert(v); improved = true; }
            for (int v = 0; v < N; v++) {
                if (!in[v]) continue;
                std::vector<int> L;
                for (int u : adj[v]) if (!in[u] && tight[u] == 1) L.push_back(u);
                if (L.size() < 2) continue;
                bool done = false;
                for (size_t a = 0; a < L.size() && !done; a++)
                    for (size_t b = a + 1; b < L.size() && !done; b++)
                        if (std::find(adj[L[a]].begin(), adj[L[a]].end(), L[b]) == adj[L[a]].end()) { remove(v); insert(L[a]); insert(L[b]); done = true; improved = true; }
            }
        }
    };
    // starts: the caller's preferred partial solutions (kept in order, conflicts skipped), completed greedily, locally improved
    std::vector<char> best_in;
    if (starts) {
        for (const auto& st : *starts) {
            reset();
            for (int v : st) if (!in[v] && tight[v] == 0) insert(v);
            greedy_fill(); local();
            if (size > best) { record(); best_in = in; }
        }
    }
    if (best_in.empty()) { reset(); greedy_fill(); local(); record(); best_in = in; }
    else { reset(); for (int i = 0; i < N; i++) if (best_in[i]) insert(i); }
    std::vector<char> saved;
    for (int it = 0; it < iters; it++) {
        saved = in;
        const int before = size;
        const int k = (rng() % 4 == 0) ? 2 + (int)(rng() % 3) : 1;
        for (int q = 0; q < k; q++) {
            const int v = (int)(rng() % N);
            if (in[v]) continue;
            std::vector<int> rm;
            for (int u : adj[v]) if (in[u]) rm.push_back(u);
            for (int u : rm) remove(u);
            insert(v);
        }
        local();
        if (size >= before) record();
        else if (rng() % 100 >= 3) { std::fill(tight.begin(), tight.end(), 0); size = 0; for (int i = 0; i < N; i++) in[i] = 0; for (int i = 0; i < N; i++) if (saved[i]) insert(i); }
    }
    return best_set;
}

inline Result upright_fill(double Sx, double Sy, const std::vector<Square>& T, double tol, int iters, std::mt19937_64& rng) {
    Result R;
    const std::vector<P2> C = candidates(Sx, Sy, T, tol);
    R.candidates = (int)C.size();
    std::vector<int> mask(C.size());
    for (size_t i = 0; i < C.size(); i++) mask[i] = lattice_mask(Sx, Sy, C[i]);
    // PCA axis of the tilted core
    double cx = 0, cy = 0;
    for (const Square& t : T) { cx += t.x; cy += t.y; }
    if (!T.empty()) { cx /= T.size(); cy /= T.size(); } else { cx = Sx / 2; cy = Sy / 2; }
    double sxx = 0, sxy = 0, syy = 0;
    for (const Square& t : T) { sxx += (t.x - cx) * (t.x - cx); sxy += (t.x - cx) * (t.y - cy); syy += (t.y - cy) * (t.y - cy); }
    const double ang = 0.5 * std::atan2(2 * sxy, sxx - syy), ux = std::cos(ang), uy = std::sin(ang);
    std::vector<std::vector<int>> starts;
    for (int l = 0; l < 4; l++) { std::vector<int> st; for (size_t i = 0; i < C.size(); i++) if (mask[i] >> l & 1) st.push_back((int)i); starts.push_back(st); }
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) if (a != b) {      // lattice a on the left of the core's axis, b on the right
        std::vector<int> st;
        for (size_t i = 0; i < C.size(); i++) {
            const double side = -(C[i].x - cx) * uy + (C[i].y - cy) * ux;
            if (mask[i] >> (side >= 0 ? a : b) & 1) st.push_back((int)i);
        }
        starts.push_back(st);
    }
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) if (a != b) {      // same split through the perpendicular direction (cut across the strip)
        std::vector<int> st;
        for (size_t i = 0; i < C.size(); i++) {
            const double side = (C[i].x - cx) * ux + (C[i].y - cy) * uy;
            if (mask[i] >> (side >= 0 ? a : b) & 1) st.push_back((int)i);
        }
        starts.push_back(st);
    }
    {   // quadrant split around the centroid: UL (x<cx,y>=cy) from the (left, top) lattice, etc.
        std::vector<int> st;
        for (size_t i = 0; i < C.size(); i++) {
            const int l = (C[i].x >= cx ? 1 : 0) + 2 * (C[i].y >= cy ? 1 : 0);
            if (mask[i] >> l & 1) st.push_back((int)i);
        }
        starts.push_back(st);
    }
    for (int i : mis_axis(C, tol, iters, rng, &starts)) R.up.push_back(C[i]);
    return R;
}

}  // namespace fill
