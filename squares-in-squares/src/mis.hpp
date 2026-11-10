// Max-independent-set composer: candidate squares from corner-anchored axis lattices and tilted lattices,
// conflict = penetration deeper than `tol`; an iterated local search (ARW-style (1,2)-swaps + perturbation)
// picks the largest conflict-free subset.  The subset is a *start* for squeeze::relax, not a certified packing.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "strip.hpp"

namespace mis {

struct Cand { double x, y, t; };

struct Graph {
    std::vector<Cand> c;
    std::vector<std::vector<int>> adj;
};

// Candidate lattice: tilted by theta, reference point (ox, oy). Adds all squares fully inside [0,S]^2.
inline void add_lattice(std::vector<Cand>& out, double S, double theta, double ox, double oy) {
    const double ct = std::cos(theta), st = std::sin(theta), hh = 0.5 * (std::fabs(ct) + std::fabs(st));
    const int R = (int)std::ceil(S * 1.5) + 2;
    for (int a = -R; a <= R; a++)
        for (int b = -R; b <= R; b++) {
            const double x = ox + a * ct - b * st, y = oy + a * st + b * ct;
            if (x < hh || y < hh || x > S - hh || y > S - hh) continue;
            out.push_back({x, y, theta});
        }
}

inline void add_axis_lattices(std::vector<Cand>& out, double S, bool four = true) {
    const int K = (int)std::floor(S);
    for (int ax = 0; ax < (four ? 2 : 1); ax++)
        for (int ay = 0; ay < 2; ay++)
            for (int i = 0; i < K; i++)
                for (int j = 0; j < K; j++) {
                    const double x = ax ? S - i - 0.5 : i + 0.5, y = ay ? S - j - 0.5 : j + 0.5;
                    bool dup = false;
                    if (std::fabs(S - K) < 1e-9) dup = (ax || ay);
                    if (!dup) out.push_back({x, y, 0.0});
                }
}

inline Graph build(std::vector<Cand> cands, double tol) {
    Graph g;
    g.c = std::move(cands);
    const int N = (int)g.c.size();
    g.adj.assign(N, {});
    // bucket by 1.5 cells
    double maxc = 0;
    for (auto& q : g.c) maxc = std::max(maxc, std::max(q.x, q.y));
    const int gx = (int)(maxc / 1.5) + 2;
    std::vector<std::vector<int>> cell(gx * gx);
    for (int i = 0; i < N; i++) cell[(int)(g.c[i].y / 1.5) * gx + (int)(g.c[i].x / 1.5)].push_back(i);
    for (int i = 0; i < N; i++) {
        const int cx = (int)(g.c[i].x / 1.5), cy = (int)(g.c[i].y / 1.5);
        for (int oy = -1; oy <= 1; oy++)
            for (int ox = -1; ox <= 1; ox++) {
                const int ax = cx + ox, ay = cy + oy;
                if (ax < 0 || ay < 0 || ax >= gx || ay >= gx) continue;
                for (int j : cell[ay * gx + ax]) {
                    if (j <= i) continue;
                    if (strip::pen(g.c[i].x, g.c[i].y, g.c[i].t, g.c[j].x, g.c[j].y, g.c[j].t) > tol) {
                        g.adj[i].push_back(j); g.adj[j].push_back(i);
                    }
                }
            }
    }
    return g;
}

// Iterated local search for a large independent set. Returns the best set found within `iters` perturbations.
inline std::vector<int> solve(const Graph& g, int iters, std::mt19937_64& rng) {
    const int N = (int)g.c.size();
    std::vector<char> in(N, 0);
    std::vector<int> tight(N, 0);   // number of neighbours inside the set
    int size = 0;
    auto insert = [&](int v) { in[v] = 1; size++; for (int u : g.adj[v]) tight[u]++; };
    auto remove = [&](int v) { in[v] = 0; size--; for (int u : g.adj[v]) tight[u]--; };
    // greedy by min degree with random tie-break
    std::vector<int> order(N);
    for (int i = 0; i < N; i++) order[i] = i;
    std::shuffle(order.begin(), order.end(), rng);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return g.adj[a].size() < g.adj[b].size(); });
    for (int v : order) if (!in[v] && tight[v] == 0) insert(v);
    std::vector<int> best_set;
    int best = -1;
    auto record = [&]() {
        if (size > best) { best = size; best_set.clear(); for (int i = 0; i < N; i++) if (in[i]) best_set.push_back(i); }
    };
    // (1,2)-swap local search: for each solution vertex v, find two non-adjacent 1-tight neighbours
    auto local = [&]() {
        bool improved = true;
        while (improved) {
            improved = false;
            for (int v = 0; v < N; v++) {   // free vertices first (tight 0)
                if (!in[v] && tight[v] == 0) { insert(v); improved = true; }
            }
            for (int v = 0; v < N; v++) {
                if (!in[v]) continue;
                std::vector<int> L;
                for (int u : g.adj[v]) if (!in[u] && tight[u] == 1) L.push_back(u);
                if (L.size() < 2) continue;
                bool done = false;
                for (size_t a = 0; a < L.size() && !done; a++)
                    for (size_t b = a + 1; b < L.size() && !done; b++) {
                        const int x = L[a], y = L[b];
                        if (std::find(g.adj[x].begin(), g.adj[x].end(), y) == g.adj[x].end()) {
                            remove(v); insert(x); insert(y); done = true; improved = true;
                        }
                    }
            }
        }
    };
    local();
    record();
    std::vector<char> saved;
    for (int it = 0; it < iters; it++) {
        saved = in;
        const int before = size;
        // perturb: force-insert k random non-solution vertices
        const int k = (rng() % 4 == 0) ? 2 + (int)(rng() % 3) : 1;
        for (int q = 0; q < k; q++) {
            int v = (int)(rng() % N);
            if (in[v]) continue;
            std::vector<int> rm;
            for (int u : g.adj[v]) if (in[u]) rm.push_back(u);
            for (int u : rm) remove(u);
            insert(v);
        }
        local();
        if (size >= before) { record(); }
        else if (rng() % 100 >= 3) {   // reject: restore
            std::fill(tight.begin(), tight.end(), 0);
            size = 0;
            for (int i = 0; i < N; i++) in[i] = 0;
            for (int i = 0; i < N; i++) if (saved[i]) insert(i);
        }
    }
    return best_set;
}

}  // namespace mis
