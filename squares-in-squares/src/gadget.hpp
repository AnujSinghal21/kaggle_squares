// Fault-line gadget DP (own implementation of the composition described in Chris Deotte's "multi-gadget fault DP" notebook;
// design cross-checked against the public description in alex chilton's starter kit).
//
//   n = K^2 + D with K = floor(sqrt n).  The lower-left continent is the integer lattice, the upper-right continent is the
//   same lattice shifted by F in x and y.  Along the staircase fault between them sits a chain of "gadgets": a w x h cell
//   block that holds w*h + dm squares inside [0, w+f] x [0, h+f].  A chain adds D = sum dm squares and the box side is K + F.
//   Gadgets are crops of donor packings (any size, 8 dihedral maps): of the squares of a donor under a map, keep the
//   w*h + dm whose overhang past the w x h window is smallest; f is the largest overhang among them.
//   Seams between consecutive gadgets: mode 0 (corners touch, guarded by the tip rule), mode 1 (one lattice column between),
//   mode 2 (one lattice row between).  A DP over (columns used, rows used, surplus, tip bin) minimises the largest f.
#pragma once
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <numeric>
#include <thread>
#include <unordered_map>
#include <vector>

#include "common.hpp"

namespace gadget {

constexpr int kWMax = 31, kDMax = 62;
constexpr double kFMax = 0.97, kE = 1e-9;
constexpr int kNB = 6;
constexpr double kThr[kNB] = {0.0, 0.3, 0.45, 0.6, 0.75, 0.9};

inline int bin_of(double v) { int b = 0; for (int i = 0; i < kNB; i++) if (v >= kThr[i]) b = i; return b; }

struct Donor { int n = 0; double L = 0; std::vector<double> x, y, t; };

struct Gadget {
    int w, h, dm;
    double f, ybr, xtl;
    int donor, sym;
    int qb, qt;
    std::vector<Square> sq;   // explicit squares (quenched gadgets); empty = crop from donor/sym
};

// 8 dihedral maps: (sx, sy, transpose)
inline void map_donor(const Donor& d, int sym, std::vector<double>& xlo, std::vector<double>& xhi, std::vector<double>& ylo,
                      std::vector<double>& yhi, std::vector<double>* mx = nullptr, std::vector<double>* my = nullptr, std::vector<double>* mt = nullptr) {
    const int sx = (sym & 1) ? -1 : 1, sy = (sym & 2) ? -1 : 1, tr = (sym & 4) ? 1 : 0;
    const int n = d.n;
    std::vector<double> x(n), y(n), t(n), he(n);
    for (int i = 0; i < n; i++) {
        x[i] = d.x[i] * sx; y[i] = d.y[i] * sy; t[i] = d.t[i] * sx * sy;
        if (tr) { std::swap(x[i], y[i]); t[i] = -t[i]; }
        he[i] = 0.5 * (std::fabs(std::cos(t[i])) + std::fabs(std::sin(t[i])));
    }
    double mnx = 1e300, mny = 1e300;
    for (int i = 0; i < n; i++) { mnx = std::min(mnx, x[i] - he[i]); mny = std::min(mny, y[i] - he[i]); }
    xlo.resize(n); xhi.resize(n); ylo.resize(n); yhi.resize(n);
    for (int i = 0; i < n; i++) {
        x[i] -= mnx; y[i] -= mny;
        xlo[i] = x[i] - he[i]; xhi[i] = x[i] + he[i]; ylo[i] = y[i] - he[i]; yhi[i] = y[i] + he[i];
    }
    if (mx) { *mx = x; *my = y; *mt = t; }
}

struct Key { uint64_t v; };

// Extract the best gadget per (w, h, dm, qb, qt) over all donors (multi-threaded).
inline std::vector<Gadget> extract(const std::vector<Donor>& donors, int threads, int wmax = kWMax) {
    std::mutex mu;
    std::unordered_map<uint64_t, Gadget> best;
    std::atomic<size_t> next{0};
    auto worker = [&]() {
        std::unordered_map<uint64_t, Gadget> local;
        std::vector<double> xlo, xhi, ylo, yhi, R;
        std::vector<int> idx;
        for (;;) {
            const size_t di = next++;
            if (di >= donors.size()) break;
            const Donor& d = donors[di];
            const int n = d.n, Lc = (int)std::ceil(d.L - 1e-9);
            const int wm = std::min(wmax, Lc);
            for (int sym = 0; sym < 8; sym++) {
                map_donor(d, sym, xlo, xhi, ylo, yhi);
                R.resize(n); idx.resize(n);
                for (int h = 1; h <= wm; h++)
                    for (int w = 1; w <= wm; w++) {
                        if (w * h + 1 > n) continue;
                        for (int i = 0; i < n; i++) { R[i] = std::max(std::max(xhi[i] - w, yhi[i] - h), 0.0); idx[i] = i; }
                        const int need = std::min(n, w * h + kDMax);
                        std::nth_element(idx.begin(), idx.begin() + (need - 1), idx.end(), [&](int a, int b) { return R[a] < R[b] || (R[a] == R[b] && a < b); });
                        std::sort(idx.begin(), idx.begin() + need, [&](int a, int b) { return R[a] < R[b] || (R[a] == R[b] && a < b); });
                        double ybr = 9.0, xtl = 9.0;
                        for (int k = 0; k < need; k++) {
                            const int i = idx[k];
                            if (xhi[i] > w + 1e-9) ybr = std::min(ybr, ylo[i]);
                            if (yhi[i] > h + 1e-9) xtl = std::min(xtl, xlo[i]);
                            const int dm = k + 1 - w * h;
                            if (dm < 1) continue;
                            if (dm > kDMax) break;
                            const double f = R[i];   // sorted ascending: f is the largest overhang among the first k+1
                            if (f > kFMax) break;
                            Gadget g{w, h, dm, f, std::min(ybr, 9.0), std::min(xtl, 9.0), (int)di, sym, bin_of(std::min(ybr, 9.0)), bin_of(std::min(xtl, 9.0))};
                            const uint64_t key = ((((uint64_t)w * 32 + h) * 64 + dm) * kNB + g.qb) * kNB + g.qt;
                            auto it = local.find(key);
                            if (it == local.end() || f < it->second.f) local[key] = g;
                        }
                    }
            }
        }
        std::lock_guard<std::mutex> lk(mu);
        for (auto& kv : local) { auto it = best.find(kv.first); if (it == best.end() || kv.second.f < it->second.f) best[kv.first] = kv.second; }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    std::vector<Gadget> G;
    G.reserve(best.size());
    for (auto& kv : best) G.push_back(kv.second);
    // dominance: gadget g is dropped if some gadget with the same (w, h, dm) or a "bigger" tip bins pattern has smaller-or-equal f.
    // Tip bins are better when larger (a larger free tip length makes seams easier): (qb, qt) dominated by (qb' >= qb, qt' >= qt, f' <= f).
    std::unordered_map<uint64_t, std::vector<int>> groups;
    for (int i = 0; i < (int)G.size(); i++) groups[((uint64_t)G[i].w * 32 + G[i].h) * 64 + G[i].dm].push_back(i);
    std::vector<char> keep(G.size(), 1);
    for (auto& kv : groups) {
        auto& v = kv.second;
        for (int a : v)
            for (int b : v)
                if (a != b && keep[b] && G[b].f <= G[a].f + 1e-15 && G[b].qb >= G[a].qb && G[b].qt >= G[a].qt &&
                    (G[b].f < G[a].f - 1e-15 || G[b].qb > G[a].qb || G[b].qt > G[a].qt || b < a)) { keep[a] = 0; break; }
    }
    // also larger (w, h) are NOT dominated by smaller ones (they use more fault); but for equal w,h larger dm is not dominated either.
    std::vector<Gadget> out;
    for (size_t i = 0; i < G.size(); i++) if (keep[i]) out.push_back(G[i]);
    return out;
}

struct Dp {
    // flat arrays dims [cw 0..W][ch 0..W][d 0..DM][q 0..NB-1]
    static constexpr int W1 = kWMax + 1, D1 = kDMax + 1;
    std::vector<double> val;
    std::vector<int32_t> pg;
    std::vector<int8_t> pm, pq;
    size_t at(int cw, int ch, int d, int q) const { return (((size_t)cw * W1 + ch) * D1 + d) * kNB + q; }
};

inline Dp run_dp(const std::vector<Gadget>& G, int threads) {
    Dp dp;
    const size_t N = (size_t)Dp::W1 * Dp::W1 * Dp::D1 * kNB;
    dp.val.assign(N, INFINITY); dp.pg.assign(N, -1); dp.pm.assign(N, 0); dp.pq.assign(N, 0);
    dp.val[dp.at(0, 0, 0, kNB - 1)] = 0;
    // process columns in increasing order; gadgets of one column layer are independent targets -> parallelise over gadgets with locks per target column
    for (int cw = 0; cw <= kWMax; cw++) {
        // is there anything finite in this column?
        bool any = false;
        for (int ch = 0; ch <= kWMax && !any; ch++) for (int d = 0; d <= kDMax && !any; d++) for (int q = 0; q < kNB; q++) if (std::isfinite(dp.val[dp.at(cw, ch, d, q)])) { any = true; break; }
        if (!any) continue;
        std::mutex mu;
        std::atomic<int> next{0};
        auto worker = [&]() {
            for (;;) {
                const int g = next++;
                if (g >= (int)G.size()) return;
                const Gadget& gd = G[g];
                const int w = gd.w, h = gd.h, dm = gd.dm;
                if (cw + w > kWMax) continue;
                // local updates collected then merged under lock
                struct Upd { size_t pos; double v; int32_t g; int8_t m, qp; };
                std::vector<Upd> ups;
                for (int ch = 0; ch + h <= kWMax; ch++)
                    for (int d = 0; d + dm <= kDMax; d++) {
                        // best over previous q for each mode
                        for (int mode = 0; mode < 3; mode++) {
                            if (cw == 0 && mode != 0) continue;
                            int tcw = cw + w + (mode == 1), tch = ch + h + (mode == 2);
                            if (tcw > kWMax || tch > kWMax) continue;
                            double bestv = INFINITY; int bq = 0;
                            for (int q = 0; q < kNB; q++) {
                                const double cur = dp.val[dp.at(cw, ch, d, q)];
                                if (!std::isfinite(cur)) continue;
                                if (mode == 0) {
                                    const bool ok0 = (kThr[q] >= gd.f + kE) || (gd.xtl >= cur + kE);
                                    if (!ok0) continue;
                                }
                                const double v = std::max(cur, gd.f);
                                if (v < bestv) { bestv = v; bq = q; }
                            }
                            if (std::isfinite(bestv)) ups.push_back({dp.at(tcw, tch, d + dm, gd.qb), bestv, (int32_t)g, (int8_t)mode, (int8_t)bq});
                        }
                    }
                std::lock_guard<std::mutex> lk(mu);
                for (const Upd& u : ups) if (u.v < dp.val[u.pos]) { dp.val[u.pos] = u.v; dp.pg[u.pos] = u.g; dp.pm[u.pos] = u.m; dp.pq[u.pos] = u.qp; }
            }
        };
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; t++) pool.emplace_back(worker);
        for (auto& t : pool) t.join();
    }
    return dp;
}

struct Link { int g; int mode; };

// Best chain for (K, D): min over cw, ch <= K and q of dp[cw][ch][D][q]; returns f or -1.
inline double best_chain(const Dp& dp, const std::vector<Gadget>& G, int K, int D, std::vector<Link>& links) {
    if (D > kDMax || K > kWMax) return -1;
    double bestv = INFINITY; int bcw = -1, bch = -1, bq = -1;
    for (int cw = 0; cw <= K; cw++) for (int ch = 0; ch <= K; ch++) for (int q = 0; q < kNB; q++) {
        const double v = dp.val[dp.at(cw, ch, D, q)];
        if (v < bestv) { bestv = v; bcw = cw; bch = ch; bq = q; }
    }
    if (!std::isfinite(bestv) || bcw < 0) return -1;
    links.clear();
    int cw = bcw, ch = bch, d = D, q = bq;
    while (cw || ch || d) {
        const size_t p = dp.at(cw, ch, d, q);
        const int g = dp.pg[p], mode = dp.pm[p], qp = dp.pq[p];
        if (g < 0) return -1;
        links.push_back({g, mode});
        cw -= G[g].w + (mode == 1); ch -= G[g].h + (mode == 2); d -= G[g].dm; q = qp;
    }
    std::reverse(links.begin(), links.end());
    return bestv;
}

// Crop: the w*h+dm squares with the smallest overhang (stable order); returns positions in the mapped frame.
inline void crop(const Donor& d, int sym, int w, int h, int dm, std::vector<Square>& out, double* fmax) {
    std::vector<double> xlo, xhi, ylo, yhi, x, y, t;
    map_donor(d, sym, xlo, xhi, ylo, yhi, &x, &y, &t);
    const int n = d.n;
    std::vector<double> R(n);
    std::vector<int> idx(n);
    for (int i = 0; i < n; i++) { R[i] = std::max(std::max(xhi[i] - w, yhi[i] - h), 0.0); idx[i] = i; }
    std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return R[a] < R[b]; });
    out.clear();
    double fm = 0;
    for (int k = 0; k < w * h + dm; k++) { out.push_back({x[idx[k]], y[idx[k]], t[idx[k]]}); fm = std::max(fm, R[idx[k]]); }
    *fmax = fm;
}

// Compose the packing for (K, D) from a chain; returns empty on failure.  Positions follow the DP's geometry:
// gadget m sits at column ox, row oy (its lower-left), the LL lattice fills below a[i], the UR lattice (shifted by F) above b[i].
inline Packing compose(int K, int D, double f, const std::vector<Link>& links, const std::vector<Gadget>& G, const std::vector<Donor>& donors) {
    const int M = (int)links.size();
    const double F = f + (M + 2) * kE;
    std::vector<int> a(K, K), b(K, K);
    int ox = 0, oy = K;
    std::vector<std::vector<Square>> gad;
    for (int m = 0; m < M; m++) {
        const Gadget& g = G[links[m].g];
        const int mode = links[m].mode;
        if (m) { ox += (mode == 1); oy -= (mode == 2); }
        oy -= g.h;
        if (oy < 0 || ox + g.w > K) return {};
        if (mode == 1 && ox - 1 >= 0) a[ox - 1] = b[ox - 1] = oy + g.h;
        for (int i = ox; i < ox + g.w; i++) { a[i] = oy; b[i] = oy + g.h; }
        std::vector<Square> sq;
        double fg = 0;
        if (g.sq.empty() && g.donor < 0) return {};   // hypothetical (what-if) gadget
        if (!g.sq.empty()) { sq = g.sq; fg = g.f; }
        else crop(donors[g.donor], g.sym, g.w, g.h, g.dm, sq, &fg);
        if (fg > f + 1e-12) return {};
        for (Square& s : sq) { s.x *= 1 + 1e-11; s.y *= 1 + 1e-11; s.x += ox + kE * (1 + m); s.y += oy + kE * (M - m); }
        gad.push_back(std::move(sq));
        ox += g.w;
    }
    for (int i = ox; i < K; i++) a[i] = b[i] = 0;
    for (int i = 0; i + 1 < K; i++) if (a[i + 1] > b[i]) return {};
    for (int i = 0; i < K; i++) if (a[i] > b[i]) return {};
    Packing p;
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < a[i]; j++) p.push_back({i + 0.5, j + 0.5, 0.0});
        for (int j = b[i]; j < K; j++) p.push_back({i + 0.5 + F, j + 0.5 + F, 0.0});
    }
    for (auto& g : gad) for (const Square& s : g) p.push_back(s);
    if ((int)p.size() != K * K + D) return {};
    for (Square& q : p) { q.x *= 1 + 2e-11; q.y *= 1 + 2e-11; }   // open the exact lattice contacts into real gaps (float64 round trip)
    return p;
}

}  // namespace gadget
