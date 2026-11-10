// Closure pass, run by main.cpp after every solver so no approach has to remember it.
//   delete:  L(n) <= L(n+1)   a packing of n+1 minus its most useful square packs n squares.
//   add-L:   a packing of n in side L plus a column and a row of upright squares packs
//            n + 2*floor(L) + 1 squares in side L + 1.
//   tile:    a packing of n in side L, repeated k x k times, packs k^2 * n squares in side k*L
//            (same score as n, so it only helps if k^2*n is currently worse than n).
// Delete runs from n=1000 down to 1 (so one gain at n ripples to every smaller n in a single backward
// pass); add-L runs upward. The two alternate until nothing changes.
#pragma once
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "common.hpp"
#include "validate.hpp"

struct ClosureChange {
    int n;
    double old_L, new_L;
    std::string how;
};

namespace closure_detail {

constexpr double kMinGain = 1e-7;   // ignore dilation-level noise (real closure gains are far larger)
constexpr double kGap = 2e-11;      // real gaps (>> 1e-15) so the new squares are never "risky" contacts

struct Ext { double lo_x, hi_x, lo_y, hi_y; };

inline Ext extent_of(const Square& s) {
    const double h = 0.5 * (std::fabs(std::cos(s.angle)) + std::fabs(std::sin(s.angle)));
    return {s.x - h, s.x + h, s.y - h, s.y + h};
}

// Index of the square whose removal shrinks the container the most, and the resulting side.
inline int best_deletion(const Packing& p, double& side_out) {
    const int m = (int)p.size();
    std::vector<Ext> e(m);
    for (int i = 0; i < m; i++) e[i] = extent_of(p[i]);
    std::vector<Ext> pre(m + 1), suf(m + 1);
    const double inf = INFINITY;
    pre[0] = suf[m] = {inf, -inf, inf, -inf};
    for (int i = 0; i < m; i++)
        pre[i + 1] = {std::min(pre[i].lo_x, e[i].lo_x), std::max(pre[i].hi_x, e[i].hi_x),
                      std::min(pre[i].lo_y, e[i].lo_y), std::max(pre[i].hi_y, e[i].hi_y)};
    for (int i = m - 1; i >= 0; i--)
        suf[i] = {std::min(suf[i + 1].lo_x, e[i].lo_x), std::max(suf[i + 1].hi_x, e[i].hi_x),
                  std::min(suf[i + 1].lo_y, e[i].lo_y), std::max(suf[i + 1].hi_y, e[i].hi_y)};
    int best = 0;
    side_out = INFINITY;
    for (int i = 0; i < m; i++) {
        const double lx = std::min(pre[i].lo_x, suf[i + 1].lo_x), hx = std::max(pre[i].hi_x, suf[i + 1].hi_x);
        const double ly = std::min(pre[i].lo_y, suf[i + 1].lo_y), hy = std::max(pre[i].hi_y, suf[i + 1].hi_y);
        const double side = std::max(hx - lx, hy - ly);
        if (side < side_out) { side_out = side; best = i; }
    }
    return best;
}

// p plus an upright column on the right and row on top: 2*floor(L)+1 more squares.
inline Packing add_frame(const Packing& p, double L) {
    double x0 = INFINITY, y0 = INFINITY;
    for (const Square& s : p) { const Ext e = extent_of(s); x0 = std::min(x0, e.lo_x); y0 = std::min(y0, e.lo_y); }
    const int k = (int)std::floor(L + 1e-9);
    Packing q = p;
    for (int j = 0; j < k; j++) q.push_back({x0 + L + kGap + 0.5, y0 + 0.5 + j * (1.0 + kGap), 0.0});          // right column
    for (int i = 0; i <= k; i++) q.push_back({x0 + 0.5 + i * (1.0 + kGap), y0 + L + (k + 1) * kGap + 0.5, 0.0});  // top row
    return q;
}

// k x k copies of p, each in its own L x L cell. The cell pitch is a multiple of 2^-30 so that
// translating dyadic coordinates (e.g. an upright grid) stays exact and creates no fake near-contacts.
inline Packing tile_copies(const Packing& p, double L, int k, double& side_out) {
    double x0 = INFINITY, y0 = INFINITY;
    for (const Square& s : p) { const Ext e = extent_of(s); x0 = std::min(x0, e.lo_x); y0 = std::min(y0, e.lo_y); }
    const double pitch = std::ceil((L + kGap) * 1073741824.0) / 1073741824.0;
    Packing q;
    q.reserve(p.size() * k * k);
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            for (const Square& s : p) q.push_back({s.x - x0 + i * pitch, s.y - y0 + j * pitch, s.angle});
    side_out = k * pitch;
    return q;
}

}  // namespace closure_detail

// Applies both closures to every puzzle in `sub`. Returns the net changes (first old side -> final new side).
inline std::vector<ClosureChange> closure_pass(Submission& sub) {
    using namespace closure_detail;
    std::vector<double> L(kMaxN + 2, 0.0);
    std::vector<char> have(kMaxN + 2, 0);
    for (int n = 1; n <= kMaxN; n++) {
        have[n] = sub.present[n] && sub.broken[n].empty() && (int)sub.pk[n].size() == n;
        L[n] = have[n] ? container_side(sub.pk[n]) : std::ceil(std::sqrt((double)n));  // missing = grid ceiling
    }
    std::map<int, ClosureChange> net;
    auto adopt = [&](int n, Packing p, double new_L, const std::string& how) {
        auto it = net.find(n);
        if (it == net.end()) net[n] = {n, L[n], new_L, how};
        else { it->second.new_L = new_L; it->second.how = how; }
        sub.pk[n] = std::move(p);
        sub.present[n] = 1;
        sub.broken[n].clear();
        have[n] = 1;
        L[n] = new_L;
    };

    for (int round = 0; round < 10; round++) {
        bool changed = false;
        for (int n = 1; n <= kMaxN; n++) {  // add-L: n -> n + 2*floor(L) + 1
            if (!have[n]) continue;
            const int m = n + 2 * (int)std::floor(L[n] + 1e-9) + 1;
            if (m > kMaxN) continue;
            Packing q = add_frame(sub.pk[n], L[n]);
            const double Lq = container_side(q);
            if (Lq >= L[m] - kMinGain) continue;
            long long risky = 0;
            if (!packing_valid(q, &risky) || risky) continue;
            adopt(m, std::move(q), Lq, "add-L from n=" + std::to_string(n));
            changed = true;
        }
        for (int n = 1; n <= kMaxN; n++) {  // tile: n -> k^2 n
            if (!have[n]) continue;
            for (int k = 2; (long long)k * k * n <= kMaxN; k++) {
                const int m = k * k * n;
                if (k * L[n] + k * kGap >= L[m] - kMinGain) continue;  // cheap pre-check before building
                double Lq = 0;
                Packing q = tile_copies(sub.pk[n], L[n], k, Lq);
                if (Lq >= L[m] - kMinGain) continue;
                long long risky = 0;
                if (!packing_valid(q, &risky) || risky) continue;
                adopt(m, std::move(q), Lq, "tile " + std::to_string(k) + "x" + std::to_string(k) + " of n=" + std::to_string(n));
                changed = true;
            }
        }
        for (int n = kMaxN - 1; n >= 1; n--) {  // delete: n+1 -> n, single backward pass
            if (!have[n + 1]) continue;
            double Lc = 0;
            const int drop = best_deletion(sub.pk[n + 1], Lc);
            if (Lc >= L[n] - kMinGain) continue;
            Packing q = sub.pk[n + 1];
            q.erase(q.begin() + drop);
            adopt(n, std::move(q), Lc, "delete from n+1");
            changed = true;
        }
        if (!changed) break;
    }
    std::vector<ClosureChange> out;
    for (auto& kv : net) out.push_back(kv.second);
    return out;
}
