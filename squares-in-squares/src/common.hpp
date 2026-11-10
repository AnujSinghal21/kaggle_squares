// Shared types and the (cheap) score maths for Squares in Squares.
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

constexpr int kMaxN = 1000;

struct Square {
    double x, y;   // center
    double angle;  // radians
};
using Packing = std::vector<Square>;

// Side length of the tightest axis-aligned container around the packing.
// A unit square rotated by a has half-extent 0.5 * (|cos a| + |sin a|) on both axes.
inline double container_side(const Packing& p) {
    double lo_x = INFINITY, hi_x = -INFINITY, lo_y = INFINITY, hi_y = -INFINITY;
    for (const Square& s : p) {
        const double h = 0.5 * (std::fabs(std::cos(s.angle)) + std::fabs(std::sin(s.angle)));
        lo_x = std::min(lo_x, s.x - h);
        hi_x = std::max(hi_x, s.x + h);
        lo_y = std::min(lo_y, s.y - h);
        hi_y = std::max(hi_y, s.y + h);
    }
    return std::max(hi_x - lo_x, hi_y - lo_y);
}

// score(n) = L / sqrt(n) - 1, lower is better.
inline double score_of(int n, double L) { return L / std::sqrt((double)n) - 1.0; }

// What a missing or invalid puzzle is charged: the trivial grid ceiling.
inline double grid_score(int n) { return score_of(n, std::ceil(std::sqrt((double)n))); }

// Baseline: ceil(sqrt(n)) x ceil(sqrt(n)) axis-aligned grid, filled row by row.
inline Packing grid_packing(int n) {
    const int k = (int)std::ceil(std::sqrt((double)n));
    Packing p;
    p.reserve(n);
    for (int i = 0; i < n; i++) p.push_back({i % k + 0.5, i / k + 0.5, 0.0});
    return p;
}

inline double score(const Packing& p) { return score_of((int)p.size(), container_side(p)); }

// ---- overlap checking (separating axis theorem), mirrors stage 1 of the Kaggle metric ----
// Pairs within +-kSatTol of touching go to a long double re-check: exact contacts (overlap == 0,
// e.g. the upright grid) pass, anything else in that band is "tight". Kaggle certifies those with
// interval arithmetic and a rounded irrational contact (rotated squares) typically fails there,
// so leave a real gap (e.g. dilate centers by ~1+1e-7).
constexpr double kSatTol = 1e-8;
// Inside that band, the long double overlap is trusted when it is clearly negative (a gap bigger
// than kCertGap survives the float64 round trip and Kaggle's interval check); exactly 0 is a clean
// contact; anything in between is "tight" = risky.
constexpr double kCertGap = 1e-15;

struct Prepared {  // per-square data for the pair test
    double cx, cy, c, s;
    double px[4], py[4];
};

inline Prepared prepare(const Square& q) {
    Prepared p;
    p.cx = q.x;
    p.cy = q.y;
    p.c = std::cos(q.angle);
    p.s = std::sin(q.angle);
    int k = 0;
    for (int sx = -1; sx <= 1; sx += 2)
        for (int sy = -1; sy <= 1; sy += 2) {
            p.px[k] = q.x + 0.5 * (p.c * sx - p.s * sy);
            p.py[k] = q.y + 0.5 * (p.s * sx + p.c * sy);
            k++;
        }
    return p;
}

// Smallest axis overlap over the 4 SAT axes: < 0 means separated by that gap, > 0 means penetrating.
inline double min_axis_overlap(const Prepared& a, const Prepared& b) {
    const double ax[4] = {a.c, -a.s, b.c, -b.s};
    const double ay[4] = {a.s, a.c, b.s, b.c};
    double best = INFINITY;
    for (int k = 0; k < 4; k++) {
        double a_lo = INFINITY, a_hi = -INFINITY, b_lo = INFINITY, b_hi = -INFINITY;
        for (int i = 0; i < 4; i++) {
            const double pa = a.px[i] * ax[k] + a.py[i] * ay[k];
            const double pb = b.px[i] * ax[k] + b.py[i] * ay[k];
            a_lo = std::min(a_lo, pa);
            a_hi = std::max(a_hi, pa);
            b_lo = std::min(b_lo, pb);
            b_hi = std::max(b_hi, pb);
        }
        best = std::min(best, std::min(a_hi, b_hi) - std::max(a_lo, b_lo));
    }
    return best;
}

inline double min_axis_overlap_ld(const Square& qa, const Square& qb) {
    using LD = long double;
    auto corners = [](const Square& q, LD* px, LD* py, LD& c, LD& s) {
        c = cosl((LD)q.angle);
        s = sinl((LD)q.angle);
        int k = 0;
        for (int sx = -1; sx <= 1; sx += 2)
            for (int sy = -1; sy <= 1; sy += 2) {
                px[k] = (LD)q.x + 0.5L * (c * sx - s * sy);
                py[k] = (LD)q.y + 0.5L * (s * sx + c * sy);
                k++;
            }
    };
    LD apx[4], apy[4], bpx[4], bpy[4], ac, as, bc, bs;
    corners(qa, apx, apy, ac, as);
    corners(qb, bpx, bpy, bc, bs);
    const LD ax[4] = {ac, -as, bc, -bs}, ay[4] = {as, ac, bs, bc};
    LD best = INFINITY;
    for (int k = 0; k < 4; k++) {
        LD a_lo = INFINITY, a_hi = -INFINITY, b_lo = INFINITY, b_hi = -INFINITY;
        for (int i = 0; i < 4; i++) {
            const LD pa = apx[i] * ax[k] + apy[i] * ay[k];
            const LD pb = bpx[i] * ax[k] + bpy[i] * ay[k];
            a_lo = std::min(a_lo, pa);
            a_hi = std::max(a_hi, pa);
            b_lo = std::min(b_lo, pb);
            b_hi = std::max(b_hi, pb);
        }
        best = std::min(best, std::min(a_hi, b_hi) - std::max(a_lo, b_lo));
    }
    return (double)best;
}

// Check one packing for overlaps. Returns the number of tight pairs, or -1 on a real overlap
// (offending pair in i, j).
inline long long check_overlaps(const Packing& p, int& bi, int& bj, double& worst) {
    const int n = (int)p.size();
    std::vector<Prepared> pre(n);
    for (int i = 0; i < n; i++) pre[i] = prepare(p[i]);
    long long tight = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            const double dx = pre[i].cx - pre[j].cx, dy = pre[i].cy - pre[j].cy;
            if (dx * dx + dy * dy > 2.0000001) continue;  // circumradius sqrt(2)/2 each: cannot touch
            const double ov = min_axis_overlap(pre[i], pre[j]);
            if (ov < -kSatTol) continue;
            const double ov_ld = ov > kSatTol ? ov : min_axis_overlap_ld(p[i], p[j]);
            if (ov_ld > kCertGap) {  // clearly penetrating
                bi = i;
                bj = j;
                worst = ov_ld;
                return -1;
            }
            if (ov_ld > -kCertGap && ov_ld != 0.0) tight++;  // exact 0 (upright grid) and real gaps pass
        }
    return tight;
}


// True if no two squares overlap. `tight` (optional) gets the number of risky near-contacts;
// `bi`/`bj` (optional) get the first offending pair when it returns false.
inline bool packing_valid(const Packing& p, long long* tight = nullptr, int* bi = nullptr, int* bj = nullptr) {
    int i = 0, j = 0;
    double worst = 0;
    const long long t = check_overlaps(p, i, j, worst);
    if (bi) *bi = i;
    if (bj) *bj = j;
    if (tight) *tight = t < 0 ? 0 : t;
    return t >= 0;
}
