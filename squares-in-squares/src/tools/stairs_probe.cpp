// Local search over which corridor-boundary grid squares to keep, objective = kept + column-DP count.
// usage: stairs_probe S theta_deg a r [iters] [tol] [step]
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include "../slab.hpp"
int main(int argc, char** argv) {
    const double S = atof(argv[1]), th = atof(argv[2]) * M_PI / 180;
    const int a = atoi(argv[3]), r = atoi(argv[4]), iters = argc > 5 ? atoi(argv[5]) : 400;
    const double tol = argc > 6 ? atof(argv[6]) : 0.01, step = argc > 7 ? atof(argv[7]) : 0.02;
    const int k = (int)std::floor(S); const double f = S - k;
    const int b = a - 1 - r;
    struct Cand { slab::Sq q; bool boundary; bool keep; };
    std::vector<Cand> cand;
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) {
        if (j - i >= a - 3 && j - i < a + 1 + 0) {}                      // (kept below by rule; boundary flagged)
        const bool ul = (j - i >= a), lr = (j - i <= b);
        const bool ul_edge = (j - i >= a - 2 && j - i <= a + 2), lr_edge = (j - i <= b + 2 && j - i >= b - 2);
        if (ul) cand.push_back({{i + .5, f + j + .5, 0}, ul_edge, true});
        else if (lr) cand.push_back({{f + i + .5, j + .5, 0}, lr_edge, true});
        // also allow the adjacent unused diagonals of each lattice as optional additions
        if (!ul && j - i >= a - 3) cand.push_back({{i + .5, f + j + .5, 0}, true, false});
        if (!lr && j - i <= b + 3 && !ul) cand.push_back({{f + i + .5, j + .5, 0}, true, false});
    }
    auto eval = [&](int* kept_out) {
        std::vector<slab::Sq> obs; int kept = 0;
        for (auto& c : cand) if (c.keep) { obs.push_back(c.q); kept++; }
        // optional squares of the two lattices could overlap each other: reject configuration if so
        for (size_t i = 0; i < obs.size(); i++) (void)i;
        auto F = slab::dp(S, th, obs, slab::Window(), step, tol, 40);
        if (kept_out) *kept_out = kept;
        return kept + F.count;
    };
    // overlap check between kept candidates (UL vs LR optional squares)
    auto valid = [&]() {
        std::vector<int> idx; for (int i = 0; i < (int)cand.size(); i++) if (cand[i].keep && cand[i].boundary) idx.push_back(i);
        for (int i = 0; i < (int)idx.size(); i++) for (int j = 0; j < (int)cand.size(); j++) {
            if (!cand[j].keep || idx[i] == j) continue;
            if (std::fabs(cand[idx[i]].q.x - cand[j].q.x) < 1 - tol && std::fabs(cand[idx[i]].q.y - cand[j].q.y) < 1 - tol) return false;
        }
        return true;
    };
    std::mt19937_64 rng(5);
    int cur = eval(nullptr), best = cur;
    std::printf("start total %d\n", cur);
    std::vector<int> bidx; for (int i = 0; i < (int)cand.size(); i++) if (cand[i].boundary) bidx.push_back(i);
    for (int it = 0; it < iters; it++) {
        const int i = bidx[rng() % bidx.size()];
        cand[i].keep = !cand[i].keep;
        int tot = -1;
        if (valid()) tot = eval(nullptr);
        if (tot >= cur) { if (tot > cur) std::printf("it %d: %d\n", it, tot); cur = tot; best = std::max(best, cur); }
        else cand[i].keep = !cand[i].keep;
    }
    std::printf("best %d\n", best);
}
