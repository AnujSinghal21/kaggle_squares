// Approach: periodic strip motif (found by tools/periodic_hop) tiled along the diagonal corridor, end caps re-solved.
// RESULT (2026-10-27): found a period-4, 10-square diagonal strip motif at f~0.7975 (tools/periodic_hop) matching the leaders' k+0.796..0.80, but tiling it + end caps (hop) only reaches 299 at f~0.835 (worse than ours 0.8229); end caps are the unsolved part.
//   Idea:    the leaders' k + ~0.797 family is a period-4, 10-square tilted motif (r=2 corridor, density 2.5 per step) that
//            needs f >= 0.7975 versus 0.8229 for the sqrt7 lattice.  Body = motif tiling (frozen), the two corner regions
//            get the missing squares added at random and are basin-hopped (src/window.hpp) at the target side.
//   Run:     MOTIF_FILE=/path/m.txt MOTIF_F0=0.7985 MOTIF_TARGETS=299 MOTIF_SECONDS=60 make run SOLVER=motif
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <thread>

#include "../motif.hpp"
#include "../solve.hpp"
#include "../squeeze.hpp"
#include "../window.hpp"

namespace {
bool wanted(int n) {
    const char* e = std::getenv("MOTIF_TARGETS");
    if (!e) return false;
    std::stringstream ss(e);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        const size_t dash = tok.find('-');
        if (dash == std::string::npos) { if (std::atoi(tok.c_str()) == n) return true; }
        else if (n >= std::atoi(tok.substr(0, dash).c_str()) && n <= std::atoi(tok.substr(dash + 1).c_str())) return true;
    }
    return false;
}
double env_d(const char* k, double d) { const char* e = std::getenv(k); return e ? std::atof(e) : d; }
}  // namespace

Packing solve(int n, const Packing& base) {
    if (base.empty() || !wanted(n)) return base;
    const char* file = std::getenv("MOTIF_FILE");
    if (!file) { std::fprintf(stderr, "set MOTIF_FILE\n"); return base; }
    motif::Motif m;
    const int r = (int)env_d("MOTIF_R", 2), p = (int)env_d("MOTIF_P", 4);
    if (!motif::load(file, m, r, p)) { std::fprintf(stderr, "cannot load %s\n", file); return base; }
    const double f0 = env_d("MOTIF_F0", m.f), seconds = env_d("MOTIF_SECONDS", 60), Rc = env_d("MOTIF_RC", 5.0);
    const int maxmiss = (int)env_d("MOTIF_MAXMISS", 8);
    const double L0 = container_side(base);
    const int k = (int)std::floor(L0 - f0 - 1e-6);
    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    std::fprintf(stderr, "n=%d: motif end-cap search, k=%d, base side %.6f (f %.4f), motif f0 %.4f\n", n, k, L0, L0 - k, f0);
    Packing best = base;
    double best_L = L0;
    // candidate sides: from just below the base down to the motif's f0
    std::vector<double> fs;
    for (double f = L0 - k - 2e-4; f >= f0 - 1e-9; f -= std::max(0.002, (L0 - k - f0) / 8)) fs.push_back(f);
    if (fs.empty() || fs.back() > f0 + 1e-9) fs.push_back(f0);
    for (double f : fs) {
        const double S = k + f;
        struct Cand { int a, s, total; };
        std::vector<Cand> cands;
        for (int a = -2; a <= 4; a++)
            for (int s = 0; s < p; s++) { auto B = motif::build(m, S, a, s); if (B.total() >= n - maxmiss && B.total() <= n + 30) cands.push_back({a, s, B.total()}); }
        std::sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) { return x.total > y.total; });
        if ((int)cands.size() > (int)env_d("MOTIF_MAXCANDS", 4)) cands.resize((int)env_d("MOTIF_MAXCANDS", 4));
        std::fprintf(stderr, "  S=%.5f (f=%.4f): %zu candidate (a,s); best body count %d\n", S, f, cands.size(), cands.empty() ? 0 : cands[0].total);
        bool any = false;
        for (const Cand& c : cands) {
            auto B = motif::build(m, S, c.a, c.s);
            const int tot = B.total();
            const int need = n - tot;
            // split the body into end A (near (0,0)), end B (near (S,S)) and the frozen rest
            std::vector<double> endv[2];
            std::vector<window::Fixed> rest;
            for (int i = 0; i < tot; i++) {
                const double x = B.v[3 * i], y = B.v[3 * i + 1], t = B.v[3 * i + 2];
                const double dA = std::hypot(x, y), dB = std::hypot(S - x, S - y);
                if (dA < Rc) { endv[0].push_back(x); endv[0].push_back(y); endv[0].push_back(t); }
                else if (dB < Rc) { endv[1].push_back(x); endv[1].push_back(y); endv[1].push_back(t); }
                else rest.push_back({x, y, t});
            }
            std::vector<double> best_sol[2];
            int gain[2] = {0, 0};
            for (int e = 0; e < 2; e++) {
                // frozen = rest + the other end as it is
                window::Env env; env.S = S + 2 * squeeze::kMargin; env.frozen = rest;
                for (size_t q = 0; q < endv[1 - e].size(); q += 3) env.frozen.push_back({endv[1 - e][q], endv[1 - e][q + 1], endv[1 - e][q + 2]});
                env.index();
                best_sol[e] = endv[e];
                const int cap = std::max(need, 0) + (need <= 0 ? 0 : 0);
                for (int extra = 1; extra <= std::max(cap, 1); extra++) {
                    std::vector<double> v = endv[e];
                    std::mt19937_64 rng(31 * n + 7 * c.a + c.s + 1000 * e + extra);
                    std::uniform_real_distribution<double> U(0, 1);
                    for (int q = 0; q < extra; q++) {
                        const double rr = Rc * std::sqrt(U(rng)), ang = U(rng) * M_PI / 2;
                        const double x = e == 0 ? rr * std::cos(ang) : S - rr * std::cos(ang), y = e == 0 ? rr * std::sin(ang) : S - rr * std::sin(ang);
                        v.push_back(std::min(std::max(x, 0.5), S - 0.5)); v.push_back(std::min(std::max(y, 0.5), S - 0.5)); v.push_back(U(rng) * M_PI / 2);
                    }
                    window::HopOptions o;
                    o.seconds = seconds; o.threads = threads; o.seed = 5 + extra + 100 * e; o.max_iter = 250;
                    std::vector<double> sol; double be = 0;
                    const bool ok = window::hop(env, v, o, sol, &be);
                    std::fprintf(stderr, "    a=%d s=%d end %c +%d: %s (E %.1e)\n", c.a, c.s, e == 0 ? 'A' : 'B', extra, ok ? "FEASIBLE" : "no", be);
                    if (!ok) break;
                    gain[e] = extra; best_sol[e] = sol;
                }
            }
            std::fprintf(stderr, "    a=%d s=%d body %d need +%d: end gains A=%d B=%d\n", c.a, c.s, tot, need, gain[0], gain[1]);
            if (gain[0] + gain[1] < need) continue;
            // assemble: rest + both ends (use only as many extras as needed: surplus squares dropped from the ends' last additions)
            std::vector<double> allv;
            for (const window::Fixed& q : rest) { allv.push_back(q.x); allv.push_back(q.y); allv.push_back(q.t); }
            for (int e = 0; e < 2; e++) for (double q : best_sol[e]) allv.push_back(q);
            while ((int)allv.size() / 3 > n) { const size_t q = allv.size() / 3 - 1; allv.resize(3 * q); }   // drop surplus from the end of the list
            double Sa = S + 2 * squeeze::kMargin + 1e-9;
            bool good = false;
            for (double slack = 0; slack <= 2e-3 && !good; slack = slack == 0 ? 1e-8 : slack * 4) {
                std::vector<double> t2 = allv;
                if (squeeze::relax(t2, Sa + slack, 4000) <= squeeze::kTol) { allv = t2; Sa += slack; good = true; }
            }
            if (!good) { std::fprintf(stderr, "      could not make it strictly feasible\n"); continue; }
            squeeze::polish(allv, Sa);
            squeeze::State st; st.v = allv; st.S = Sa;
            Packing pk = squeeze::to_packing(st);
            const double Lp = container_side(pk);
            std::fprintf(stderr, "      packing at side %.9f\n", Lp);
            if (Lp < best_L) { best_L = Lp; best = pk; }
            any = true;
            break;
        }
        if (!any) break;   // could not realise at this side: lower sides are harder
    }
    std::fprintf(stderr, "n=%d: -> side %.9f (was %.9f)\n", n, best_L, L0);
    return best;
}
