// Approach: column-DP strip composition (src/fault.hpp) + L-BFGS relax + squeeze polish.
// RESULT (2026-10-27): dead end for improving, useful as a model. Column DP with rigid UL/LR stairs reproduces n=299 exactly (256+43) but its tolerance-based counts are optimistic; real relaxable side is ~the existing 17.8241. Leaders' k+0.796 structure is NOT in this model class.
//   Idea:    best packings for n > ~250 are two corner-anchored grids separated by a diagonal corridor filled with tilted
//            columns (stacks of h unit squares) that slide independently to nest in the stairs.  For a given container
//            side S, tilt theta and corridor cut (a, r) the most squares such a strip holds is a 1-D dynamic program.
//            We bisect S for the smallest side whose best count reaches n, build the packing, relax away the DP
//            tolerance and polish.  Related: Deotte's fault-line gadget DP (Kaggle discussion), Friedman's strip families.
//   Run:     FAULT_TARGETS=300 FAULT_SECONDS=120 make run SOLVER=fault
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <sstream>
#include <thread>

#include "../fault.hpp"
#include "../solve.hpp"
#include "../squeeze.hpp"

namespace {
bool wanted(int n) {
    const char* e = std::getenv("FAULT_TARGETS");
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

struct Hit { fault::Params P; int total; };

// All (theta, a, r) whose DP count at side S is >= n, best first.
std::vector<Hit> hits_at(double S, int n, int threads, double th_lo, double th_hi, double th_step, double tol) {
    struct Job { double th; int a, r; };
    std::vector<Job> jobs;
    const int k = (int)std::floor(S);
    (void)k;
    // corridor position: the strip runs near the main diagonal; a (UL cut) in a small window around 0
    for (double th = th_lo; th <= th_hi + 1e-9; th += th_step)
        for (int a = 0; a <= 3; a++)
            for (int r = 1; r <= 3; r++) jobs.push_back({th, a, r});
    std::atomic<size_t> next{0};
    std::mutex mu;
    std::vector<Hit> out;
    auto worker = [&]() {
        for (;;) {
            const size_t i = next++;
            if (i >= jobs.size()) return;
            fault::Params P;
            P.S = S; P.theta = jobs[i].th * M_PI / 180; P.a = jobs[i].a; P.b = jobs[i].a - 1 - jobs[i].r; P.tol = tol;
            const fault::Result R = fault::solve(P, false);
            if (R.total() >= n) { std::lock_guard<std::mutex> lk(mu); out.push_back({P, R.total()}); }
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    std::sort(out.begin(), out.end(), [&](const Hit& x, const Hit& y) {
        if (x.total != y.total) return x.total < y.total;  // fewest surplus squares first (least trimming)
        return x.P.theta < y.P.theta;
    });
    return out;
}

// Build the full state (strip squares first) for a hit and trim down to n squares.
bool build(const Hit& h, int n, std::vector<double>& v, std::mt19937_64& rng) {
    const fault::Result R = fault::solve(h.P, true);
    const double S = h.P.S, ct = std::cos(h.P.theta), st = std::sin(h.P.theta);
    const int k = (int)std::floor(S);
    const double f = S - k;
    v.clear();
    int ns = 0;
    for (const fault::Column& c : R.cols)
        for (int q = 0; q < c.h; q++) {
            const double w = c.w + (q - (c.h - 1) * 0.5), u = c.u;
            v.push_back(u * ct - w * st); v.push_back(u * st + w * ct); v.push_back(h.P.theta);
            ns++;
        }
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++) {
            if (j - i >= h.P.a) { v.push_back(i + .5); v.push_back(f + j + .5); v.push_back(0); }
            else if (j - i <= h.P.b) { v.push_back(f + i + .5); v.push_back(j + .5); v.push_back(0); }
        }
    return strip::trim_to(v, ns, n, rng);
}
}  // namespace

Packing solve(int n, const Packing& base) {
    if (base.empty() || !wanted(n)) return base;
    const double seconds = env_d("FAULT_SECONDS", 120), tol = env_d("FAULT_TOL", 0.015);
    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    const double L0 = container_side(base);
    const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    const double th_lo = env_d("FAULT_TH_LO", 14), th_hi = env_d("FAULT_TH_HI", 42), th_step = env_d("FAULT_TH_STEP", 0.5);
    std::fprintf(stderr, "n=%d: fault search below %.9f (budget %.0fs)\n", n, L0, seconds);
    const int k = (int)std::floor(L0 - 0.3);  // the grid side; f in (0.3, ...]
    // bisect S in [k + 0.3, L0 - 1e-4] on "some config reaches n"
    double lo = k + 0.3, hi = L0 - 1e-4;
    if (hits_at(hi, n, threads, th_lo, th_hi, th_step, tol).empty()) {
        std::fprintf(stderr, "n=%d: DP finds nothing below the base side\n", n);
        return base;
    }
    if (std::getenv("FAULT_SSTART")) hi = std::atof(std::getenv("FAULT_SSTART"));
    else for (int it = 0; it < 14 && std::chrono::steady_clock::now() < t_end; it++) {
        const double mid = 0.5 * (lo + hi);
        if (!hits_at(mid, n, threads, th_lo, th_hi, th_step, tol).empty()) hi = mid; else lo = mid;
    }
    std::fprintf(stderr, "  DP count reaches %d at S ~ %.5f (base %.5f)\n", n, hi, L0);
    // realise: try sides from hi upward until a hit relaxes to zero energy
    Packing best = base;
    double best_S = L0;
    std::mt19937_64 rng(4242 + n);
    for (double S = hi; S < L0 - 1e-6 && std::chrono::steady_clock::now() < t_end; S += 0.0005) {
        std::vector<Hit> hs = hits_at(S, n, threads, th_lo, th_hi, th_step, tol);
        int tried = 0;
        for (const Hit& h : hs) {
            if (++tried > 12) break;
            std::vector<double> v;
            if (!build(h, n, v, rng)) continue;
            const double Sm = S + 2 * squeeze::kMargin + 1e-9;
            const double E = squeeze::relax(v, Sm, 4000);
            if (std::getenv("FAULT_VERBOSE")) std::fprintf(stderr, "  S=%.5f theta %.2f a=%d r=%d total %d: relax E=%.3e\n", S, h.P.theta * 180 / M_PI, h.P.a, h.P.a - 1 - h.P.b, h.total, E);
            if (E > squeeze::kTol) continue;
            double Sp = Sm;
            squeeze::polish(v, Sp);
            squeeze::State st; st.v = v; st.S = Sp;
            Packing p = squeeze::to_packing(st);
            const double Lp = container_side(p);
            std::fprintf(stderr, "  realised at S=%.6f: theta %.2f a=%d b=%d total %d -> L=%.9f\n", S, h.P.theta * 180 / M_PI, h.P.a, h.P.b, h.total, Lp);
            if (Lp < best_S) { best_S = Lp; best = p; }
            if (Lp < L0) goto done;
        }
    }
done:
    if (best_S < L0) {
        // squeeze the realised packing down further
        squeeze::Options o;
        o.seconds = std::max(5.0, std::chrono::duration<double>(t_end - std::chrono::steady_clock::now()).count());
        o.threads = threads; o.seed = 99 + n;
        Packing q = squeeze::run(best, o);
        if (container_side(q) < best_S) best = q;
        std::fprintf(stderr, "n=%d: -> side %.9f (was %.9f)\n", n, container_side(best), L0);
    }
    return best;
}
