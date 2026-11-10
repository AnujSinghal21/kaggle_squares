// Densest periodic packing of unit squares in a vertical strip of width W (flat walls), period p along y, m squares per period.
//   usage: strip_periodic W p m seconds [seed]   -> FEASIBLE / no
#include <cstdio>
#include <cstdlib>
#include <random>
#include "../window.hpp"
int main(int argc, char** argv) {
    const double W = atof(argv[1]); const double p = atof(argv[2]); const int m = atoi(argv[3]); const double secs = atof(argv[4]);
    const uint64_t seed = argc > 5 ? atoll(argv[5]) : 1;
    window::Env env; env.S = W + 2 * squeeze::kMargin; env.walls = true; env.wall_y = false; env.px = 0; env.py = p; env.wrap_a0 = 0; env.index();
    // wrap along y: reuse wrap logic via a custom approach: positions are wrapped manually by restricting y to [0,p) after relax (done in hop through wrap_a0 using x+y: not valid here), so wrap in a thin loop below
    std::mt19937_64 rng(seed); std::uniform_real_distribution<double> U(0, 1);
    std::vector<double> v;
    for (int i = 0; i < m; i++) { v.push_back(0.5 + (W - 1) * U(rng)); v.push_back(p * U(rng)); v.push_back(U(rng) * M_PI / 2); }
    window::HopOptions o; o.seconds = secs; o.threads = 8; o.seed = seed; o.max_iter = 300;
    std::vector<double> out; double be = 0;
    const bool ok = window::hop(env, v, o, out, &be);
    std::printf("W=%.4f p=%d/1000 m=%d (density %.4f/row, extras %.4f/row over %d): %s (E %.2e)\n", W, (int)(p*1000), m, (double)m / p, (double)m / p - std::floor(W), (int)std::floor(W), ok ? "FEASIBLE" : "no", be);
}
