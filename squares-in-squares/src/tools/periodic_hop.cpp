// Periodic strip cell with hopping search and optionally freed stair layers.
//   usage: periodic_hop r p m near f seconds [seed] [threads]
//   Fixed: UL squares (i+.5, f+j+.5) with j-i >= near, LR squares (f+i+.5, j+.5) with j-i <= -1-r-near.
//   Free : m squares per period p*(1,1): the `near` stair layers on each side (start on-lattice) + the rest random in the corridor.
#include <cstdio>
#include <cstdlib>
#include <random>
#include "../window.hpp"
int main(int argc, char** argv) {
    if (argc < 7) { std::fprintf(stderr, "usage: r p m near f seconds [seed] [threads]\n"); return 1; }
    const int r = atoi(argv[1]), p = atoi(argv[2]), m = atoi(argv[3]), near = atoi(argv[4]);
    const double f = atof(argv[5]), secs = atof(argv[6]);
    const uint64_t seed = argc > 7 ? atoll(argv[7]) : 1;
    const int threads = argc > 8 ? atoi(argv[8]) : 8;
    window::Env env; env.S = 1e9; env.walls = false; env.px = p; env.py = p;
    for (int i = -4 * p; i < 5 * p; i++)
        for (int d = -r - near - 6; d <= near + 6; d++) {
            const int j = i + d;
            if (d >= near) env.frozen.push_back({i + .5, f + j + .5, 0});
            if (d <= -1 - r - near) env.frozen.push_back({f + i + .5, j + .5, 0});
        }
    // frozen hash covers coordinates in [-1, S+1]; shift everything by +off so coordinates stay positive
    const double off = 4.0 * p + 8;
    for (auto& q : env.frozen) { q.x += off; q.y += off; }
    env.S = 10 * p + 40; env.index(); env.wrap_a0 = off; env.band_lo = -r - near - f - 1.2; env.band_hi = near + f + 1.2;
    std::vector<double> v;
    for (int L = 0; L < 2; L++)     // near layers on both sides start on-lattice
        for (int k = 0; k < near; k++)
            for (int i = 0; i < p; i++) {
                if (L == 0) { const int d = k, j = i + d; v.push_back(i + .5 + off); v.push_back(f + j + .5 + off); v.push_back(0); }
                else { const int d = -r - near + k, j = i + d; v.push_back(f + i + .5 + off); v.push_back(j + .5 + off); v.push_back(0); }
            }
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0, 1);
    // random corridor squares: along-strip position a in [0,p), across y-x in [-r-0.5, 0.5]
    while ((int)v.size() / 3 < m) {
        const double a = p * U(rng), dd = -r - 0.5 + (r + 1) * U(rng) * 0.9;
        v.push_back(a - dd / 2 + off); v.push_back(a + dd / 2 + off); v.push_back(U(rng) * M_PI / 2);
    }
    window::HopOptions o; o.seconds = secs; o.threads = threads; o.seed = seed; o.max_iter = 300;
    std::vector<double> out; double be = 0;
    const bool ok = window::hop(env, v, o, out, &be);
    std::printf("r=%d p=%d m=%d near=%d f=%.4f: %s (best E %.2e)\n", r, p, m, near, f, ok ? "FEASIBLE" : "no", be);
    if (ok) for (size_t i = 0; i < out.size() / 3; i++) std::printf("  %.12f %.12f %.10f\n", out[3 * i] - off, out[3 * i + 1] - off, std::fmod(out[3 * i + 2] * 180 / M_PI, 90.0));
}
