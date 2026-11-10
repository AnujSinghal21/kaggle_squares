// usage: endcap_probe motif.txt r p k f a s corner(0=BL,1=TR) Rc extras seconds [seed]
// Frees the body squares within Rc of one container corner, adds `extras` random squares there, hops. Prints FEASIBLE or best E.
#include <cstdio>
#include <cstdlib>
#include <random>
#include "../motif.hpp"
#include "../squeeze.hpp"
#include "../window.hpp"
int main(int argc, char** argv) {
    if (argc < 12) { std::fprintf(stderr, "usage: motif r p k f a s corner Rc extras seconds [seed]\n"); return 1; }
    motif::Motif m;
    const int r = atoi(argv[2]), p = atoi(argv[3]), k = atoi(argv[4]);
    if (!motif::load(argv[1], m, r, p)) return 1;
    const double f = atof(argv[5]); const int a = atoi(argv[6]), s = atoi(argv[7]), corner = atoi(argv[8]);
    const double Rc = atof(argv[9]); const int extras = atoi(argv[10]); const double secs = atof(argv[11]);
    const uint64_t seed = argc > 12 ? atoll(argv[12]) : 1;
    const double S = k + f;
    auto B = motif::build(m, S, a, s);
    window::Env env; env.S = S + 2 * squeeze::kMargin;
    std::vector<double> fv;
    const double cx = corner ? S : 0, cy = corner ? S : 0;
    const bool slide = std::getenv("SLIDE") != nullptr;   // free the whole strip (not only the corner) and allow global slides
    for (int i = 0; i < B.total(); i++) {
        const double x = B.v[3 * i], y = B.v[3 * i + 1], t = B.v[3 * i + 2];
        if (std::hypot(x - cx, y - cy) < Rc || (slide && i < B.n_strip)) { fv.push_back(x); fv.push_back(y); fv.push_back(t); }
        else env.frozen.push_back({x, y, t});
    }
    env.index();
    const int nfree0 = (int)fv.size() / 3;
    std::mt19937_64 rng(seed); std::uniform_real_distribution<double> U(0, 1);
    for (int e = 0; e < extras; e++) {
        const double rr = Rc * std::sqrt(U(rng)), ang = U(rng) * M_PI / 2;
        const int cc = corner == 2 ? (int)(rng() % 2) : corner;
        double x = cc ? S - rr * std::cos(ang) : rr * std::cos(ang), y = cc ? S - rr * std::sin(ang) : rr * std::sin(ang);
        fv.push_back(std::min(std::max(x, 0.5), S - 0.5)); fv.push_back(std::min(std::max(y, 0.5), S - 0.5)); fv.push_back(U(rng) * M_PI / 2);
    }
    window::HopOptions o; o.seconds = secs; o.threads = 8; o.seed = seed; o.max_iter = 250;
    if (slide) { o.shift_dx = 1; o.shift_dy = 1; }
    std::vector<double> sol; double be = 0;
    const bool ok = window::hop(env, fv, o, sol, &be);
    std::printf("body %d, corner %d Rc %.1f free %d (+%d): %s (best E %.2e)\n", B.total(), corner, Rc, nfree0, extras, ok ? "FEASIBLE" : "no", be);
}
