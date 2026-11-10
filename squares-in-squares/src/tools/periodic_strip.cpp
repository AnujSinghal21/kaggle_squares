// Experiment: minimum grid offset f for a periodic diagonal strip.
//   UL grid squares (i+.5, f+j+.5) with j-i >= 0, LR grid squares (f+i+.5, j+.5) with j-i <= -1-r (both rigid),
//   m free squares per period p*(1,1).  Minimise f + mu*E over free squares and f.  Prints the best f found.
//   build: g++ -O3 -march=native -std=c++17 -pthread src/tools/periodic_strip.cpp -o build/periodic_strip
//   usage: periodic_strip r p m seconds f [threads]
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <random>
#include <thread>
#include <vector>
#include "../strip.hpp"
#include "../squeeze.hpp"

struct Prob { int r, p, m; };

// energy over variables z = [x,y,t]*m + f
static double E_of(const Prob& P, const std::vector<double>& z) {
    const int m = P.m;
    const double f = z[3 * m];
    double E = 0;
    const double kMargin = 0;
    // fixed stairs near the strip
    for (int q = 0; q < m; q++) {
        const double x = z[3 * q], y = z[3 * q + 1], t = z[3 * q + 2];
        for (int tt = -1; tt <= 1; tt++) {
            const double X = x + tt * P.p, Y = y + tt * P.p;
            // grid squares within reach
            const int i0 = (int)std::floor(X) - 2, i1 = (int)std::floor(X) + 2, j0 = (int)std::floor(Y - f) - 2, j1 = (int)std::floor(Y - f) + 2;
            for (int i = i0; i <= i1; i++)
                for (int j = j0; j <= j1; j++) {
                    if (j - i >= 0) {  // UL
                        const double d = strip::pen(i + .5, f + j + .5, 0, X, Y, t) + kMargin;
                        if (d > 0) E += d * d;
                    }
                }
            const int i2 = (int)std::floor(X - f) - 2, i3 = (int)std::floor(X - f) + 2, j2 = (int)std::floor(Y) - 2, j3 = (int)std::floor(Y) + 2;
            for (int i = i2; i <= i3; i++)
                for (int j = j2; j <= j3; j++) {
                    if (j - i <= -1 - P.r) {  // LR
                        const double d = strip::pen(f + i + .5, j + .5, 0, X, Y, t) + kMargin;
                        if (d > 0) E += d * d;
                    }
                }
            for (int q2 = 0; q2 < m; q2++) {
                if (tt == 0 && q2 <= q) continue;
                const double d = strip::pen(z[3 * q2], z[3 * q2 + 1], z[3 * q2 + 2], X, Y, t);
                if (d > 0) E += d * d;
            }
            if (tt != 0) {  // image of q vs q itself handled by q2 loop only for q2>q at tt==0; add self image
                const double d = strip::pen(x, y, t, X, Y, t);
                if (d > 0) E += d * d;
            }
        }
    }
    return E;
}

static double grad(const Prob& P, std::vector<double>& z, double mu, std::vector<double>& g) {
    const int N = (int)z.size();
    g.assign(N, 0);
    const double h = 1e-6;
    for (int k = 0; k < N; k++) {
        const double s = z[k];
        z[k] = s + h; const double a = E_of(P, z);
        z[k] = s - h; const double b = E_of(P, z);
        z[k] = s;
        g[k] = mu * (a - b) / (2 * h);
    }
    g[N - 1] = 0;  // f is held fixed
    return mu * E_of(P, z);
}

int main(int argc, char** argv) {
    if (argc < 6) { std::fprintf(stderr, "usage: r p m seconds f [threads]\n"); return 1; }
    Prob P{std::atoi(argv[1]), std::atoi(argv[2]), std::atoi(argv[3])};
    const double secs = std::atof(argv[4]);
    const double fixed_f = std::atof(argv[5]);
    const int threads = argc > 6 ? std::atoi(argv[6]) : 8;
    std::mutex mu;
    double best = 1e9;
    std::vector<double> bz;
    std::atomic<bool> found{false};
    const auto t_end = std::chrono::steady_clock::now() + std::chrono::duration<double>(secs);
    long runs = 0;
    auto worker = [&](int tid) {
        std::mt19937_64 rng(1234 + tid);
        std::uniform_real_distribution<double> U(0, 1);
        while (std::chrono::steady_clock::now() < t_end && !found) {
            std::vector<double> z(3 * P.m + 1);
            const double f0 = fixed_f;
            for (int q = 0; q < P.m; q++) {
                // along-strip position in [0,p), across offset: y-x centre around -(r+1)/2
                const double a = P.p * U(rng), cx = -(P.r + 1) / 2.0 + (U(rng) - 0.5) * (P.r + 1);   // y - x
                z[3 * q] = a - cx / 2; z[3 * q + 1] = a + cx / 2; z[3 * q + 2] = U(rng) * M_PI / 2;
            }
            z[3 * P.m] = f0;
            for (double mu : {1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8}) {
                auto F = [&](const std::vector<double>& x, std::vector<double>& g) { std::vector<double> y = x; return grad(P, y, mu, g); };
                squeeze::lbfgs(z, F, 400, [](double) { return false; });
                for (int q = 0; q < P.m; q++) {  // keep every free square inside the reference period along the strip
                    const double a = 0.5 * (z[3 * q] + z[3 * q + 1]), sh = std::floor(a / P.p) * P.p;
                    z[3 * q] -= sh; z[3 * q + 1] -= sh;
                }
            }
            const double E = E_of(P, z);
            std::lock_guard<std::mutex> lk(mu);
            runs++;
            if (E < 1e-14 && (bz.empty() || z[3 * P.m] < best)) {
                best = z[3 * P.m]; bz = z; found = true;
                std::fprintf(stderr, "[t%d] f=%.6f E=%.2e\n", tid, best, E);
            }
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    std::printf("r=%d p=%d m=%d runs=%ld best f=%.6f\n", P.r, P.p, P.m, runs, best);
    for (int q = 0; q < P.m && !bz.empty(); q++) std::printf("  sq %d: x=%.5f y=%.5f deg=%.3f\n", q, bz[3 * q], bz[3 * q + 1], bz[3 * q + 2] * 180 / M_PI);
}
