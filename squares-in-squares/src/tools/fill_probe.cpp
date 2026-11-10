// usage: fill_probe file.csv n [S_override] [iters]  : keep the tilted squares of packing n, refill uprights, report count vs n
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "../book.hpp"
#include "../fill.hpp"
#include "../recompose.hpp"
int main(int argc, char** argv) {
    const int n = atoi(argv[2]);
    Packing p = book_packing(argv[1], n);
    auto sp = recompose::split(p);
    const double S = argc > 3 && atof(argv[3]) > 0 ? atof(argv[3]) : sp.L;
    const int iters = argc > 4 ? atoi(argv[4]) : 200;
    std::vector<Square> T;
    for (auto& t : sp.tilt) T.push_back({t.x, t.y, t.t});
    std::mt19937_64 rng(1);
    auto t0 = std::chrono::steady_clock::now();
    auto R = fill::upright_fill(S, S, T, 1e-9, iters, rng);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("n=%d S=%.6f tilted %zu, candidates %d, uprights found %zu (original %zu), total %zu vs n %d  [%.2fs]\n", n, S, T.size(), R.candidates, R.up.size(), sp.axis.size(), T.size() + R.up.size(), n, sec);
}
