// usage: fill_motif motif.txt k f   : for each (a,s) tile the motif, keep the strip squares as the tilted core, MIS-refill uprights
#include <cstdio>
#include <cstdlib>
#include "../motif.hpp"
#include "../fill.hpp"
int main(int argc, char** argv) {
    motif::Motif m; motif::load(argv[1], m, 2, 4);
    const int k = atoi(argv[2]); const double f = atof(argv[3]), S = k + f;
    std::mt19937_64 rng(3);
    int best = 0;
    for (int a = -2; a <= 4; a++) for (int s = 0; s < 4; s++) {
        auto B = motif::build(m, S, a, s);
        std::vector<Square> T;
        for (int i = 0; i < B.n_strip; i++) T.push_back({B.v[3 * i], B.v[3 * i + 1], B.v[3 * i + 2]});
        auto R = fill::upright_fill(S, S, T, 1e-4, 100, rng);
        const int tot = (int)T.size() + (int)R.up.size();
        std::printf("a=%d s=%d strip %zu fill %zu total %d (lattice-body total %d)\n", a, s, T.size(), R.up.size(), tot, B.total());
        best = std::max(best, tot);
    }
    std::printf("best %d\n", best);
}
