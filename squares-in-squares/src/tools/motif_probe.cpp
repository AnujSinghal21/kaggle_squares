// usage: motif_probe motif.txt r p f_lo f_hi [k_min k_max]   -> for each k the best count over corridor offset a and phase s
#include <cstdio>
#include <cstdlib>
#include "../motif.hpp"
int main(int argc, char** argv) {
    motif::Motif m;
    if (!motif::load(argv[1], m, atoi(argv[2]), atoi(argv[3]))) { std::fprintf(stderr, "cannot load motif\n"); return 1; }
    const double f = atof(argv[4]);
    const int k0 = argc > 5 ? atoi(argv[5]) : 14, k1 = argc > 6 ? atoi(argv[6]) : 31;
    std::printf("motif: %zu squares, period %d, f=%.4f\n", m.sq.size() / 3, m.p, f);
    for (int k = k0; k <= k1; k++) {
        int best = 0, ba = 0, bs = 0;
        for (int a = -4; a <= 8; a++) for (int s = 0; s < m.p; s++) {
            auto B = motif::build(m, k + f, a, s);
            if (B.total() > best) { best = B.total(); ba = a; bs = s; }
        }
        std::printf("k=%d S=%.4f best count %d (a=%d s=%d)\n", k, k + f, best, ba, bs);
    }
}
