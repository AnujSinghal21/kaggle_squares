// Write motif-tiled body packings (strictly feasible after a relax) as donor packings for the gadget DP.
//   usage: motif_donors motif.txt r p f out.csv [k_min k_max]
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "../motif.hpp"
#include "../squeeze.hpp"
int main(int argc, char** argv) {
    if (argc < 6) return 1;
    motif::Motif m;
    if (!motif::load(argv[1], m, atoi(argv[2]), atoi(argv[3]))) return 1;
    const double f = atof(argv[4]);
    const int k0 = argc > 6 ? atoi(argv[6]) : 12, k1 = argc > 7 ? atoi(argv[7]) : 31;
    std::ofstream out(argv[5]);
    out << "n,square_id,x,y,angle\n";
    out.precision(17);
    int written = 0;
    for (int k = k0; k <= k1; k++)
        for (int a = -2; a <= 3; a++)
            for (int s = 0; s < m.p; s++) {
                const double S = k + f;
                auto B = motif::build(m, S, a, s);
                std::vector<double> v = B.v;
                double Sx = S + 2 * squeeze::kMargin + 1e-9;
                bool ok = false;
                for (double slack = 1e-8; slack <= 4e-3 && !ok; slack *= 4) {
                    std::vector<double> t = v;
                    if (squeeze::relax(t, Sx + slack, 4000) <= squeeze::kTol) { v = t; Sx += slack; ok = true; }
                }
                if (!ok) continue;
                const int n = (int)v.size() / 3;
                for (int i = 0; i < n; i++) out << n << "," << i << "," << v[3 * i] << "," << v[3 * i + 1] << "," << v[3 * i + 2] << "\n";
                written++;
            }
    std::fprintf(stderr, "wrote %d donor packings\n", written);
}
