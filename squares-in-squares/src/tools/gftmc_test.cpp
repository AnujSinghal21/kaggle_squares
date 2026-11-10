// usage: gftmc_test gadgets.txt index seconds [seed]   (file format: "G w h dm f nsq" + squares) : run FTMC on gadget #index
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "../gftmc.hpp"
int main(int argc, char** argv) {
    std::ifstream in(argv[1]); std::string line; int idx = atoi(argv[2]), cur = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] != 'G') continue;
        int w, h, dm, nsq; double f;
        std::sscanf(line.c_str(), "G %d %d %d %lf %d", &w, &h, &dm, &f, &nsq);
        std::vector<Square> sq;
        for (int i = 0; i < nsq; i++) { std::getline(in, line); double x, y, t; std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &t); sq.push_back({x, y, t}); }
        if (cur++ != idx) continue;
        gadget::Gadget g{w, h, dm, f, 9, 9, -1, 0, 0, 0};
        // normalise to origin and recompute f
        gadget::Donor d; d.n = nsq; for (auto& q : sq) { d.x.push_back(q.x); d.y.push_back(q.y); d.t.push_back(q.angle); }
        std::vector<double> xlo, xhi, ylo, yhi, mx, my, mt; gadget::map_donor(d, 0, xlo, xhi, ylo, yhi, &mx, &my, &mt);
        double f2 = 0; for (int i = 0; i < nsq; i++) { f2 = std::max(f2, std::max(xhi[i] - w, yhi[i] - h)); g.sq.push_back({mx[i], my[i], mt[i]}); }
        g.f = f2;
        std::printf("gadget %dx%d+%d  f=%.6f (declared %.6f), %d squares\n", w, h, dm, f2, f, nsq);
        gadget::Donor none; gadget::Gadget out;
        std::vector<gftmc::CoreSeed> alts;
        if (argc > 5) {   // neighbour file: gadgets (after/before lists) whose key is adjacent to this one
            std::ifstream in2(argv[5]); std::string l2;
            while (std::getline(in2, l2)) {
                if (l2.empty() || l2[0] != 'G') continue;
                int w2, h2, dm2, n2; double f2d; std::sscanf(l2.c_str(), "G %d %d %d %lf %d", &w2, &h2, &dm2, &f2d, &n2);
                std::vector<Square> s2; for (int i = 0; i < n2; i++) { std::getline(in2, l2); double x, y, t; std::sscanf(l2.c_str(), "%lf %lf %lf", &x, &y, &t); s2.push_back({x, y, t}); }
                if (std::abs(w2 - w) + std::abs(h2 - h) + std::abs(dm2 - dm) > 2 || (w2 == w && h2 == h && dm2 == dm)) continue;
                gadget::Donor d2; d2.n = n2; for (auto& q : s2) { d2.x.push_back(q.x); d2.y.push_back(q.y); d2.t.push_back(q.angle); }
                std::vector<double> a1, a2, a3, a4, m1, m2, m3; gadget::map_donor(d2, 0, a1, a2, a3, a4, &m1, &m2, &m3);
                gftmc::CoreSeed cs; double fa = 0; for (int i = 0; i < n2; i++) { fa = std::max(fa, std::max(a2[i] - w2, a4[i] - h2)); Square q{m1[i], m2[i], m3[i]}; if (gftmc::is_tilted(q)) cs.T.push_back(q); }
                cs.W = w2 + fa; cs.H = h2 + fa; if (!cs.T.empty()) alts.push_back(cs);
            }
            std::printf("%zu neighbour cores\n", alts.size());
        }
        const bool ok = gftmc::ftmc_gadget(g, none, atof(argv[3]), argc > 4 ? atoll(argv[4]) : 1, out, &alts);
        std::printf("%s f -> %.6f\n", ok ? "IMPROVED" : "no improvement", ok ? out.f : f2);
        return 0;
    }
}
