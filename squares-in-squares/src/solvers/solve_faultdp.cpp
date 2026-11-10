// Approach: fault-line gadget DP (src/gadget.hpp): n = K^2 + D, two lattices shifted by F, a chain of gadgets cropped from donor
// RESULT (2026-10-27): reproduces the public fault-line DP (our packings are its optimum); gains come only with improved gadgets (solve_gquench.cpp). Compose must dilate by 1+2e-11 or exact lattice contacts fail validation.
//   packings along the seam, DP minimising the largest overhang f.  Own C++ implementation of the composition described in
//   Chris Deotte's "multi-gadget fault DP" notebook and in alex chilton's starter kit (CC0 code, read as documentation).
//   Donors: every book + the newest submission + FDP_DONORS (comma-separated extra csv files, e.g. couzo_big 1000 < n <= 3500).
//   Run:    FDP_DONORS=books/raw/starter_kit_full/couzo_big.csv make run SOLVER=faultdp ARGS="325 1000"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <thread>

#include "../gadget.hpp"
#include "../solve.hpp"

namespace {
namespace fs = std::filesystem;
using gadget::Donor;

struct State {
    bool ready = false;
    std::vector<Donor> donors;
    std::vector<gadget::Gadget> G;
    gadget::Dp dp;
};
State& S() { static State s; return s; }

// header-driven csv loader: columns n, square_id (optional), x, y, angle
void load_donors(const std::string& path, std::map<std::pair<int, long long>, Donor>& out) {
    std::ifstream in(path);
    if (!in) { std::fprintf(stderr, "faultdp: cannot open %s\n", path.c_str()); return; }
    std::string line;
    std::getline(in, line);
    std::vector<std::string> cols;
    { std::stringstream ss(line); std::string c; while (std::getline(ss, c, ',')) cols.push_back(c); }
    int cn = -1, cx = -1, cy = -1, ca = -1;
    for (int i = 0; i < (int)cols.size(); i++) { if (cols[i] == "n") cn = i; if (cols[i] == "x") cx = i; if (cols[i] == "y") cy = i; if (cols[i] == "angle") ca = i; }
    if (cn < 0 || cx < 0 || cy < 0 || ca < 0) { std::fprintf(stderr, "faultdp: unexpected header in %s\n", path.c_str()); return; }
    std::map<int, Donor> cur;
    std::vector<std::string> f(cols.size());
    while (std::getline(in, line)) {
        size_t pos = 0; int k = 0;
        while (k < (int)cols.size()) { size_t q = line.find(',', pos); f[k++] = line.substr(pos, q == std::string::npos ? std::string::npos : q - pos); if (q == std::string::npos) break; pos = q + 1; }
        if (k < (int)cols.size()) continue;
        const int n = std::atoi(f[cn].c_str());
        Donor& d = cur[n];
        d.n = n; d.x.push_back(std::atof(f[cx].c_str())); d.y.push_back(std::atof(f[cy].c_str())); d.t.push_back(std::atof(f[ca].c_str()));
    }
    for (auto& kv : cur) {
        Donor& d = kv.second;
        if ((int)d.x.size() != d.n) continue;
        Packing p;
        for (int i = 0; i < d.n; i++) p.push_back({d.x[i], d.y[i], d.t[i]});
        d.L = container_side(p);
        if (d.L >= std::ceil(std::sqrt((double)d.n)) - 1e-9) continue;     // plain grids carry no gadget
        out.emplace(std::make_pair(d.n, (long long)std::llround(d.L * 1e9)), std::move(d));
    }
}


// external gadget file: "G w h dm f nsq" followed by nsq lines "x y angle" (any frame); all 8 dihedral variants are added
std::vector<Donor>& ext_donors() { static std::vector<Donor> v; return v; }

void load_gadget_file(const std::string& path, std::vector<gadget::Gadget>& out) {
    std::ifstream in(path);
    if (!in) { std::fprintf(stderr, "faultdp: cannot open gadget file %s\n", path.c_str()); return; }
    std::string line;
    int added = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] != 'G') continue;
        int w, h, dm, nsq; double fdecl;
        if (std::sscanf(line.c_str(), "G %d %d %d %lf %d", &w, &h, &dm, &fdecl, &nsq) != 5) continue;
        Donor d; d.n = nsq;
        for (int i = 0; i < nsq; i++) { std::getline(in, line); double x, y, t; std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &t); d.x.push_back(x); d.y.push_back(y); d.t.push_back(t); }
        for (int sym = 0; sym < 8; sym++) {
            std::vector<double> xlo, xhi, ylo, yhi, mx, my, mt;
            gadget::map_donor(d, sym, xlo, xhi, ylo, yhi, &mx, &my, &mt);
            const int ww = (sym & 4) ? h : w, hh = (sym & 4) ? w : h;
            double f = 0, ybr = 9, xtl = 9;
            for (int i = 0; i < nsq; i++) {
                f = std::max(f, std::max(xhi[i] - ww, yhi[i] - hh));
                if (xhi[i] > ww + 1e-9) ybr = std::min(ybr, ylo[i]);
                if (yhi[i] > hh + 1e-9) xtl = std::min(xtl, xlo[i]);
            }
            if (f > gadget::kFMax) continue;
            gadget::Gadget g{ww, hh, dm, f, std::min(ybr, 9.0), std::min(xtl, 9.0), -1, sym, gadget::bin_of(std::min(ybr, 9.0)), gadget::bin_of(std::min(xtl, 9.0))};
            for (int i = 0; i < nsq; i++) g.sq.push_back({mx[i], my[i], mt[i]});
            out.push_back(std::move(g));
            added++;
        }
        ext_donors().push_back(d);
        ext_donors().back().L = std::max(w, h) + fdecl;
    }
    std::fprintf(stderr, "faultdp: %d external gadgets (8 variants each) from %s\n", added, path.c_str());
}

void init() {
    State& s = S();
    if (s.ready) return;
    const auto t0 = std::chrono::steady_clock::now();
    auto el = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    std::vector<std::string> files;
    for (const auto& e : fs::directory_iterator("books")) if (e.path().extension() == ".csv") files.push_back(e.path().string());
    std::string newest; fs::file_time_type nt;
    for (const auto& e : fs::directory_iterator("submissions")) if (e.path().extension() == ".csv" && (newest.empty() || fs::last_write_time(e) > nt)) { newest = e.path().string(); nt = fs::last_write_time(e); }
    if (!newest.empty()) files.push_back(newest);
    if (const char* ex = std::getenv("FDP_DONORS")) { std::stringstream ss(ex); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) files.push_back(t); }
    std::map<std::pair<int, long long>, Donor> uniq;
    for (const std::string& f : files) { load_donors(f, uniq); std::fprintf(stderr, "faultdp: loaded %s -> %zu distinct donors [%.0fs]\n", f.c_str(), uniq.size(), el()); }
    for (auto& kv : uniq) { Donor d = std::move(kv.second); s.donors.push_back(std::move(d)); }
    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    s.G = gadget::extract(s.donors, threads);
    std::fprintf(stderr, "faultdp: %zu donors -> %zu gadgets after dominance [%.0fs]\n", s.donors.size(), s.G.size(), el());
    if (const char* gf = std::getenv("FDP_GADGETS")) { std::stringstream ss(gf); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) load_gadget_file(t, s.G); }
    if (const char* wi = std::getenv("FDP_WHATIF")) {   // "w,h,dm,f;...": pretend those keys (and transposes) have overhang f
        std::stringstream ks(wi); std::string tok;
        while (std::getline(ks, tok, ';')) {
            int w, h, dm; double f; if (std::sscanf(tok.c_str(), "%d,%d,%d,%lf", &w, &h, &dm, &f) != 4) continue;
            int cnt = 0;
            for (auto& g : s.G) if (g.dm == dm && ((g.w == w && g.h == h) || (g.w == h && g.h == w)) && g.f > f) { g.f = f; cnt++; }
            if (!cnt) { gadget::Gadget g{w, h, dm, f, 0.9, 0.9, -1, 0, gadget::bin_of(0.9), gadget::bin_of(0.9)}; s.G.push_back(g); g.w = h; g.h = w; s.G.push_back(g); }
            std::fprintf(stderr, "whatif %s: %d gadgets lowered\n", tok.c_str(), cnt);
        }
    }
    if (!ext_donors().empty()) {
        std::vector<gadget::Gadget> crops = gadget::extract(ext_donors(), threads);
        const int off = (int)s.donors.size();
        for (auto& d : ext_donors()) s.donors.push_back(std::move(d));
        for (auto& g : crops) { g.donor += off; s.G.push_back(std::move(g)); }
        std::fprintf(stderr, "faultdp: crop closure of external gadgets added %zu gadgets\n", crops.size());
    }
    s.dp = gadget::run_dp(s.G, threads);
    std::fprintf(stderr, "faultdp: DP done [%.0fs]\n", el());
    s.ready = true;
}
}  // namespace

Packing solve(int n, const Packing& base) {
    const int K = (int)std::floor(std::sqrt((double)n) + 1e-9), D = n - K * K;
    if (D == 0 || K > gadget::kWMax || D > gadget::kDMax || base.empty()) return base;
    init();
    State& s = S();
    std::vector<gadget::Link> links;
    const double f = gadget::best_chain(s.dp, s.G, K, D, links);
    const double L0 = container_side(base);
    if (std::getenv("FDP_REPORT")) std::fprintf(stderr, "REPORT n=%d K=%d D=%d dp=%.6f ours=%.6f\n", n, K, D, f < 0 ? -1.0 : K + f, L0);
    if (std::getenv("FDP_CHAIN") && f >= 0) {
        std::fprintf(stderr, "CHAIN n=%d K=%d D=%d f=%.6f ours=%.6f:", n, K, D, f, L0 - K);
        for (const auto& l : links) { const auto& g = s.G[l.g]; std::fprintf(stderr, " [%dx%d+%d f=%.4f m%d donor n=%d L=%.4f sym%d tips %.2f/%.2f]", g.w, g.h, g.dm, g.f, l.mode, g.donor >= 0 ? s.donors[g.donor].n : -1, g.donor >= 0 ? s.donors[g.donor].L : 0.0, g.sym, g.ybr, g.xtl); }
        std::fprintf(stderr, "\n");
    }
    if (f < 0 || K + f + 2e-9 >= L0 - 1e-7) return base;
    Packing p = gadget::compose(K, D, f, links, s.G, s.donors);
    if (p.empty()) { std::fprintf(stderr, "n=%d: DP says %.9f < %.9f but compose failed\n", n, K + f, L0); return base; }
    std::fprintf(stderr, "n=%d: DP side %.9f vs ours %.9f (%zu gadgets)\n", n, container_side(p), L0, links.size());
    return p;
}
