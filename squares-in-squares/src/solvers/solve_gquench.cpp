// Approach (gadget quench): fault-line gadget DP + squeeze of the bottleneck gadgets (src/gadget.hpp) -- base of faultdp: n = K^2 + D, two lattices shifted by F, a chain of gadgets cropped from donor
// RESULT (2026-10-27): the productive approach so far. Plain-squeeze sweeps over high-overhang gadget keys (GQ_SWEEP, tau 10 s) + crop closure + DP compose took the score 22.2987572 -> 22.2924 (-0.0063); saturated at ~1e-5/cycle. FTMC modes (GQ_FTMC) improve some gadgets (3x14+4 .9726->.8156, 6x20+3 .7467->.6925) but not the rigid chain bottlenecks. Run one heavy process at a time (OOM kills).
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
#include <set>
#include <array>
#include <unordered_map>
#include <atomic>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

#include "../gadget.hpp"
#include "../solve.hpp"
#include "../squeeze.hpp"
#include "../gftmc.hpp"
#include "../window.hpp"

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


// All 8 dihedral variants of an explicit-square gadget (squares in any frame; frame is renormalised per variant).
std::vector<gadget::Gadget> variants_of(int w, int h, int dm, const std::vector<Square>& sq, int tag_donor = -1) {
    std::vector<gadget::Gadget> out;
    Donor d; d.n = (int)sq.size();
    for (const Square& q : sq) { d.x.push_back(q.x); d.y.push_back(q.y); d.t.push_back(q.angle); }
    for (int sym = 0; sym < 8; sym++) {
        std::vector<double> xlo, xhi, ylo, yhi, mx, my, mt;
        gadget::map_donor(d, sym, xlo, xhi, ylo, yhi, &mx, &my, &mt);
        const int ww = (sym & 4) ? h : w, hh = (sym & 4) ? w : h;
        double f = 0, ybr = 9, xtl = 9;
        for (int i = 0; i < d.n; i++) {
            f = std::max(f, std::max(xhi[i] - ww, yhi[i] - hh));
            if (xhi[i] > ww + 1e-9) ybr = std::min(ybr, ylo[i]);
            if (yhi[i] > hh + 1e-9) xtl = std::min(xtl, xlo[i]);
        }
        if (f > gadget::kFMax) continue;
        gadget::Gadget g{ww, hh, dm, f, std::min(ybr, 9.0), std::min(xtl, 9.0), tag_donor, sym, gadget::bin_of(std::min(ybr, 9.0)), gadget::bin_of(std::min(xtl, 9.0))};
        for (int i = 0; i < d.n; i++) g.sq.push_back({mx[i], my[i], mt[i]});
        out.push_back(std::move(g));
    }
    return out;
}

void save_gadgets(const std::vector<gadget::Gadget>& v, const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return;
    for (const auto& g : v) {
        std::fprintf(f, "G %d %d %d %.17g %zu\n", g.w, g.h, g.dm, g.f, g.sq.size());
        for (const Square& q : g.sq) std::fprintf(f, "%.17g %.17g %.17g\n", q.x, q.y, q.angle);
    }
    std::fclose(f);
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
        std::vector<Square> sq;
        for (int i = 0; i < nsq; i++) { std::getline(in, line); double x, y, t; std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &t); sq.push_back({x, y, t}); }
        for (auto& g : variants_of(w, h, dm, sq)) { out.push_back(std::move(g)); added++; }
        { Donor d; d.n = (int)sq.size(); for (const Square& q : sq) { d.x.push_back(q.x); d.y.push_back(q.y); d.t.push_back(q.angle); }
          d.L = std::max(w, h) + fdecl; ext_donors().push_back(std::move(d)); }
    }
    std::fprintf(stderr, "faultdp: %d external gadget variants from %s\n", added, path.c_str());
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
    if (const char* gf = std::getenv("FDP_GADGETS")) { std::stringstream ss(gf); std::string tt; while (std::getline(ss, tt, ',')) if (!tt.empty()) load_gadget_file(tt, s.G); }
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


namespace {
double env_d(const char* k, double d) { const char* e = std::getenv(k); return e ? std::atof(e) : d; }
std::map<int, double>& leaders() {
    static std::map<int, double> m;
    static bool loaded = false;
    if (!loaded) {
        loaded = true;
        const char* path = std::getenv("GQ_LEADER");
        if (path) { std::FILE* f = std::fopen(path, "r"); if (f) { char line[512]; while (std::fgets(line, sizeof line, f)) { int k; double L; if (std::sscanf(line, "%d,%lf", &k, &L) == 2) m[k] = L; } std::fclose(f); } }
    }
    return m;
}

using gftmc::tips_of;
// Squeeze one gadget (single thread). Returns true and fills `out` if its overhang f dropped.
bool quench(const gadget::Gadget& g, const Donor& d, double seconds, uint64_t seed, gadget::Gadget& out) {
    std::vector<Square> sq;
    double f0 = g.f;
    if (!g.sq.empty()) sq = g.sq; else gadget::crop(d, g.sym, g.w, g.h, g.dm, sq, &f0);
    squeeze::State st;
    for (const Square& q : sq) { st.v.push_back(q.x); st.v.push_back(q.y); st.v.push_back(q.angle); }
    st.yoff = g.h - g.w;
    bool feasible = false;
    for (double slack = 1e-9; slack <= 2e-5 && !feasible; slack *= 10) {
        squeeze::State t = st;
        t.S = g.w + f0 + 2 * squeeze::kMargin + slack;
        if (squeeze::relax(t.v, t.S, 4000, t.yoff) <= squeeze::kTol) { st = t; feasible = true; }
    }
    if (!feasible) return false;
    squeeze::Options o; o.seconds = seconds; o.threads = 1; o.seed = seed; o.max_iter = 300;
    squeeze::State r = squeeze::run_state(st, o);
    const double f1 = r.S - g.w;
    if (f1 >= g.f - 1e-6) return false;
    out = g;
    out.f = f1;
    out.sq.clear();
    for (int i = 0; i < squeeze::count(r); i++) out.sq.push_back({r.v[3 * i], r.v[3 * i + 1], r.v[3 * i + 2]});
    tips_of(out.sq, g.w, g.h, &out.ybr, &out.xtl);
    out.qb = gadget::bin_of(out.ybr); out.qt = gadget::bin_of(out.xtl);
    return true;
}

// neighbour-key cores from the current gadget table (best gadget of each adjacent key) used as FTMC seeds
std::vector<gftmc::CoreSeed> alts_for(const gadget::Gadget& g, const std::vector<gadget::Gadget>& G, const std::vector<Donor>& donors,
                                      const std::unordered_map<uint64_t, int>& best_of) {
    std::vector<gftmc::CoreSeed> out;
    for (int dw = -1; dw <= 1; dw++) for (int dh = -1; dh <= 1; dh++) for (int dd = -1; dd <= 1; dd++) {
        if (!dw && !dh && !dd) continue;
        const int w = g.w + dw, h = g.h + dh, dm = g.dm + dd;
        if (w < 1 || h < 1 || dm < 1) continue;
        auto it = best_of.find(((uint64_t)w * 64 + h) * 128 + dm);
        if (it == best_of.end()) continue;
        const gadget::Gadget& n = G[it->second];
        std::vector<Square> sq;
        double fn = n.f;
        if (!n.sq.empty()) sq = n.sq; else if (n.donor >= 0) gadget::crop(donors[n.donor], n.sym, n.w, n.h, n.dm, sq, &fn); else continue;
        gftmc::CoreSeed cs;
        for (const Square& q : sq) if (gftmc::is_tilted(q)) cs.T.push_back(q);
        if (cs.T.empty()) continue;
        cs.W = n.w + fn; cs.H = n.h + fn;
        out.push_back(std::move(cs));
    }
    return out;
}

// ---- "what is missing" analysis: single gadget keys that would complete a chain at the leaders' overhang ----
struct KeyNeed { double weight = 0; double goal = 1e9; int puzzles = 0; };

std::map<uint64_t, KeyNeed> missing_keys(const std::vector<gadget::Gadget>& G, const gadget::Dp& dpt) {
    std::map<uint64_t, KeyNeed> need;
    // unique (w,h,dm) -> min f over all gadgets
    std::map<uint64_t, double> keyf;
    for (const auto& g : G) { const uint64_t k = ((uint64_t)g.w * 64 + g.h) * 128 + g.dm; auto it = keyf.find(k); if (it == keyf.end() || g.f < it->second) keyf[k] = g.f; }
    constexpr int W1 = gadget::kWMax + 1, D1 = gadget::kDMax + 1;
    for (int n = 200; n <= 1000; n++) {
        const int K = (int)std::floor(std::sqrt((double)n) + 1e-9), D = n - K * K;
        if (D == 0 || K > gadget::kWMax || D > gadget::kDMax) continue;
        auto lit = leaders().find(n);
        if (lit == leaders().end()) continue;
        const double ft = lit->second - K + 1e-6;
        if (ft <= 0.05) continue;
        {   std::vector<gadget::Link> lk; const double fd = gadget::best_chain(dpt, G, K, D, lk);
            if (fd < 0 || K + fd - lit->second < env_d("GQ_NEED_MINGAP", 0.0)) continue;
            if (fd > env_d("GQ_NEED_FMAX", 2.0) || fd < env_d("GQ_NEED_FMIN", 0.0)) continue; }
        // gadget keys usable at this goal
        std::vector<std::array<int, 3>> A;
        for (auto& kv : keyf) if (kv.second <= ft) { const int w = (int)(kv.first / 64 / 128), h = (int)((kv.first / 128) % 64), dm = (int)(kv.first % 128); if (w <= K && h <= K) A.push_back({w, h, dm}); }
        auto idx = [&](int cw, int ch, int d) { return ((size_t)cw * W1 + ch) * D1 + d; };
        std::vector<char> R((size_t)W1 * W1 * D1, 0), B((size_t)W1 * W1 * D1, 0);
        R[idx(0, 0, 0)] = 1;
        for (int cw = 0; cw <= K; cw++) for (int ch = 0; ch <= K; ch++) for (int d = 0; d <= D; d++) {
            if (!R[idx(cw, ch, d)]) continue;
            for (const auto& a : A) for (int mode = 0; mode < 3; mode++) {
                if (cw == 0 && mode) continue;
                const int tw = cw + a[0] + (mode == 1), th = ch + a[1] + (mode == 2), td = d + a[2];
                if (tw > K || th > K || td > D) continue;
                R[idx(tw, th, td)] = 1;
            }
        }
        bool solved = false;
        for (int cw = 0; cw <= K && !solved; cw++) for (int ch = 0; ch <= K; ch++) if (R[idx(cw, ch, D)]) { solved = true; break; }
        if (solved) continue;   // already completable at this overhang (the DP will take it)
        for (int cw = K; cw >= 0; cw--) for (int ch = K; ch >= 0; ch--) for (int d = D; d >= 0; d--) {
            bool b = d == D;
            if (!b) for (const auto& a : A) { for (int mode = 0; mode < 3 && !b; mode++) { if (cw == 0 && mode) continue; const int tw = cw + a[0] + (mode == 1), th = ch + a[1] + (mode == 2), td = d + a[2]; if (tw > K || th > K || td > D) continue; if (B[idx(tw, th, td)]) b = true; } if (b) break; }
            B[idx(cw, ch, d)] = b;
        }
        std::vector<size_t> rs;
        for (int cw = 0; cw <= K; cw++) for (int ch = 0; ch <= K; ch++) for (int d = 0; d <= D; d++) if (R[idx(cw, ch, d)]) rs.push_back(idx(cw, ch, d));
        std::set<uint64_t> suff;
        for (const auto& kv : keyf) {
            if (kv.second <= ft) continue;                       // already usable (is in A)
            const int w = (int)(kv.first / 64 / 128), h = (int)((kv.first / 128) % 64), dm = (int)(kv.first % 128);
            if (w > K || h > K || dm > D) continue;
            if (std::min(w, h) > (int)env_d("GQ_NEED_MINDIM", 99)) continue;
            if (w > (int)env_d("GQ_NEED_MAXW", 99) || h > (int)env_d("GQ_NEED_MAXH", 99) || w < (int)env_d("GQ_NEED_MINW", 1)) continue;
            bool ok = false;
            for (size_t r : rs) {
                const int d = (int)(r % D1), ch = (int)((r / D1) % W1), cw = (int)(r / D1 / W1);
                for (int mode = 0; mode < 3 && !ok; mode++) {
                    if (cw == 0 && mode) continue;
                    const int tw = cw + w + (mode == 1), th = ch + h + (mode == 2), td = d + dm;
                    if (tw > K || th > K || td > D) continue;
                    if (B[idx(tw, th, td)]) ok = true;
                }
                if (ok) break;
            }
            if (ok) suff.insert(kv.first);
        }
        if (suff.empty()) continue;
        const double gap = (lit->second < 0 ? 0 : 1.0) * 1.0;
        (void)gap;
        for (uint64_t k : suff) { KeyNeed& kn = need[k]; kn.weight += 1.0 / suff.size() * (std::sqrt((double)n) > 0 ? 1.0 : 0); kn.goal = std::min(kn.goal, ft); kn.puzzles++; }
    }
    return need;
}

// ---- lever scan: for each candidate key, how much puzzle score would we gain if its overhang dropped by `delta`? ----
void lever_scan(State& s) {
    const double delta = env_d("GQ_LEVER_DELTA", 0.03);
    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    std::map<int, double> ours;   // n -> our side, from the newest submission (loaded as donors with ordinary books: take min side per n from GQ_OURS)
    {
        std::map<std::pair<int, long long>, Donor> uniq;
        if (const char* path = std::getenv("GQ_OURS")) load_donors(path, uniq);
        for (auto& kv : uniq) { auto it = ours.find(kv.first.first); if (it == ours.end() || kv.second.L < it->second) ours[kv.first.first] = kv.second.L; }
    }
    auto total_gain = [&](const gadget::Dp& dp) {
        double g = 0;
        for (int n = 200; n <= 1000; n++) {
            const int K = (int)std::floor(std::sqrt((double)n) + 1e-9), D = n - K * K;
            if (D == 0 || K > gadget::kWMax || D > gadget::kDMax) continue;
            auto it = ours.find(n); if (it == ours.end()) continue;
            std::vector<gadget::Link> lk; const double f = gadget::best_chain(dp, s.G, K, D, lk);
            if (f >= 0 && K + f < it->second - 1e-7) g += (it->second - (K + f)) / std::sqrt((double)n);
        }
        return g;
    };
    const double base_gain = total_gain(s.dp);
    std::map<uint64_t, double> bestf; std::map<uint64_t, std::vector<int>> members;
    for (int i = 0; i < (int)s.G.size(); i++) { const auto& g = s.G[i]; const int w = std::min(g.w, g.h), h = std::max(g.w, g.h); const uint64_t k = ((uint64_t)w * 64 + h) * 128 + g.dm; members[k].push_back(i); auto it = bestf.find(k); if (it == bestf.end() || g.f < it->second) bestf[k] = g.f; }
    std::vector<std::pair<double, uint64_t>> res;
    const int wmax = (int)env_d("GQ_SW_W", 10), hmax = (int)env_d("GQ_SW_H", 24), dmax = (int)env_d("GQ_SW_DM", 14);
    int done = 0;
    for (auto& kv : bestf) {
        const uint64_t k = kv.first; const int w = (int)(k / 64 / 128), h = (int)((k / 128) % 64), dm = (int)(k % 128);
        if (w > wmax || h > hmax || dm > dmax || kv.second < env_d("GQ_LEVER_FMIN", 0.5) || kv.second > 0.99) continue;
        std::vector<double> saved;
        for (int i : members[k]) { saved.push_back(s.G[i].f); s.G[i].f = std::max(0.05, s.G[i].f - delta); }
        const gadget::Dp dp = gadget::run_dp(s.G, threads);
        const double gain = total_gain(dp) - base_gain;
        size_t q = 0; for (int i : members[k]) s.G[i].f = saved[q++];
        res.push_back({gain, k});
        if (++done % 25 == 0) std::fprintf(stderr, "lever: %d keys scanned\n", done);
    }
    std::sort(res.rbegin(), res.rend());
    std::fprintf(stderr, "lever scan (delta %.3f, baseline gain %.5f): top keys by puzzle-score gain\n", delta, base_gain);
    for (size_t i = 0; i < std::min<size_t>(60, res.size()); i++) { const uint64_t k = res[i].second; std::fprintf(stderr, "LEVER %dx%d+%d f=%.4f gain %.5f\n", (int)(k / 64 / 128), (int)((k / 128) % 64), (int)(k % 128), bestf[k], res[i].first); }
}

void rounds() {
    State& s = S();
    const double budget = env_d("GQ_SECONDS", 1800), tau = env_d("GQ_TAU", 45);
    const int batch = (int)env_d("GQ_BATCH", 40);
    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    const auto t0 = std::chrono::steady_clock::now();
    auto el = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    std::set<int> tried;
    std::vector<gadget::Gadget> saved_gadgets;
    std::map<int, double> goal_of;
    bool pass_improved = false;
    for (int round = 0; el() < budget; round++) {
        std::unordered_map<uint64_t, int> best_of;
        for (int i = 0; i < (int)s.G.size(); i++) { const auto& g = s.G[i]; const uint64_t k = ((uint64_t)g.w * 64 + g.h) * 128 + g.dm; auto it = best_of.find(k); if (it == best_of.end() || g.f < s.G[it->second].f) best_of[k] = i; }
        std::map<int, double> weight;
        double total_gap = 0;
        std::vector<std::pair<double, int>> order;
        if (std::getenv("GQ_NEED")) {
            auto need = missing_keys(s.G, s.dp);
            std::vector<std::pair<double, uint64_t>> ranked;
            for (auto& kv : need) if (!tried.count((int)kv.first)) ranked.push_back({kv.second.weight, kv.first});
            std::sort(ranked.rbegin(), ranked.rend());
            std::fprintf(stderr, "gquench: %zu single-gadget keys would complete some chain; top: ", need.size());
            for (size_t i = 0; i < std::min<size_t>(std::getenv("GQ_NEED_ONLY") ? 60 : 8, ranked.size()); i++) { const uint64_t k = ranked[i].second; std::fprintf(stderr, "%dx%d+%d(w%.1f,f<=%.3f,%dp) ", (int)(k / 64 / 128), (int)((k / 128) % 64), (int)(k % 128), ranked[i].first, need[k].goal, need[k].puzzles); }
            std::fprintf(stderr, "\n");
            if (std::getenv("GQ_NEED_ONLY")) break;
            for (auto& r : ranked) {
                if ((int)order.size() >= batch) break;
                auto it = best_of.find(r.second);
                if (it == best_of.end()) continue;
                order.push_back({r.first, it->second});
                goal_of[it->second] = need[r.second].goal;
                tried.insert((int)r.second);
            }
        } else if (std::getenv("GQ_KEYS")) {
            std::stringstream ks(std::getenv("GQ_KEYS")); std::string tok;
            const int repl = (int)env_d("GQ_REPL", 8);
            while (std::getline(ks, tok, ';')) {
                int w, h, dm; double gf = -1; if (std::sscanf(tok.c_str(), "%d,%d,%d,%lf", &w, &h, &dm, &gf) < 3) continue;
                auto it = best_of.find(((uint64_t)w * 64 + h) * 128 + dm);
                if (it == best_of.end()) { std::fprintf(stderr, "key %s not in table\n", tok.c_str()); continue; }
                if (gf > 0) goal_of[it->second] = gf;
                for (int r = 0; r < repl; r++) order.push_back({s.G[it->second].f, it->second});
            }
        } else if (std::getenv("GQ_SWEEP")) {
            const int wmax = (int)env_d("GQ_SW_W", 8), hmax = (int)env_d("GQ_SW_H", 24), dmax = (int)env_d("GQ_SW_DM", 8), wmin = (int)env_d("GQ_SW_WMIN", 2);
            for (auto& kv : best_of) {
                const auto& g = s.G[kv.second];
                if (g.w < wmin || g.w > wmax || g.h > hmax || g.dm > dmax || g.h < (int)env_d("GQ_SW_HMIN", 3)) continue;
                if (g.f < env_d("GQ_SW_FMIN", 0.45) || g.f > env_d("GQ_SW_FMAX", 1.0)) continue;                 // already very tight
                if (tried.count((int)kv.first)) continue;
                order.push_back({g.f, kv.second});
            }
            std::sort(order.rbegin(), order.rend());
            if ((int)order.size() > batch) order.resize(batch);
            for (auto& o : order) { const auto& g = s.G[o.second]; tried.insert((int)(((uint64_t)g.w * 64 + g.h) * 128 + g.dm)); }
        } else {
        for (int n = 100; n <= 1000; n++) {
            const int K = (int)std::floor(std::sqrt((double)n) + 1e-9), D = n - K * K;
            if (D == 0 || K > gadget::kWMax || D > gadget::kDMax) continue;
            auto it = leaders().find(n);
            if (it == leaders().end()) continue;
            std::vector<gadget::Link> links;
            const double f = gadget::best_chain(s.dp, s.G, K, D, links);
            if (f < 0) continue;
            const double gap = (K + f - it->second) / std::sqrt((double)n);
            if (gap <= 1e-6) continue;
            if (f < env_d("GQ_CHAIN_FMIN", 0.0) || f > env_d("GQ_CHAIN_FMAX", 2.0)) continue;
            total_gap += gap;
            int nb = 0;
            auto small = [&](const gadget::Gadget& g) { return !(g.w >= K - 1 && g.h >= K - 1); };
            for (const auto& l : links) if (s.G[l.g].f >= f - 0.03 && small(s.G[l.g])) nb++;
            for (const auto& l : links) if (s.G[l.g].f >= f - 0.03 && small(s.G[l.g])) weight[l.g] += gap / std::max(1, nb);
        }
        for (auto& kv : weight) if (!tried.count(kv.first)) order.push_back({kv.second, kv.first});
        std::sort(order.rbegin(), order.rend());
        if ((int)order.size() > batch) order.resize(batch);
        }
        if (std::getenv("GQ_KEYS")) { tried.clear(); }
        if (order.empty()) {
            if (std::getenv("GQ_SWEEP") && pass_improved) { tried.clear(); pass_improved = false; std::fprintf(stderr, "gquench: sweep pass complete, starting another\n"); continue; }
            std::fprintf(stderr, "gquench: nothing left to try\n"); break;
        }
        if ((int)order.size() > batch) order.resize(batch);
        std::fprintf(stderr, "gquench round %d [%.0fs]: total DP gap to leaders %.5f, quenching %zu gadgets (top weight %.5f)\n", round, el(), total_gap, order.size(), order[0].first);
        std::atomic<int> next{0};
        std::mutex mu;
        std::vector<gadget::Gadget> improved;
        auto worker = [&](int tid) {
            for (;;) {
                const int i = next++;
                if (i >= (int)order.size()) return;
                const gadget::Gadget& g = s.G[order[i].second];
                gadget::Gadget out;
                static const Donor kNoDonor;
                if (std::getenv("GQ_TRACE")) { std::lock_guard<std::mutex> lk(mu); std::fprintf(stderr, "  start %dx%d+%d f=%.5f donor %d sym %d sq %zu\n", g.w, g.h, g.dm, g.f, g.donor, g.sym, g.sq.size()); }
                const Donor& dd = g.donor >= 0 ? s.donors[g.donor] : kNoDonor;
                bool ok = false;
                try { ok = std::getenv("GQ_FTMC") ? false : quench(g, dd, tau, 1000 + i * 17 + tid, out); } catch (const std::exception& e) { std::fprintf(stderr, "  quench of %dx%d+%d failed: %s\n", g.w, g.h, g.dm, e.what()); }
                if (!ok && std::getenv("GQ_FTMC")) try { auto alts = alts_for(g, s.G, s.donors, best_of); ok = gftmc::ftmc_gadget(g, dd, tau, 5000 + i * 31 + tid, out, &alts, goal_of.count(order[i].second) ? goal_of[order[i].second] : -1.0); } catch (const std::exception& e) { std::fprintf(stderr, "  ftmc failed: %s\n", e.what()); }
                std::lock_guard<std::mutex> lk(mu);
                std::fprintf(stderr, "  gadget %dx%d+%d f=%.5f (weight %.5f): %s\n", g.w, g.h, g.dm, g.f, order[i].first, ok ? ("-> " + std::to_string(out.f)).c_str() : "unchanged");
                if (ok) improved.push_back(out);
            }
        };
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; t++) pool.emplace_back(worker, t);
        for (auto& t : pool) t.join();
        if (!std::getenv("GQ_SWEEP")) for (const auto& o : order) tried.insert(o.second);
        if (!improved.empty() && round % 3 == 2 && !std::getenv("GQ_SWEEP")) tried.clear();   // neighbours changed: let earlier targets try again
        if (improved.empty()) continue;
        pass_improved = true;
        {   std::map<uint64_t, gadget::Gadget> bestk;
            for (auto& g : improved) { const uint64_t k = ((uint64_t)g.w * 64 + g.h) * 128 + g.dm; auto it = bestk.find(k); if (it == bestk.end() || g.f < it->second.f) bestk[k] = g; }
            improved.clear(); for (auto& kv : bestk) improved.push_back(kv.second); }
        {   // improved gadgets: keep them (8 variants) AND use them as donors so every smaller corner window becomes a gadget too
            std::vector<Donor> nd;
            for (auto& g : improved) {
                saved_gadgets.push_back(g);
                for (auto& v : variants_of(g.w, g.h, g.dm, g.sq)) s.G.push_back(std::move(v));
                Donor d; d.n = (int)g.sq.size(); d.L = std::max(g.w, g.h) + g.f;
                for (const Square& q : g.sq) { d.x.push_back(q.x); d.y.push_back(q.y); d.t.push_back(q.angle); }
                nd.push_back(std::move(d));
            }
            std::vector<gadget::Gadget> crops = gadget::extract(nd, threads);
            const int off = (int)s.donors.size();
            for (auto& d : nd) s.donors.push_back(std::move(d));
            for (auto& g : crops) { g.donor += off; s.G.push_back(std::move(g)); }
            std::fprintf(stderr, "gquench: crop closure added %zu gadgets\n", crops.size());
        }
        save_gadgets(saved_gadgets, std::getenv("GQ_SAVE") ? std::getenv("GQ_SAVE") : "work/gq_gadgets.txt");
        s.dp = gadget::run_dp(s.G, threads);
        std::fprintf(stderr, "gquench: %zu gadgets improved, DP recomputed [%.0fs]\n", improved.size(), el());
    }
}
}  // namespace

Packing solve(int n, const Packing& base) {
    const int K = (int)std::floor(std::sqrt((double)n) + 1e-9), D = n - K * K;
    if (D == 0 || K > gadget::kWMax || D > gadget::kDMax || base.empty()) return base;
    static bool quenched = false;
    init();
    if (!quenched) {
        quenched = true;
        if (const char* dk = std::getenv("GQ_DUMP")) {   // GQ_DUMP="w,h,dm;..." -> work/dump_<w>x<h>+<dm>.txt (best gadget of the key, explicit squares)
            std::stringstream ks(dk); std::string tok;
            while (std::getline(ks, tok, ';')) {
                int w, h, dm; if (std::sscanf(tok.c_str(), "%d,%d,%d", &w, &h, &dm) != 3) continue;
                int bi = -1; for (int i = 0; i < (int)S().G.size(); i++) { const auto& g = S().G[i]; if (g.w == w && g.h == h && g.dm == dm && (bi < 0 || g.f < S().G[bi].f)) bi = i; }
                if (bi < 0) { std::fprintf(stderr, "dump: key %s missing\n", tok.c_str()); continue; }
                gadget::Gadget g = S().G[bi]; std::vector<Square> sq; double f = g.f;
                if (!g.sq.empty()) sq = g.sq; else gadget::crop(S().donors[g.donor], g.sym, g.w, g.h, g.dm, sq, &f);
                char name[128]; std::snprintf(name, sizeof name, "work/dump_%dx%d+%d.txt", w, h, dm);
                std::FILE* o = std::fopen(name, "w"); std::fprintf(o, "%d %d %d %.12f\\n", w, h, dm, g.f);
                for (const Square& q : sq) std::fprintf(o, "%.12f %.12f %.12f\\n", q.x, q.y, q.angle); std::fclose(o);
                std::fprintf(stderr, "dump: %s f=%.5f\\n", name, g.f);
            }
            return base;
        }
        if (std::getenv("GQ_LEVER")) lever_scan(S()); else rounds();
    }
    State& s = S();
    std::vector<gadget::Link> links;
    const double f = gadget::best_chain(s.dp, s.G, K, D, links);
    const double L0 = container_side(base);
    if (f < 0 || K + f + 2e-9 >= L0 - 1e-7) return base;
    Packing p = gadget::compose(K, D, f, links, s.G, s.donors);
    if (p.empty()) { std::fprintf(stderr, "n=%d: DP says %.9f < %.9f but compose failed\n", n, K + f, L0); return base; }
    std::fprintf(stderr, "n=%d: DP side %.9f vs ours %.9f (%zu gadgets)\n", n, container_side(p), L0, links.size());
    { long long tight = 0; int bi = -1, bj = -1; if (!packing_valid(p, &tight, &bi, &bj)) { std::fprintf(stderr, "   INVALID pair %d,%d: (%.9f,%.9f,%.5f) vs (%.9f,%.9f,%.5f)\n", bi, bj, bi >= 0 ? p[bi].x : 0, bi >= 0 ? p[bi].y : 0, bi >= 0 ? p[bi].angle : 0, bj >= 0 ? p[bj].x : 0, bj >= 0 ? p[bj].y : 0, bj >= 0 ? p[bj].angle : 0); for (const auto& l : links) { const auto& g = s.G[l.g]; std::fprintf(stderr, "   link %dx%d+%d f=%.6f mode%d sym%d\n", g.w, g.h, g.dm, g.f, l.mode, g.sym); } } }
    return p;
}
