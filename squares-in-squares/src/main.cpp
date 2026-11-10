// Squares in Squares: pack n unit squares into the smallest axis-aligned square.
// Score per puzzle: L / sqrt(n) - 1 (lower is better). L is derived from the squares.
#include <cctype>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <string>
#include <vector>

#include "closure.hpp"
#include "common.hpp"
#include "solve.hpp"
#include "validate.hpp"

namespace fs = std::filesystem;

static std::string timestamped_path(const std::string& dir) {
    std::time_t t = std::time(nullptr);
    char buf[32];
    std::strftime(buf, sizeof buf, "%Y%m%d_%H%M%S", std::localtime(&t));
    std::string path = dir + "/submission_" + buf + ".csv";
    for (int i = 2; fs::exists(path); i++)  // never overwrite a same-second run
        path = dir + "/submission_" + buf + "_" + std::to_string(i) + ".csv";
    return path;
}

// Newest *.csv in dir by modification time, or "" if there is none.
static std::string newest_csv(const std::string& dir) {
    std::string best;
    fs::file_time_type best_t;
    if (!fs::exists(dir)) return best;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (e.path().extension() != ".csv") continue;
        const auto t = fs::last_write_time(e);
        if (best.empty() || t > best_t) { best = e.path().string(); best_t = t; }
    }
    return best;
}

static int usage(const char* prog) {
    std::fprintf(stderr,
        "usage: %s [n_min [n_max]] [--base file.csv | --fresh] [--replace]\n"
        "       %s validate <file.csv>\n"
        "  solve n_min..n_max (default 1..%d) and write a timestamped CSV into submissions/.\n"
        "  Starts from the newest file in submissions/ by default: only n_min..n_max is recomputed and\n"
        "  the better of old/new is kept per puzzle, so a worse solve() can never regress you.\n"
        "  --base     start from this file instead of the newest submission.\n"
        "  --fresh    ignore existing submissions; the output holds only what is solved now.\n"
        "  --replace  always take the newly solved packing, even if it scores worse.\n",
        prog, prog, kMaxN);
    return 1;
}

int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "validate") {
        if (argc < 3) return usage(argv[0]);
        const Report rep = validate_csv(argv[2]);
        print_report(rep, argv[2]);
        return rep.fatal.empty() && rep.invalid == 0 ? 0 : 2;
    }

    std::string base_path;
    bool replace = false, fresh = false;
    int nums[2], num_count = 0;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        if (a == "--base" && i + 1 < argc) base_path = argv[++i];
        else if (a == "--replace") replace = true;
        else if (a == "--fresh") fresh = true;
        else if (num_count < 2 && !a.empty() && std::isdigit((unsigned char)a[0])) nums[num_count++] = std::atoi(a.c_str());
        else return usage(argv[0]);
    }
    const int n_min = num_count > 0 ? nums[0] : 1;
    const int n_max = num_count > 1 ? nums[1] : (num_count == 1 ? n_min : kMaxN);
    if (n_min < 1 || n_max > kMaxN || n_min > n_max) return usage(argv[0]);

    if (fresh && !base_path.empty()) return usage(argv[0]);
    if (!fresh && base_path.empty()) base_path = newest_csv("submissions");  // may stay empty

    Submission sub;  // stays empty without a base: the output holds only what we solve
    if (!base_path.empty()) {
        sub = load_submission(base_path);
        if (!sub.fatal.empty()) { std::fprintf(stderr, "base: %s\n", sub.fatal.c_str()); return 1; }
        const Report before = validate_submission(sub);
        std::printf("base %s: score %.10f (%d valid, %d invalid, %d missing)\n", base_path.c_str(), before.total,
                    before.valid, before.invalid, before.missing);
    }

    struct Adopted { int n; double old_L, new_L, delta; };  // delta = score drop (positive = better)
    std::vector<Adopted> adopted;
    const auto t0 = std::chrono::steady_clock::now();
    int kept = 0, rejected = 0, none = 0;
    for (int n = n_min; n <= n_max; n++) {
        const bool have_old = sub.present[n] && sub.broken[n].empty() && packing_valid(sub.pk[n]);
        const Packing no_base;
        Packing p = solve(n, have_old ? sub.pk[n] : no_base);
        if (p.empty()) { none++; continue; }  // this approach has no answer for n
        long long risky = 0;
        if ((int)p.size() != n || !packing_valid(p, &risky)) { std::fprintf(stderr, "n=%d: solve() returned an invalid packing, skipped\n", n); rejected++; continue; }
        if (risky) { std::fprintf(stderr, "n=%d: %lld risky near-contacts (Kaggle may reject), skipped\n", n, risky); rejected++; continue; }
        const double new_L = container_side(p);
        const double old_L = have_old ? container_side(sub.pk[n]) : NAN;
        if (have_old && !replace && score_of(n, old_L) <= score_of(n, new_L)) { kept++; continue; }
        adopted.push_back({n, old_L, new_L, have_old ? score_of(n, old_L) - score_of(n, new_L) : NAN});
        sub.pk[n] = std::move(p);
        sub.present[n] = 1;
        sub.broken[n].clear();
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::printf("solved n=%d..%d in %.2fs: %zu adopted, %d kept from base, %d no answer, %d rejected\n", n_min, n_max, secs,
                adopted.size(), kept, none, rejected);
    if (!adopted.empty()) {  // the n where this approach's packing was taken, biggest gain first
        std::sort(adopted.begin(), adopted.end(), [](const Adopted& x, const Adopted& y) {
            return std::isnan(x.delta) != std::isnan(y.delta) ? std::isnan(y.delta) : x.delta > y.delta;
        });
        double total = 0;
        std::printf("  %-6s %-14s %-14s %s\n", "n", "L before", "L after", "score delta (negative = better)");
        for (size_t i = 0; i < adopted.size(); i++) {
            const Adopted& a = adopted[i];
            if (!std::isnan(a.delta)) total += a.delta;
            if (i < 40) {
                if (std::isnan(a.delta)) std::printf("  %-6d %-14s %-14.9f new (score %.10f)\n", a.n, "-", a.new_L, a.new_L / std::sqrt((double)a.n) - 1);
                else std::printf("  %-6d %-14.9f %-14.9f -%.10f\n", a.n, a.old_L, a.new_L, a.delta);
            }
        }
        if (adopted.size() > 40) std::printf("  ... and %zu more\n", adopted.size() - 40);
        if (total > 0) std::printf("  total score delta vs base: -%.10f\n", total);
    }
    // Closure pass over all 1..1000, so a gain at n also reaches every smaller n (and add-L reaches larger ones).
    const std::vector<ClosureChange> closed = closure_pass(sub);
    if (!closed.empty()) {
        double total = 0;
        std::printf("closure: %zu puzzles improved by deleting from n+1 or adding an L-frame\n", closed.size());
        std::printf("  %-6s %-14s %-14s %-12s %s\n", "n", "L before", "L after", "delta", "how");
        size_t shown = 0;
        for (const ClosureChange& c : closed) {
            const double d = score_of(c.n, c.old_L) - score_of(c.n, c.new_L);
            total += d;
            if (shown++ < 40) std::printf("  %-6d %-14.9f %-14.9f -%-11.8f %s\n", c.n, c.old_L, c.new_L, d, c.how.c_str());
        }
        if (closed.size() > 40) std::printf("  ... and %zu more\n", closed.size() - 40);
        std::printf("  closure score delta: -%.10f\n", total);
    }
    if (adopted.empty() && closed.empty() && !base_path.empty()) {
        std::printf("no change over the base, nothing written\n");
        return 0;
    }

    fs::create_directories("submissions");
    const std::string path = timestamped_path("submissions");
    write_submission(path, sub);

    const Report rep = validate_csv(path);  // re-read from disk: same float64 round trip as Kaggle
    std::printf("wrote %s\n", path.c_str());
    print_report(rep, path);
    return rep.fatal.empty() && rep.invalid == 0 ? 0 : 2;
}
