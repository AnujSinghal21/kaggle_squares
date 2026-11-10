// Submission CSV I/O and a local re-implementation of the Kaggle metric: parse a CSV the way pandas
// would (float64), check overlaps (common.hpp), and report the score the leaderboard would give.
#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "common.hpp"

// All puzzles of a submission, indexed by n (1..kMaxN). A puzzle is absent if no row mentioned it.
struct Submission {
    std::vector<Packing> pk = std::vector<Packing>(kMaxN + 1);
    std::vector<char> present = std::vector<char>(kMaxN + 1, 0);
    std::vector<std::string> broken = std::vector<std::string>(kMaxN + 1);  // non-empty = parse error reason
    long long rows = 0, bad_rows = 0;
    std::string fatal;
};

inline Submission load_submission(const std::string& path) {
    Submission sub;
    std::ifstream in(path);
    if (!in) { sub.fatal = "cannot open " + path; return sub; }
    std::string line;
    std::getline(in, line);
    if (line.rfind("id,n,square_id,x,y,angle", 0) != 0) { sub.fatal = "unexpected header: " + line; return sub; }

    std::vector<std::vector<char>> seen(kMaxN + 1);
    auto fail = [&](int n, const char* msg) { if (sub.broken[n].empty()) sub.broken[n] = msg; };

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        sub.rows++;
        char* e = nullptr;
        std::strtoll(line.c_str(), &e, 10);  // id (uniqueness is not checked locally)
        if (*e != ',') { sub.bad_rows++; continue; }
        const long n = std::strtol(e + 1, &e, 10);
        if (*e != ',' || n < 1 || n > kMaxN) { sub.bad_rows++; continue; }
        sub.present[n] = 1;
        const long sid = std::strtol(e + 1, &e, 10);
        if (*e != ',') { fail((int)n, "unparseable row"); continue; }
        const double x = std::strtod(e + 1, &e);
        if (*e != ',') { fail((int)n, "unparseable row"); continue; }
        const double y = std::strtod(e + 1, &e);
        if (*e != ',') { fail((int)n, "unparseable row"); continue; }
        const double a = std::strtod(e + 1, &e);
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(a)) { fail((int)n, "non-finite value"); continue; }
        if (sid < 0 || sid >= n) { fail((int)n, "square_id out of range"); continue; }
        if (seen[n].empty()) {
            seen[n].assign(n, 0);
            sub.pk[n].assign(n, {0, 0, 0});
        }
        if (seen[n][sid]) { fail((int)n, "duplicate square_id"); continue; }
        seen[n][sid] = 1;
        sub.pk[n][sid] = {x, y, a};
    }
    for (int n = 1; n <= kMaxN; n++)
        if (sub.present[n] && sub.broken[n].empty() &&
            (seen[n].empty() || std::find(seen[n].begin(), seen[n].end(), 0) != seen[n].end()))
            sub.broken[n] = "wrong number of squares";
    return sub;
}

// Rows are written grouped by n, ids sequential; 17 digits round-trip a double exactly.
inline void write_submission(const std::string& path, const Submission& sub) {
    std::ofstream out(path);
    out << std::setprecision(17) << "id,n,square_id,x,y,angle\n";
    long long id = 0;
    for (int n = 1; n <= kMaxN; n++) {
        if (!sub.present[n] || !sub.broken[n].empty()) continue;
        for (int i = 0; i < (int)sub.pk[n].size(); i++) {
            const Square& s = sub.pk[n][i];
            out << id++ << ',' << n << ',' << i << ',' << s.x << ',' << s.y << ',' << s.angle << '\n';
        }
    }
}

struct PuzzleResult {
    enum Status { Missing, Invalid, Valid } status = Missing;
    double L = 0.0;
    long long tight = 0;  // risky near-contacts, see common.hpp
    std::string reason;
};

struct Report {
    std::vector<PuzzleResult> puzzles = std::vector<PuzzleResult>(kMaxN + 1);  // index n
    double total = 0.0, grid_total = 0.0;
    int valid = 0, invalid = 0, missing = 0;
    long long rows = 0, bad_rows = 0;
    std::string fatal;
};

inline Report validate_submission(const Submission& sub) {
    Report r;
    r.rows = sub.rows;
    r.bad_rows = sub.bad_rows;
    r.fatal = sub.fatal;
    for (int n = 1; n <= kMaxN; n++) r.grid_total += grid_score(n);
    if (!r.fatal.empty()) return r;

    for (int n = 1; n <= kMaxN; n++) {
        PuzzleResult& pr = r.puzzles[n];
        if (!sub.present[n]) { pr.status = PuzzleResult::Missing; r.missing++; r.total += grid_score(n); continue; }
        std::string why = sub.broken[n];
        if (why.empty()) {
            int i = 0, j = 0;
            if (!packing_valid(sub.pk[n], &pr.tight, &i, &j))
                why = "squares " + std::to_string(i) + " and " + std::to_string(j) + " overlap";
        }
        if (!why.empty()) { pr.status = PuzzleResult::Invalid; pr.reason = why; r.invalid++; r.total += grid_score(n); continue; }
        pr.status = PuzzleResult::Valid;
        pr.L = container_side(sub.pk[n]);
        r.valid++;
        r.total += score_of(n, pr.L);
    }
    return r;
}

inline Report validate_csv(const std::string& path) { return validate_submission(load_submission(path)); }

inline void print_report(const Report& r, const std::string& path) {
    if (!r.fatal.empty()) { std::printf("validation failed: %s\n", r.fatal.c_str()); return; }
    std::printf("validated %s\n", path.c_str());
    std::printf("  puzzles: %d valid, %d invalid, %d missing (charged grid ceiling)\n", r.valid, r.invalid, r.missing);
    int shown = 0;
    long long tight_puzzles = 0, tight_pairs = 0;
    for (int n = 1; n <= kMaxN; n++) {
        const PuzzleResult& pr = r.puzzles[n];
        if (pr.status == PuzzleResult::Invalid && shown++ < 10) std::printf("  INVALID n=%d: %s\n", n, pr.reason.c_str());
        if (pr.status == PuzzleResult::Valid && pr.tight) { tight_puzzles++; tight_pairs += pr.tight; }
    }
    if (r.invalid > 10) std::printf("  ... and %d more invalid\n", r.invalid - 10);
    if (r.bad_rows) std::printf("  warning: %lld rows with a bad n/id were ignored\n", r.bad_rows);
    if (tight_pairs)
        std::printf("  warning: %lld puzzles have %lld risky near-contacts (|gap| < 1e-15: Kaggle may reject these; dilate to open a real gap)\n", tight_puzzles, tight_pairs);
    std::printf("  grid baseline total : %.10f\n", r.grid_total);
    std::printf("  SCORE (lower=better): %.10f  (%+.10f vs grid)\n", r.total, r.total - r.grid_total);
}
