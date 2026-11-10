// Approach: smooth squeeze + basin hopping (src/squeeze.hpp) started from the best packing so far.
// RESULT (2026-10-27): polishes only. Rigid strip/slab packings do not move (n=299/978 unchanged); 8e-7..1e-6 gains on numeric puzzles; polish sweep over 237 trailing puzzles gave ~1e-4 total. Benchmark from a plain grid (SQ_GRID=1): reaches most known optima for n<=37 given ~150 s, misses n=40 (6.88 vs 6.83).
//   Idea:    shrink the box with a penetration-energy L-BFGS, kick a few squares when stuck, polish smoothly.
//            Same family as Evan Daniel's `packer hop` / `fq` (github.com/evand/square-packing, MIT) but our own
//            implementation. Edit kTargets / kSeconds below, then:  make run SOLVER=squeeze
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <thread>

#include "../solve.hpp"
#include "../squeeze.hpp"

namespace {
const int kTargets[] = {304, 305};      // the n to work on; everything else keeps the base
const double kSeconds = 60;        // wall-clock budget per n, using all hardware threads
const double kFixedDrop = 0;       // 0 = gradual shrinking; > 0 = hop at start side minus this (e.g. 0.0129 for n=299)
}  // namespace

Packing solve(int n, const Packing& base) {
    bool wanted = false;
    for (int t : kTargets) wanted |= (t == n);
    if (const char* e = std::getenv("SQ_TARGETS")) { wanted = false; std::string tok; std::stringstream ss(e); while (std::getline(ss, tok, ',')) if (std::atoi(tok.c_str()) == n) wanted = true; }
    if (!wanted || base.empty()) return base;
    squeeze::Options o;
    o.seconds = std::getenv("SQ_SECONDS") ? std::atof(std::getenv("SQ_SECONDS")) : kSeconds;
    o.threads = (int)std::max(1u, std::thread::hardware_concurrency());
    o.seed = 12345 + n;
    o.fixed_drop = kFixedDrop;
    o.verbose = kFixedDrop > 0;
    std::fprintf(stderr, "n=%d: squeezing from side %.9f for %.0fs on %d threads\n", n, container_side(base), o.seconds, o.threads);
    Packing startp = std::getenv("SQ_GRID") ? grid_packing(n) : base;   // SQ_GRID=1: benchmark from a plain grid (result is printed, not adopted)
    Packing p = squeeze::run(startp, o);
    if (std::getenv("SQ_GRID")) { std::fprintf(stderr, "n=%d: from grid -> side %.9f\n", n, container_side(p)); return base; }
    std::fprintf(stderr, "n=%d: -> side %.9f\n", n, container_side(p));
    return p;
}
