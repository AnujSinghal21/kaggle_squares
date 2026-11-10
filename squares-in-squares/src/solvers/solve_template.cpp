// Template for a new approach: copy this file to src/solvers/solve_<name>.cpp and edit solve().
//   Source/idea: <cite the paper, repo, notebook or describe your own idea here>
//   License:     <of any data or code you borrow>
//
// Build and run:  make run SOLVER=<name> ARGS="50 60"      (n=50..60 on top of the newest submission)
//                 make run SOLVER=<name> ARGS="--fresh 5"   (ignore the base, see what you produce alone)
// main.cpp adopts a returned packing only if it is valid and scores strictly better than `base`,
// then lists every n it adopted with the score delta, so it is safe to return `base` or {} for the rest.
#include "../solve.hpp"

Packing solve(int n, const Packing& base) {
    // 1. Only touch the n you have an idea for; everything else keeps the base.
    //    if (n < 50 || n > 60) return base;     // or return {} for "no answer"

    // 2. Build a candidate. `base` is empty if the submission has nothing for n; fall back to a grid.
    Packing p = base.empty() ? grid_packing(n) : base;

    // 3. ... improve p here. Helpers in common.hpp: score(p), container_side(p), packing_valid(p).
    //    Rotated squares that touch exactly need a real gap: dilate centres by ~1 + 1e-11 about the middle.

    return p;  // returning the base unchanged is a no-op: nothing is adopted
}
