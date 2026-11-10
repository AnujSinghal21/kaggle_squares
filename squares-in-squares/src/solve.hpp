#pragma once
#include "common.hpp"

// Return a packing of n unit squares (n in 1..kMaxN), or an empty Packing for "no answer for this n".
// `base` is the best packing of n so far (the newest submission), or empty if there is none.
// Returning `base` (or anything not strictly better) leaves the submission unchanged for that n.
// Each approach lives in its own src/solvers/solve_<name>.cpp; exactly one is linked per binary.
Packing solve(int n, const Packing& base);
