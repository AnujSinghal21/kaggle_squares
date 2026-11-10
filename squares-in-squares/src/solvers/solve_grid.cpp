// Baseline: ceil(sqrt(n)) x ceil(sqrt(n)) upright grid. Identical to the sample submission.
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return grid_packing(n); }
