// usage: block_probe S theta_deg cx cy W H [tol]   -> count
#include <cstdio>
#include <cstdlib>
#include "../block.hpp"
int main(int argc, char** argv) {
    block::Params P; P.S = atof(argv[1]); P.theta = atof(argv[2]) * M_PI / 180; P.cx = atof(argv[3]); P.cy = atof(argv[4]); P.W = atoi(argv[5]); P.H = atoi(argv[6]);
    if (argc > 7) P.tol = atof(argv[7]);
    auto B = block::build(P);
    printf("block %d total %d\n", B.n_block, B.total());
}
