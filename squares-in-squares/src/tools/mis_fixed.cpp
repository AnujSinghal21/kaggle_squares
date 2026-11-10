#include <cstdio>
#include <cstdlib>
#include "../mis.hpp"
int main(int argc,char**argv){
  double S=atof(argv[1]),th=atof(argv[2]),ox=atof(argv[3]),oy=atof(argv[4]),tol=atof(argv[5]);
  std::vector<mis::Cand> c; mis::add_axis_lattices(c,S); mis::add_lattice(c,S,th,ox,oy);
  auto g=mis::build(c,tol); std::mt19937_64 rng(1); int best=0;
  for(int r=0;r<20;r++){ auto s=mis::solve(g,300,rng); best=std::max<int>(best,s.size()); }
  printf("cands %zu best %d\n",g.c.size(),best);
}
