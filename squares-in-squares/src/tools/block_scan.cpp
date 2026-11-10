// usage: block_scan S [th_lo th_hi th_step] [wh_lo wh_hi] [diag_only=1]
// scans block tilt, size (W x H) and centre; prints the best total count and parameters (top few)
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <thread>
#include "../block.hpp"
int main(int argc, char** argv) {
    const double S = atof(argv[1]);
    const double th_lo = argc > 2 ? atof(argv[2]) : 40, th_hi = argc > 3 ? atof(argv[3]) : 50, th_step = argc > 4 ? atof(argv[4]) : 0.25;
    const int wl = argc > 5 ? atoi(argv[5]) : 8, wh = argc > 6 ? atoi(argv[6]) : 14;
    const int diag = argc > 7 ? atoi(argv[7]) : 1;
    const double tol = argc > 8 ? atof(argv[8]) : 0.02;
    struct Job { double th; int W, H; };
    std::vector<Job> jobs;
    for (double th = th_lo; th <= th_hi + 1e-9; th += th_step) for (int W = wl; W <= wh; W++) for (int H = W; H <= wh; H++) jobs.push_back({th, W, H});
    std::atomic<size_t> next{0}; std::mutex mu; int best = 0;
    auto worker = [&]() {
        for (;;) {
            size_t i = next++; if (i >= jobs.size()) return;
            const Job& j = jobs[i];
            const double ext = (j.W + j.H) * 0.5 * 0.72;
            const double cstep = diag ? 0.02 : 0.04; const double c0 = diag ? 0.25 : 0.35, c1 = diag ? 0.75 : 0.65;
            for (double cx = S * c0; cx < S * c1; cx += cstep)
                for (double cy = diag ? cx : S * c0; cy < (diag ? cx + 1e-9 : S * c1); cy += cstep) {
                    (void)ext;
                    block::Params P; P.S = S; P.theta = j.th * M_PI / 180; P.cx = cx; P.cy = cy; P.W = j.W; P.H = j.H; P.tol = tol;
                    auto B = block::build(P);
                    if (B.n_block < 0) continue;
                    std::lock_guard<std::mutex> lk(mu);
                    if (B.total() > best) { best = B.total(); std::fprintf(stderr, "count %d theta %.2f c=(%.2f,%.2f) %dx%d\n", best, j.th, cx, cy, j.W, j.H); }
                }
        }
    };
    std::vector<std::thread> pool; for (int t = 0; t < 8; t++) pool.emplace_back(worker); for (auto& t : pool) t.join();
    printf("S=%.5f best %d\n", S, best);
}
