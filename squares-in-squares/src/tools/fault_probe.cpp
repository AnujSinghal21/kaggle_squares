// usage: fault_probe S [theta_deg|-1 scan] [a_min a_max] [r_min r_max]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <thread>
#include <atomic>
#include "../fault.hpp"
int main(int argc, char** argv) {
    const double S = std::atof(argv[1]);
    const double th_fixed = argc > 2 ? std::atof(argv[2]) : -1;
    const int a0 = argc > 3 ? std::atoi(argv[3]) : 0, a1 = argc > 4 ? std::atoi(argv[4]) : 4;
    const int r0 = argc > 5 ? std::atoi(argv[5]) : 2, r1 = argc > 6 ? std::atoi(argv[6]) : 2;
    struct Job { double th; int a, r; };
    std::vector<Job> jobs;
    if (th_fixed > 0) jobs.push_back({th_fixed, a0, r0});
    for (double th = 12; th <= 45 && th_fixed <= 0; th += 0.25)
        for (int a = a0; a <= a1; a++) for (int r = r0; r <= r1; r++) jobs.push_back({th, a, r});
    std::atomic<size_t> next{0}; std::mutex mu; int best = 0; fault::Result bestR; Job bj{};
    auto worker = [&]() {
        for (;;) {
            size_t i = next++; if (i >= jobs.size()) return;
            fault::Params P; P.S = S; P.max_h = 14; P.theta = jobs[i].th * M_PI / 180; P.a = jobs[i].a; P.b = jobs[i].a - 1 - jobs[i].r;
            auto R = fault::solve(P, false);
            std::lock_guard<std::mutex> lk(mu);
            if (R.total() > best) { best = R.total(); bj = jobs[i]; std::fprintf(stderr, "theta %.2f a=%d r=%d grid %d strip %d total %d\n", bj.th, bj.a, bj.r, R.grid, R.strip, R.total()); }
        }
    };
    std::vector<std::thread> pool; for (int t = 0; t < 8; t++) pool.emplace_back(worker); for (auto& t : pool) t.join();
    std::printf("S=%.5f best total %d (theta %.2f a=%d r=%d)\n", S, best, bj.th, bj.a, bj.r);
}
