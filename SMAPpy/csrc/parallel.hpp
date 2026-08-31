// Splitting per-frame work over threads.
//
// Every stage here treats frames independently, so parallelising means handing
// each thread a contiguous range of frames.  No locks and no shared state: the
// only coordination is the join at the end.
#pragma once

#include <algorithm>
#include <functional>
#include <thread>
#include <vector>

namespace smappy {

// Resolve a thread count: 0 or negative means "one per core".
inline int resolve_threads(int requested, long long items) {
    if (requested <= 0) {
        requested = static_cast<int>(std::thread::hardware_concurrency());
        if (requested <= 0) requested = 1;
    }
    if (items > 0 && requested > items) requested = static_cast<int>(items);
    return std::max(requested, 1);
}

// Run `body(begin, end, thread_index)` over contiguous chunks of [0, n).
template <class Body>
void parallel_ranges(long long n, int n_threads, Body body) {
    n_threads = resolve_threads(n_threads, n);
    if (n_threads == 1 || n <= 1) {
        body(0LL, n, 0);
        return;
    }

    const long long step = (n + n_threads - 1) / n_threads;
    std::vector<std::thread> pool;
    pool.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
        const long long begin = t * step;
        const long long end = std::min(begin + step, n);
        if (begin < end) pool.emplace_back(body, begin, end, t);
    }
    for (auto& th : pool) th.join();
}

}  // namespace smappy
