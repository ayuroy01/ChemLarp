#pragma once
// Shared-memory parallelism with results that do not depend on the number of threads.
//
// parallel_for(count, body) runs body(i) for i in [0, count) on a fixed pool of worker threads (the calling
// thread works too). Indices are handed out dynamically, so uneven work balances, but each index must write
// only its own outputs: then the result is the same for any thread count. Reductions use a fixed number of
// partitions chosen by the caller (never by the thread count) and combine the partial results in partition
// order, so floating-point sums are bit-identical from 1 to N threads.
//
// Thread count: CHEMLARP_THREADS in the environment, else set_thread_count(), else the hardware concurrency.
// A nested call (from inside a body) runs serially on the calling thread. The first exception thrown by a
// body is rethrown by parallel_for after all indices have finished. Builds without threads (WebAssembly
// without pthreads, or CHEMLARP_NO_THREADS) run everything serially with the same results.

#include <cstddef>
#include <functional>
#include <vector>

namespace chemlarp {

// Number of threads parallel_for uses (at least 1). Changing it waits for running work to finish.
void set_thread_count(int n);
int thread_count();

void parallel_for(std::size_t count, const std::function<void(std::size_t)>& body);

// Splits [0, n) into `parts` contiguous ranges of nearly equal size; range k is [begin(k), begin(k+1)).
inline std::size_t partition_begin(std::size_t n, std::size_t parts, std::size_t k) { return parts ? n * k / parts : 0; }

// Deterministic sum of per-partition vectors: body(k, partial) fills `partial` (zeroed, length `size`) for
// partition k; the partials are added in partition order into the returned vector.
std::vector<double> parallel_sum(std::size_t parts, std::size_t size,
                                 const std::function<void(std::size_t, std::vector<double>&)>& body);

}  // namespace chemlarp
