#include "chemlarp/parallel.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <exception>

#if defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)
#define CHEMLARP_NO_THREADS 1
#endif

#ifndef CHEMLARP_NO_THREADS
#include <condition_variable>
#include <mutex>
#include <thread>
#endif

namespace chemlarp {

#ifdef CHEMLARP_NO_THREADS

void set_thread_count(int) {}
int thread_count() { return 1; }
void parallel_for(std::size_t count, const std::function<void(std::size_t)>& body) {
  for (std::size_t i = 0; i < count; ++i) body(i);
}

#else

namespace {

thread_local bool tl_inside = false;

class Pool {
 public:
  static Pool& get() {
    static Pool p;
    return p;
  }
  ~Pool() { stop(); }

  int size() {
    std::lock_guard<std::mutex> g(run_);
    return threads_;
  }

  void resize(int n) {
    std::lock_guard<std::mutex> g(run_);
    n = std::max(1, n);
    if (n == threads_) return;
    stop();
    start(n);
  }

  void run(std::size_t count, const std::function<void(std::size_t)>& body) {
    if (count == 0) return;
    if (tl_inside || count == 1) {
      for (std::size_t i = 0; i < count; ++i) body(i);
      return;
    }
    std::lock_guard<std::mutex> g(run_);  // one parallel region at a time
    if (threads_ <= 1) {
      tl_inside = true;
      try {
        for (std::size_t i = 0; i < count; ++i) body(i);
      } catch (...) {
        tl_inside = false;
        throw;
      }
      tl_inside = false;
      return;
    }
    {
      std::lock_guard<std::mutex> l(m_);
      body_ = &body;
      count_ = count;
      next_.store(0);
      error_ = nullptr;
      active_ = static_cast<int>(workers_.size());
      ++generation_;
    }
    cv_.notify_all();
    work();  // the caller takes part
    std::unique_lock<std::mutex> l(m_);
    done_cv_.wait(l, [&] { return active_ == 0; });
    body_ = nullptr;
    if (error_) std::rethrow_exception(error_);
  }

 private:
  Pool() {
    int n = static_cast<int>(std::thread::hardware_concurrency());
    if (const char* e = std::getenv("CHEMLARP_THREADS")) n = std::atoi(e);
    start(std::max(1, std::min(n, 256)));
  }

  void start(int n) {
    threads_ = n;
    quit_ = false;
    for (int k = 1; k < n; ++k) workers_.emplace_back([this] { loop(); });
  }

  void stop() {
    {
      std::lock_guard<std::mutex> l(m_);
      quit_ = true;
    }
    cv_.notify_all();
    for (auto& t : workers_) t.join();
    workers_.clear();
    threads_ = 1;
  }

  void work() {
    tl_inside = true;
    for (;;) {
      const std::size_t i = next_.fetch_add(1);
      if (i >= count_) break;
      try {
        (*body_)(i);
      } catch (...) {
        std::lock_guard<std::mutex> l(m_);
        if (!error_) error_ = std::current_exception();
      }
    }
    tl_inside = false;
  }

  void loop() {
    std::uint64_t seen = 0;
    for (;;) {
      {
        std::unique_lock<std::mutex> l(m_);
        cv_.wait(l, [&] { return quit_ || generation_ != seen; });
        if (quit_) return;
        seen = generation_;
      }
      work();
      {
        std::lock_guard<std::mutex> l(m_);
        if (--active_ == 0) done_cv_.notify_one();
      }
    }
  }

  std::mutex run_;  // serialises parallel regions and resizes
  std::mutex m_;
  std::condition_variable cv_, done_cv_;
  std::vector<std::thread> workers_;
  const std::function<void(std::size_t)>* body_ = nullptr;
  std::size_t count_ = 0;
  std::atomic<std::size_t> next_{0};
  std::exception_ptr error_;
  int active_ = 0, threads_ = 1;
  std::uint64_t generation_ = 0;
  bool quit_ = false;
};

}  // namespace

void set_thread_count(int n) { Pool::get().resize(n); }
int thread_count() { return Pool::get().size(); }
void parallel_for(std::size_t count, const std::function<void(std::size_t)>& body) { Pool::get().run(count, body); }

#endif

std::vector<double> parallel_sum(std::size_t parts, std::size_t size,
                                 const std::function<void(std::size_t, std::vector<double>&)>& body) {
  std::vector<std::vector<double>> partial(parts, std::vector<double>(size, 0.0));
  parallel_for(parts, [&](std::size_t k) { body(k, partial[k]); });
  std::vector<double> out(size, 0.0);
  for (std::size_t k = 0; k < parts; ++k)
    for (std::size_t i = 0; i < size; ++i) out[i] += partial[k][i];
  return out;
}

}  // namespace chemlarp
