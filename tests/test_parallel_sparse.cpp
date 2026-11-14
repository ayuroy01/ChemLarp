// Thread pool determinism and the sparse LU (against the dense LU and in the BDF integrator).
#include <atomic>
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "chemlarp/chemlarp_c.h"
#include "chemlarp/ode.hpp"
#include "chemlarp/parallel.hpp"
#include "chemlarp/sparse.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

TEST(parallel_for_covers_every_index_and_rethrows) {
  for (int threads : {1, 3, 8}) {
    set_thread_count(threads);
    CHECK_EQ(thread_count(), threads);
    std::vector<int> hit(1000, 0);
    parallel_for(hit.size(), [&](std::size_t i) { hit[i] += 1; });
    for (int h : hit) CHECK_EQ(h, 1);
    // Nested regions run serially inside the body.
    std::atomic<int> inner{0};
    parallel_for(10, [&](std::size_t) { parallel_for(10, [&](std::size_t) { ++inner; }); });
    CHECK_EQ(inner.load(), 100);
    bool thrown = false;
    try {
      parallel_for(50, [&](std::size_t i) {
        if (i == 17) throw std::runtime_error("x");
      });
    } catch (const std::runtime_error&) {
      thrown = true;
    }
    CHECK(thrown);
  }
  // parallel_sum adds the partitions in order: the same bits for any thread count.
  auto run = [] {
    return parallel_sum(16, 3, [](std::size_t k, std::vector<double>& p) {
      for (std::size_t i = 0; i < 1000; ++i) p[i % 3] += 1.0 / static_cast<double>(k * 1000 + i + 1);
    });
  };
  set_thread_count(1);
  const auto a = run();
  set_thread_count(7);
  const auto b = run();
  CHECK(a == b);
  set_thread_count(4);
}

TEST(qc_results_do_not_depend_on_thread_count) {
  std::string ref;
  for (int threads : {1, 2, 5}) {
    set_thread_count(threads);
    char* s = cl_qc_json("CC(=O)O", "b3lyp", "3-21G", 0, 1, 0);
    const std::string j(s);
    cl_free_string(s);
    CHECK(j.rfind("{\"ok\":true", 0) == 0);
    if (ref.empty()) ref = j;
    CHECK(j == ref);  // the whole JSON: energy, orbitals, charges, dipole
  }
  set_thread_count(4);
}

TEST(sparse_lu_matches_dense) {
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  for (std::size_t n : {1, 5, 40, 200}) {
    for (int trial = 0; trial < 5; ++trial) {
      // Random pattern (about 4 entries per column), a modest diagonal so that pivoting happens.
      std::vector<double> a(n * n, 0.0);
      for (std::size_t c = 0; c < n; ++c) {
        a[c * n + c] = (trial % 2 ? 0.01 : 2.0) * u(rng);
        for (int k = 0; k < 4; ++k) a[(rng() % n) * n + c] = u(rng);
      }
      std::vector<double> b(n);
      for (auto& x : b) x = u(rng);
      LuFactor dense;
      if (!dense.factor(a, n)) continue;  // singular draw
      std::vector<double> xd = b;
      dense.solve(xd);
      SparseLu sp;
      CHECK(sp.factor(SparseMatrix::from_dense(a.data(), n)));
      std::vector<double> xs = b;
      sp.solve(xs);
      double res = 0.0, scale = 0.0, diff = 0.0;
      for (std::size_t r = 0; r < n; ++r) {
        double ax = 0.0;
        for (std::size_t c = 0; c < n; ++c) ax += a[r * n + c] * xs[c];
        res = std::max(res, std::fabs(ax - b[r]));
        scale = std::max(scale, std::fabs(xd[r]));
        diff = std::max(diff, std::fabs(xs[r] - xd[r]));
      }
      CHECK(res < 1e-9);
      CHECK(diff < 1e-7 * std::max(1.0, scale));
    }
  }
  // Singular matrix: reported.
  std::vector<double> z = {1, 2, 2, 4};
  SparseLu sp;
  CHECK(!sp.factor(SparseMatrix::from_dense(z.data(), 2)));
  // A tridiagonal matrix in natural order has no fill: L + U has 3n - 2 + n entries (unit diagonal of L).
  const std::size_t n = 100;
  std::vector<double> t(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    t[i * n + i] = 4.0;
    if (i) t[i * n + i - 1] = t[(i - 1) * n + i] = -1.0;
  }
  CHECK(sp.factor(SparseMatrix::from_dense(t.data(), n)));
  CHECK_EQ(sp.nnz_l() + sp.nnz_u(), 4 * n - 2);
}

TEST(bdf_sparse_and_dense_agree_on_brusselator) {
  // 1-D Brusselator with diffusion (Hairer & Wanner, Solving ODEs II, problem BRUS), N = 60 grid points:
  // 120 stiff equations with a banded Jacobian.
  const std::size_t N = 60, n = 2 * N;
  const double alpha = 0.02 * (N + 1) * (N + 1);
  OdeRhs f = [&](double, const std::vector<double>& y, std::vector<double>& d) {
    d.assign(n, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
      const double u = y[2 * i], v = y[2 * i + 1];
      const double ul = i ? y[2 * i - 2] : 1.0, ur = i + 1 < N ? y[2 * i + 2] : 1.0;
      const double vl = i ? y[2 * i - 1] : 3.0, vr = i + 1 < N ? y[2 * i + 3] : 3.0;
      d[2 * i] = 1.0 + u * u * v - 4.0 * u + alpha * (ul - 2.0 * u + ur);
      d[2 * i + 1] = 3.0 * u - u * u * v + alpha * (vl - 2.0 * v + vr);
    }
  };
  OdeJacobian jac = [&](double, const std::vector<double>& y, std::vector<double>& J) {
    J.assign(n * n, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
      const double u = y[2 * i], v = y[2 * i + 1];
      const std::size_t a = 2 * i, b = 2 * i + 1;
      J[a * n + a] = 2.0 * u * v - 4.0 - 2.0 * alpha;
      J[a * n + b] = u * u;
      J[b * n + a] = 3.0 - 2.0 * u * v;
      J[b * n + b] = -u * u - 2.0 * alpha;
      if (i) J[a * n + a - 2] = J[b * n + b - 2] = alpha;
      if (i + 1 < N) J[a * n + a + 2] = J[b * n + b + 2] = alpha;
    }
  };
  std::vector<std::vector<double>> finals;
  for (int solver : {1, 2}) {
    std::vector<double> y(n);
    for (std::size_t i = 0; i < N; ++i) {
      const double x = static_cast<double>(i + 1) / (N + 1);
      y[2 * i] = 1.0 + std::sin(2.0 * 3.14159265358979323846 * x);
      y[2 * i + 1] = 3.0;
    }
    BdfOptions o;
    o.rtol = 1e-8;
    o.atol = {1e-10};
    o.linear_solver = solver;
    BdfStats st;
    Diagnostics d;
    CHECK(integrate_bdf(f, &jac, 0.0, y, 10.0, {}, [](double, const std::vector<double>&) {}, nullptr, o, st, d));
    CHECK_EQ(st.sparse, solver == 2);
    if (solver == 2) CHECK(st.factor_nnz < n * n / 10);
    finals.push_back(y);
  }
  double diff = 0.0;
  for (std::size_t i = 0; i < n; ++i) diff = std::max(diff, std::fabs(finals[0][i] - finals[1][i]));
  CHECK(diff < 1e-6);
}
