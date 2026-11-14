#pragma once
// Sparse LU factorisation for unsymmetric systems (A x = b).
//
// Column ordering: minimum degree on the pattern of A + A^T (Tinney–Walker scheme 2, with exact degrees
// on the elimination graph), which keeps fill low for the nearly symmetric patterns of reaction-network
// and discretised-PDE Jacobians. Factorisation: left-looking Gilbert–Peierls LU (Gilbert & Peierls, SIAM J.
// Sci. Stat. Comput. 9 (1988) 862; the formulation of Davis, "Direct Methods for Sparse Linear Systems",
// SIAM 2006, cs_lu) with threshold partial pivoting: the diagonal entry is kept when it is at least `tol`
// times the largest entry of its column, else the largest is taken. Result: P A Q = L U.

#include <cstddef>
#include <vector>

namespace chemlarp {

// Compressed sparse column matrix.
struct SparseMatrix {
  std::size_t n = 0;               // square
  std::vector<std::size_t> colptr;  // n + 1
  std::vector<std::size_t> rowind;  // nnz
  std::vector<double> values;       // nnz

  // From a dense row-major n x n array; entries with |a| <= drop are left out.
  static SparseMatrix from_dense(const double* a, std::size_t n, double drop = 0.0);
  std::size_t nnz() const { return rowind.size(); }
};

// Fill-reducing column permutation (q[k] = column eliminated k-th).
std::vector<std::size_t> minimum_degree_order(const SparseMatrix& a);

class SparseLu {
 public:
  // Factorises with ordering `q` (empty: computed by minimum_degree_order). Returns false if singular.
  bool factor(const SparseMatrix& a, std::vector<std::size_t> q = {}, double tol = 0.1);
  void solve(std::vector<double>& b) const;  // in place
  std::size_t nnz_l() const { return li_.size(); }
  std::size_t nnz_u() const { return ui_.size(); }
  const std::vector<std::size_t>& column_order() const { return q_; }

 private:
  std::size_t n_ = 0;
  std::vector<std::size_t> lp_, li_, up_, ui_;  // L unit lower (CSC, diagonal first), U upper (CSC, diagonal last)
  std::vector<double> lx_, ux_;
  std::vector<std::size_t> pinv_, q_;           // row permutation (inverse) and column order
};

}  // namespace chemlarp
