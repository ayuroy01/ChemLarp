#include "chemlarp/sparse.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace chemlarp {

SparseMatrix SparseMatrix::from_dense(const double* a, std::size_t n, double drop) {
  SparseMatrix m;
  m.n = n;
  m.colptr.assign(n + 1, 0);
  for (std::size_t c = 0; c < n; ++c) {
    for (std::size_t r = 0; r < n; ++r) {
      const double v = a[r * n + c];
      if (std::fabs(v) > drop || r == c) {  // keep the diagonal in the pattern
        m.rowind.push_back(r);
        m.values.push_back(v);
      }
    }
    m.colptr[c + 1] = m.rowind.size();
  }
  return m;
}

std::vector<std::size_t> minimum_degree_order(const SparseMatrix& a) {
  const std::size_t n = a.n;
  // Elimination graph of A + A^T (no self loops).
  std::vector<std::set<std::size_t>> adj(n);
  for (std::size_t c = 0; c < n; ++c)
    for (std::size_t p = a.colptr[c]; p < a.colptr[c + 1]; ++p) {
      const std::size_t r = a.rowind[p];
      if (r == c) continue;
      adj[r].insert(c);
      adj[c].insert(r);
    }
  std::vector<char> done(n, 0);
  // Buckets by degree; ties broken by the lowest index, so the order is deterministic.
  std::set<std::pair<std::size_t, std::size_t>> queue;
  for (std::size_t v = 0; v < n; ++v) queue.insert({adj[v].size(), v});
  std::vector<std::size_t> order;
  order.reserve(n);
  while (!queue.empty()) {
    const std::size_t v = queue.begin()->second;
    queue.erase(queue.begin());
    done[v] = 1;
    order.push_back(v);
    const std::vector<std::size_t> nb(adj[v].begin(), adj[v].end());
    for (std::size_t u : nb) {
      queue.erase({adj[u].size(), u});
      adj[u].erase(v);
    }
    // The neighbours become a clique.
    for (std::size_t i = 0; i < nb.size(); ++i)
      for (std::size_t j = i + 1; j < nb.size(); ++j) {
        adj[nb[i]].insert(nb[j]);
        adj[nb[j]].insert(nb[i]);
      }
    for (std::size_t u : nb) queue.insert({adj[u].size(), u});
    adj[v].clear();
  }
  return order;
}

bool SparseLu::factor(const SparseMatrix& a, std::vector<std::size_t> q, double tol) {
  const std::size_t n = a.n;
  n_ = n;
  q_ = q.empty() ? minimum_degree_order(a) : std::move(q);
  constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();
  pinv_.assign(n, kNone);
  lp_.assign(n + 1, 0);
  up_.assign(n + 1, 0);
  li_.clear(), lx_.clear(), ui_.clear(), ux_.clear();
  std::vector<double> x(n, 0.0);
  std::vector<std::size_t> xi(2 * n), stack(n), mark(n, kNone), pstack(n);
  double scale = 0.0;
  for (double v : a.values) scale = std::max(scale, std::fabs(v));
  if (n == 0) return true;
  if (scale == 0.0 || !std::isfinite(scale)) return false;
  for (std::size_t k = 0; k < n; ++k) {
    lp_[k] = li_.size();
    up_[k] = ui_.size();
    const std::size_t col = q_[k];
    // Reach of A(:,col) in the graph of L (columns of rows already pivoted): DFS, topological order.
    std::size_t top = n;
    for (std::size_t p = a.colptr[col]; p < a.colptr[col + 1]; ++p) {
      const std::size_t r0 = a.rowind[p];
      if (mark[r0] == k) continue;
      std::size_t head = 0;
      stack[0] = r0;
      while (true) {
        const std::size_t j = stack[head];
        const std::size_t jn = pinv_[j];  // column of L for row j, if pivoted
        if (mark[j] != k) {
          mark[j] = k;
          pstack[head] = jn == kNone ? 0 : lp_[jn] + 1;  // skip the unit diagonal stored first
        }
        bool pushed = false;
        if (jn != kNone) {
          const std::size_t end = lp_[jn + 1];
          for (std::size_t p2 = pstack[head]; p2 < end; ++p2) {
            const std::size_t i = li_[p2];
            if (mark[i] == k) continue;
            pstack[head] = p2 + 1;
            stack[++head] = i;
            pushed = true;
            break;
          }
        }
        if (pushed) continue;
        xi[--top] = j;
        if (head == 0) break;
        --head;
      }
    }
    // x = L \\ A(:,col) on the reach.
    for (std::size_t p = top; p < n; ++p) x[xi[p]] = 0.0;
    for (std::size_t p = a.colptr[col]; p < a.colptr[col + 1]; ++p) x[a.rowind[p]] = a.values[p];
    for (std::size_t p = top; p < n; ++p) {
      const std::size_t j = xi[p];
      const std::size_t jn = pinv_[j];
      if (jn == kNone) continue;
      const double xj = x[j];
      for (std::size_t p2 = lp_[jn] + 1; p2 < lp_[jn + 1]; ++p2) x[li_[p2]] -= lx_[p2] * xj;
    }
    // Pivot: the diagonal row (col) when large enough, else the largest unpivoted entry.
    std::size_t ipiv = kNone;
    double amax = -1.0;
    for (std::size_t p = top; p < n; ++p) {
      const std::size_t i = xi[p];
      if (pinv_[i] == kNone) {
        const double t = std::fabs(x[i]);
        if (t > amax) amax = t, ipiv = i;
      } else {
        ui_.push_back(pinv_[i]);
        ux_.push_back(x[i]);
      }
    }
    if (ipiv == kNone || amax <= 1e-300 + 1e-15 * scale) return false;
    if (pinv_[col] == kNone && mark[col] == k && std::fabs(x[col]) >= tol * amax) ipiv = col;
    const double pivot = x[ipiv];
    ui_.push_back(k);
    ux_.push_back(pivot);
    pinv_[ipiv] = k;
    li_.push_back(ipiv);
    lx_.push_back(1.0);
    for (std::size_t p = top; p < n; ++p) {
      const std::size_t i = xi[p];
      if (pinv_[i] == kNone) {
        li_.push_back(i);
        lx_.push_back(x[i] / pivot);
      }
      x[i] = 0.0;
    }
  }
  lp_[n] = li_.size();
  up_[n] = ui_.size();
  // Row indices of L in pivot order.
  for (auto& i : li_) i = pinv_[i];
  return true;
}

void SparseLu::solve(std::vector<double>& b) const {
  const std::size_t n = n_;
  std::vector<double> x(n);
  for (std::size_t i = 0; i < n; ++i) x[pinv_[i]] = b[i];  // x = P b
  for (std::size_t j = 0; j < n; ++j) {                      // L y = x (unit diagonal first)
    const double xj = x[j];
    for (std::size_t p = lp_[j] + 1; p < lp_[j + 1]; ++p) x[li_[p]] -= lx_[p] * xj;
  }
  for (std::size_t j = n; j-- > 0;) {  // U z = y (diagonal last)
    x[j] /= ux_[up_[j + 1] - 1];
    const double xj = x[j];
    for (std::size_t p = up_[j]; p < up_[j + 1] - 1; ++p) x[ui_[p]] -= ux_[p] * xj;
  }
  for (std::size_t k = 0; k < n; ++k) b[q_[k]] = x[k];  // b = Q z
}

}  // namespace chemlarp
