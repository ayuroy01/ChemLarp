#include "chemlarp/column.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <string>

namespace chemlarp {
namespace {

std::string num(double v) {
  char b[32];
  std::snprintf(b, sizeof b, "%g", v);
  return b;
}

// Dense LU with partial pivoting for the small diagonal blocks (row-major n × n).
bool lu_factor(std::vector<double>& a, std::vector<std::size_t>& piv, std::size_t n) {
  piv.resize(n);
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t p = k;
    for (std::size_t r = k + 1; r < n; ++r)
      if (std::fabs(a[r * n + k]) > std::fabs(a[p * n + k])) p = r;
    if (!(std::fabs(a[p * n + k]) > 1e-300)) return false;
    piv[k] = p;
    if (p != k)
      for (std::size_t c = 0; c < n; ++c) std::swap(a[k * n + c], a[p * n + c]);
    for (std::size_t r = k + 1; r < n; ++r) {
      const double f = a[r * n + k] /= a[k * n + k];
      if (f != 0.0)
        for (std::size_t c = k + 1; c < n; ++c) a[r * n + c] -= f * a[k * n + c];
    }
  }
  return true;
}

// Solves for `cols` right-hand sides stored column-wise in b (n × cols, row-major).
void lu_solve(const std::vector<double>& a, const std::vector<std::size_t>& piv, std::size_t n, std::vector<double>& b,
              std::size_t cols) {
  for (std::size_t k = 0; k < n; ++k)
    if (piv[k] != k)
      for (std::size_t c = 0; c < cols; ++c) std::swap(b[k * cols + c], b[piv[k] * cols + c]);
  for (std::size_t c = 0; c < cols; ++c) {
    for (std::size_t r = 1; r < n; ++r) {
      double s = b[r * cols + c];
      for (std::size_t k = 0; k < r; ++k) s -= a[r * n + k] * b[k * cols + c];
      b[r * cols + c] = s;
    }
    for (std::size_t r = n; r-- > 0;) {
      double s = b[r * cols + c];
      for (std::size_t k = r + 1; k < n; ++k) s -= a[r * n + k] * b[k * cols + c];
      b[r * cols + c] = s / a[r * n + r];
    }
  }
}

struct StageProps {
  std::vector<double> K, hL, hV;
};

struct Column {
  const ProcessComponents& c;
  const ColumnSpec& s;
  std::size_t N = 0, C = 0, m = 0;
  bool total = true;
  std::vector<double> P;                // stage pressures
  std::vector<std::vector<double>> f;   // feed component flows per stage
  std::vector<double> Hf;               // feed enthalpy per stage
  double F = 0.0, Fs = 1.0, Hs = 1.0, D = 0.0, B = 0.0;
  Diagnostics notes;                    // property warnings, merged once at the end

  Column(const ProcessComponents& comps, const ColumnSpec& spec) : c(comps), s(spec) {}

  double T(const std::vector<double>& X, std::size_t j) const { return X[j * m]; }
  const double* l(const std::vector<double>& X, std::size_t j) const { return &X[j * m + 1]; }
  const double* v(const std::vector<double>& X, std::size_t j) const { return &X[j * m + 1 + C]; }

  std::optional<StageProps> props(const std::vector<double>& X, std::size_t j) {
    const double* lj = l(X, j);
    double L = 0.0;
    for (std::size_t i = 0; i < C; ++i) L += std::max(lj[i], 0.0);
    std::vector<double> x(C);
    for (std::size_t i = 0; i < C; ++i) x[i] = L > 0.0 ? std::max(lj[i], 0.0) / L : 1.0 / double(C);
    StageProps p;
    auto K = k_values(c, T(X, j), P[j], x, notes);
    auto hL = component_enthalpies(c, T(X, j), P[j], false, notes);
    auto hV = component_enthalpies(c, T(X, j), P[j], true, notes);
    if (!K || !hL || !hV) return std::nullopt;
    p.K = std::move(*K);
    p.hL = std::move(*hL);
    p.hV = std::move(*hV);
    return p;
  }

  // Residuals of stage j (m values) from the variables and the properties of stages j−1, j, j+1.
  void residual(const std::vector<double>& X, std::size_t j, const std::vector<StageProps>& p, double* r) const {
    const double* lj = l(X, j);
    const double* vj = v(X, j);
    double L = 0.0, V = 0.0;
    for (std::size_t i = 0; i < C; ++i) L += lj[i], V += vj[i];
    const double Lg = std::max(L, 1e-12 * Fs), Vg = std::max(V, 1e-12 * Fs);
    for (std::size_t i = 0; i < C; ++i) {
      double in = f[j][i];
      if (j > 0) in += l(X, j - 1)[i];
      if (j + 1 < N) in += v(X, j + 1)[i];
      r[1 + i] = (lj[i] + vj[i] - in) / Fs;
    }
    const StageProps& pj = p[j];
    if (j == 0 && total) {
      // vj holds the distillate: same composition as the reflux, at its bubble point; L = R D.
      double sk = 0.0;
      for (std::size_t i = 0; i < C; ++i) sk += pj.K[i] * lj[i] / Lg;
      for (std::size_t i = 0; i + 1 < C; ++i) r[1 + C + i] = vj[i] / Vg - lj[i] / Lg;
      r[1 + C + (C - 1)] = sk - 1.0;
      r[0] = (L - s.reflux_ratio * V) / Fs;
      return;
    }
    for (std::size_t i = 0; i < C; ++i) r[1 + C + i] = vj[i] / Vg - pj.K[i] * lj[i] / Lg;
    if (j == 0) {
      r[0] = (L - s.reflux_ratio * V) / Fs;  // partial condenser: reflux ratio
    } else if (j + 1 == N) {
      r[0] = (L - B) / Fs;  // reboiler: bottoms flow (the distillate flow by the overall balance)
    } else {
      double H = -Hf[j];
      for (std::size_t i = 0; i < C; ++i) {
        H += lj[i] * pj.hL[i] + vj[i] * pj.hV[i];
        H -= l(X, j - 1)[i] * p[j - 1].hL[i] + v(X, j + 1)[i] * p[j + 1].hV[i];
      }
      r[0] = H / Hs;
    }
  }

  bool all_props(const std::vector<double>& X, std::vector<StageProps>& p) {
    p.resize(N);
    for (std::size_t j = 0; j < N; ++j) {
      auto q = props(X, j);
      if (!q) return false;
      p[j] = std::move(*q);
    }
    return true;
  }

  bool residuals(const std::vector<double>& X, std::vector<double>& R, std::vector<StageProps>& p) {
    if (!all_props(X, p)) return false;
    R.assign(N * m, 0.0);
    for (std::size_t j = 0; j < N; ++j) residual(X, j, p, &R[j * m]);
    return true;
  }

  // Newton step from the block-tridiagonal Jacobian (finite differences).
  bool newton_step(std::vector<double>& X, const std::vector<double>& R, std::vector<StageProps>& p, std::vector<double>& dx) {
    const std::size_t mm = m * m;
    std::vector<double> A(N * mm, 0.0), Bk(N * mm, 0.0), Ck(N * mm, 0.0);
    std::vector<double> r(m);
    for (std::size_t k = 0; k < N; ++k) {
      const StageProps saved = p[k];
      for (std::size_t q = 0; q < m; ++q) {
        const std::size_t idx = k * m + q;
        const double x0 = X[idx];
        const double h = q == 0 ? 1e-6 * std::max(x0, 1.0) : 1e-7 * std::max(std::fabs(x0), 1e-4 * Fs);
        X[idx] = x0 + h;
        auto pk = props(X, k);
        if (!pk) return false;
        p[k] = std::move(*pk);
        for (std::size_t j = k == 0 ? 0 : k - 1; j <= std::min(k + 1, N - 1); ++j) {
          residual(X, j, p, r.data());
          double* blk = j + 1 == k ? &Ck[j * mm] : j == k ? &Bk[j * mm] : &A[j * mm];
          for (std::size_t e = 0; e < m; ++e) blk[e * m + q] = (r[e] - R[j * m + e]) / h;
        }
        X[idx] = x0;
      }
      p[k] = saved;
    }
    // Block Thomas: Γ_j = M_j⁻¹ C_j, ψ_j = M_j⁻¹ (−R_j − A_j ψ_{j−1}), M_j = B_j − A_j Γ_{j−1}.
    std::vector<double> G(N * mm, 0.0), psi(N * m, 0.0), M(mm);
    std::vector<std::size_t> piv;
    for (std::size_t j = 0; j < N; ++j) {
      std::copy(&Bk[j * mm], &Bk[j * mm] + mm, M.begin());
      std::vector<double> rhs(m);
      for (std::size_t e = 0; e < m; ++e) rhs[e] = -R[j * m + e];
      if (j > 0) {
        const double* a = &A[j * mm];
        const double* g = &G[(j - 1) * mm];
        for (std::size_t r1 = 0; r1 < m; ++r1)
          for (std::size_t k = 0; k < m; ++k) {
            const double av = a[r1 * m + k];
            if (av == 0.0) continue;
            for (std::size_t c1 = 0; c1 < m; ++c1) M[r1 * m + c1] -= av * g[k * m + c1];
            rhs[r1] -= av * psi[(j - 1) * m + k];
          }
      }
      if (!lu_factor(M, piv, m)) return false;
      lu_solve(M, piv, m, rhs, 1);
      std::copy(rhs.begin(), rhs.end(), &psi[j * m]);
      if (j + 1 < N) {
        std::vector<double> g(&Ck[j * mm], &Ck[j * mm] + mm);
        lu_solve(M, piv, m, g, m);
        std::copy(g.begin(), g.end(), &G[j * mm]);
      }
    }
    dx.assign(N * m, 0.0);
    for (std::size_t j = N; j-- > 0;) {
      for (std::size_t e = 0; e < m; ++e) {
        double v1 = psi[j * m + e];
        if (j + 1 < N)
          for (std::size_t k = 0; k < m; ++k) v1 -= G[j * mm + e * m + k] * dx[(j + 1) * m + k];
        dx[j * m + e] = v1;
      }
    }
    return true;
  }
};

double max_abs(const std::vector<double>& v) {
  double r = 0.0;
  for (double x : v) r = std::max(r, std::fabs(x));
  return std::isfinite(r) ? r : 1e300;
}

}  // namespace

std::optional<ColumnResult> distillation_column(const ProcessComponents& c, const ColumnSpec& s, Diagnostics& diag) {
  auto bad = [&](const std::string& msg) {
    diag.error(DiagCode::InvalidMixtureInput, msg);
    return std::nullopt;
  };
  Column col(c, s);
  col.N = s.stages;
  col.C = c.size();
  col.m = 2 * col.C + 1;
  col.total = s.condenser == Condenser::Total;
  const std::size_t N = col.N, C = col.C, m = col.m;
  if (N < 3) return bad("a column needs at least 3 stages (condenser, one stage, reboiler)");
  if (c.condensable.size() != C) return bad("every component of a column needs a vapour-pressure equation");
  if (!(s.reflux_ratio > 0.0) || !std::isfinite(s.reflux_ratio)) return bad("the reflux ratio must be positive");
  if (!(s.P_top_Pa > 0.0)) return bad("the column needs a top pressure");
  if ((s.distillate_rate > 0.0) == (s.bottoms_rate > 0.0)) return bad("give either the distillate or the bottoms flow");
  if (s.feeds.empty()) return bad("the column has no feed");
  const double Pb = s.P_bottom_Pa > 0.0 ? s.P_bottom_Pa : s.P_top_Pa;
  for (std::size_t j = 0; j < N; ++j) col.P.push_back(s.P_top_Pa + (Pb - s.P_top_Pa) * double(j) / double(N - 1));
  col.f.assign(N, std::vector<double>(C, 0.0));
  col.Hf.assign(N, 0.0);
  std::vector<double> fz(C, 0.0), q(N, 0.0), Fj(N, 0.0);
  double Tf = 0.0;
  for (const auto& fd : s.feeds) {
    if (fd.stage < 1 || fd.stage > N) return bad("feed stage " + std::to_string(fd.stage) + " is outside the column");
    if (fd.stream.n.size() != C) return bad("a feed needs one flow per component");
    const std::size_t j = fd.stage - 1;
    const double Ff = fd.stream.total();
    for (std::size_t i = 0; i < C; ++i) col.f[j][i] += fd.stream.n[i], fz[i] += fd.stream.n[i];
    col.Hf[j] += fd.stream.H_W;
    q[j] += (1.0 - fd.stream.vapor_fraction) * Ff;
    Fj[j] += Ff;
    col.F += Ff;
    Tf += fd.stream.T_K * Ff;
  }
  if (!(col.F > 0.0)) return bad("the feed has no flow");
  Tf /= col.F;
  col.D = s.distillate_rate > 0.0 ? s.distillate_rate : col.F - s.bottoms_rate;
  col.B = col.F - col.D;
  if (!(col.D > 0.0) || !(col.B > 0.0)) return bad("the distillate and bottoms flows must both lie between 0 and the feed flow");
  col.Fs = col.F;
  col.Hs = col.F * 1e4;  // 10 kJ/mol

  // ---- Starting profile ----
  // Product split by volatility at the feed temperature: the lightest components fill the distillate.
  std::vector<std::size_t> order(C);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return c.vle.antoine[a].pressure(Tf) > c.vle.antoine[b].pressure(Tf);
  });
  std::vector<double> xd(C), xb(C);
  double rest = col.D;
  for (std::size_t i : order) {
    const double d = std::min(fz[i], rest);
    rest -= d;
    xd[i] = d + 1e-6 * col.F;
    xb[i] = fz[i] - d + 1e-6 * col.F;
  }
  auto normalize = [](std::vector<double>& x) {
    double t = 0.0;
    for (double v : x) t += v;
    for (double& v : x) v /= t;
  };
  normalize(xd);
  normalize(xb);
  Diagnostics& notes = col.notes;
  double Ttop = Tf, Tbot = Tf;
  if (auto t = col.total ? bubble_temperature(c.vle, col.P[0], xd, notes) : dew_temperature(c.vle, col.P[0], xd, notes)) Ttop = t->T_K;
  if (auto t = bubble_temperature(c.vle, col.P[N - 1], xb, notes)) Tbot = t->T_K;
  if (Tbot < Ttop) std::swap(Tbot, Ttop);
  std::vector<double> Ts(N), L(N), V(N);
  std::vector<std::vector<double>> x(N, std::vector<double>(C));
  for (std::size_t j = 0; j < N; ++j) {
    const double w = double(j) / double(N - 1);
    Ts[j] = Ttop + (Tbot - Ttop) * w;
    for (std::size_t i = 0; i < C; ++i) x[j][i] = xd[i] + (xb[i] - xd[i]) * w;
  }
  // Constant molar overflow: V[0] is the distillate (total condenser: liquid, carried in the vapour slot).
  L[0] = s.reflux_ratio * col.D;
  V[0] = col.D;
  double Vin = L[0] + V[0] - Fj[0];  // vapour from stage 2
  for (std::size_t j = 1; j + 1 < N; ++j) {
    L[j] = L[j - 1] + q[j];
    V[j] = Vin;
    Vin = V[j] - (Fj[j] - q[j]);
  }
  L[N - 1] = col.B;
  V[N - 1] = std::max(L[N - 2] + Fj[N - 1] - col.B, 1e-3 * col.F);
  for (std::size_t j = 1; j + 1 < N; ++j) {
    L[j] = std::max(L[j], 1e-3 * col.F);
    V[j] = std::max(V[j], 1e-3 * col.F);
  }
  // Bubble-point (Wang–Henke) passes with these flows: tridiagonal component balances, then stage
  // temperatures from the bubble point of each stage's liquid.
  for (int pass = 0; pass < 12; ++pass) {
    std::vector<std::vector<double>> K(N);
    for (std::size_t j = 0; j < N; ++j) {
      auto k = k_values(c, Ts[j], col.P[j], x[j], notes);
      if (!k) {
        diag.error(DiagCode::ColumnNotConverged, "the starting profile could not be evaluated");
        return std::nullopt;
      }
      K[j] = *k;
    }
    for (std::size_t i = 0; i < C; ++i) {
      // a_j x_{j−1} + b_j x_j + c_j x_{j+1} = f_ij
      std::vector<double> a(N, 0.0), b(N), cc(N, 0.0), d(N);
      for (std::size_t j = 0; j < N; ++j) {
        const bool tc = j == 0 && col.total;
        b[j] = L[j] + (tc ? V[0] : V[j] * K[j][i]);
        if (j > 0) a[j] = -L[j - 1];
        if (j + 1 < N) cc[j] = -V[j + 1] * K[j + 1][i];
        d[j] = col.f[j][i];
      }
      for (std::size_t j = 1; j < N; ++j) {
        const double w = a[j] / b[j - 1];
        b[j] -= w * cc[j - 1];
        d[j] -= w * d[j - 1];
      }
      std::vector<double> xi(N);
      xi[N - 1] = d[N - 1] / b[N - 1];
      for (std::size_t j = N - 1; j-- > 0;) xi[j] = (d[j] - cc[j] * xi[j + 1]) / b[j];
      for (std::size_t j = 0; j < N; ++j) x[j][i] = std::max(xi[j], 1e-12);
    }
    double change = 0.0;
    for (std::size_t j = 0; j < N; ++j) {
      normalize(x[j]);
      if (auto t = bubble_temperature(c.vle, col.P[j], x[j], notes)) {
        change = std::max(change, std::fabs(t->T_K - Ts[j]));
        Ts[j] = t->T_K;
      }
    }
    if (change < 0.01) break;
  }
  std::vector<double> X(N * m);
  for (std::size_t j = 0; j < N; ++j) {
    X[j * m] = Ts[j];
    auto K = k_values(c, Ts[j], col.P[j], x[j], notes);
    if (!K) return std::nullopt;
    std::vector<double> y(C);
    for (std::size_t i = 0; i < C; ++i) y[i] = (j == 0 && col.total) ? x[j][i] : (*K)[i] * x[j][i];
    normalize(y);
    for (std::size_t i = 0; i < C; ++i) {
      X[j * m + 1 + i] = L[j] * x[j][i];
      X[j * m + 1 + C + i] = V[j] * y[i];
    }
  }

  // ---- Newton (Naphtali–Sandholm) ----
  std::vector<double> R, dx;
  std::vector<StageProps> p;
  // Property warnings repeat for every stage and iteration; the Antoine-range warning names a temperature,
  // so those are reported once with the range of temperatures.
  auto finish_notes = [&]() {
    std::vector<std::string> seen;
    double tlo = 1e300, thi = 0.0;
    const std::string tag = "outside its Antoine range at ";
    for (const auto& it : notes.items()) {
      if (it.severity == Severity::Error) continue;
      const auto pos = it.message.find(tag);
      if (pos != std::string::npos) {
        const double t = std::atof(it.message.c_str() + pos + tag.size());
        tlo = std::min(tlo, t);
        thi = std::max(thi, t);
        continue;
      }
      if (std::find(seen.begin(), seen.end(), it.message) != seen.end()) continue;
      seen.push_back(it.message);
      diag.warn(it.code, it.message);
    }
    if (thi > 0.0)
      diag.warn(DiagCode::NoVaporPressureData, "vapour pressures were extrapolated outside their Antoine ranges between " +
                                                   num(std::round(tlo * 10) / 10) + " and " + num(std::round(thi * 10) / 10) + " K");
  };
  if (!col.residuals(X, R, p)) {
    finish_notes();
    diag.error(DiagCode::ColumnNotConverged, "the column equations could not be evaluated at the starting profile");
    return std::nullopt;
  }
  ColumnResult res;
  double norm = max_abs(R);
  bool converged = norm < 1e-9;
  int it = 0;
  for (; it < s.max_iterations && !converged; ++it) {
    if (!col.newton_step(X, R, p, dx)) break;
    // Backtracking line search on the largest residual; temperatures and flows are kept physical.
    double alpha = 1.0, best = 1e300;
    std::vector<double> Xbest, Rbest;
    std::vector<StageProps> pbest;
    for (int ls = 0; ls < 10; ++ls, alpha *= 0.5) {
      std::vector<double> Xt = X;
      for (std::size_t j = 0; j < N; ++j) {
        double& T = Xt[j * m];
        T = std::clamp(T + alpha * dx[j * m], T - 40.0, T + 40.0);
        T = std::max(T, 20.0);
        for (std::size_t k = 1; k < m; ++k) {
          double& v1 = Xt[j * m + k];
          const double nv = v1 + alpha * dx[j * m + k];
          v1 = nv >= 0.0 ? nv : 0.1 * v1;  // a flow that would turn negative drops to a tenth
        }
      }
      std::vector<double> Rt;
      std::vector<StageProps> pt;
      if (!col.residuals(Xt, Rt, pt)) continue;
      const double nt = max_abs(Rt);
      if (nt < best) best = nt, Xbest = Xt, Rbest = Rt, pbest = pt;
      if (nt < (1.0 - 1e-4 * alpha) * norm) break;
    }
    if (Xbest.empty()) break;
    X = std::move(Xbest);
    R = std::move(Rbest);
    p = std::move(pbest);
    norm = best;
    converged = norm < 1e-9;
  }
  finish_notes();
  res.iterations = it;
  res.residual = norm;
  if (!converged) {
    diag.error(DiagCode::ColumnNotConverged, "the column did not converge in " + std::to_string(it) +
                                                 " Newton iterations (largest residual " + num(norm) +
                                                 "); check that the specifications are feasible");
    return std::nullopt;
  }

  // ---- Results ----
  auto stage_vec = [&](std::size_t j, bool vap) {
    std::vector<double> n(C);
    for (std::size_t i = 0; i < C; ++i) n[i] = X[j * m + 1 + (vap ? C : 0) + i];
    return n;
  };
  auto frac = [](std::vector<double> n) {
    double t = 0.0;
    for (double v : n) t += v;
    for (double& v : n) v = t > 0.0 ? v / t : 0.0;
    return n;
  };
  auto H = [&](const std::vector<double>& n, const std::vector<double>& h) {
    double r = 0.0;
    for (std::size_t i = 0; i < C; ++i) r += n[i] * h[i];
    return r;
  };
  for (std::size_t j = 0; j < N; ++j) {
    ColumnStage st;
    st.T_K = X[j * m];
    st.P_Pa = col.P[j];
    const auto lj = stage_vec(j, false), vj = stage_vec(j, true);
    st.L = std::accumulate(lj.begin(), lj.end(), 0.0);
    st.x = frac(lj);
    if (j == 0 && col.total) {
      st.V = 0.0;
      st.y = st.x;
    } else {
      st.V = std::accumulate(vj.begin(), vj.end(), 0.0);
      st.y = frac(vj);
    }
    res.stages.push_back(std::move(st));
  }
  const auto l0 = stage_vec(0, false), v0 = stage_vec(0, true), v1 = stage_vec(1, true);
  const auto lN = stage_vec(N - 1, false), vN = stage_vec(N - 1, true), lN1 = stage_vec(N - 2, false);
  // Condenser: out − in.
  res.condenser_duty_W = H(l0, p[0].hL) + (col.total ? H(v0, p[0].hL) : H(v0, p[0].hV)) - H(v1, p[1].hV) - col.Hf[0];
  res.reboiler_duty_W = H(lN, p[N - 1].hL) + H(vN, p[N - 1].hV) - H(lN1, p[N - 2].hL) - col.Hf[N - 1];
  Stream& d = res.distillate;
  d.T_K = X[0];
  d.P_Pa = col.P[0];
  d.n = v0;
  d.vapor_fraction = col.total ? 0.0 : 1.0;
  d.x = d.y = frac(v0);
  d.H_W = H(v0, col.total ? p[0].hL : p[0].hV);
  Stream& b = res.bottoms;
  b.T_K = X[(N - 1) * m];
  b.P_Pa = col.P[N - 1];
  b.n = lN;
  b.vapor_fraction = 0.0;
  b.x = b.y = frac(lN);
  b.H_W = H(lN, p[N - 1].hL);
  return res;
}

}  // namespace chemlarp
