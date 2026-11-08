#include "omnichem/gibbs.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "omnichem/elements.hpp"
#include "omnichem/thermo.hpp"

namespace omnichem {

bool solve_dense(std::vector<double>& m, std::vector<double>& rhs, std::size_t n) {
  double scale = 0.0;
  for (double v : m) scale = std::max(scale, std::fabs(v));
  if (scale == 0.0) return false;
  for (std::size_t col = 0; col < n; ++col) {
    std::size_t piv = col;
    for (std::size_t r = col + 1; r < n; ++r)
      if (std::fabs(m[r * n + col]) > std::fabs(m[piv * n + col])) piv = r;
    if (std::fabs(m[piv * n + col]) < 1e-14 * scale) return false;
    if (piv != col) {
      for (std::size_t c = 0; c < n; ++c) std::swap(m[col * n + c], m[piv * n + c]);
      std::swap(rhs[col], rhs[piv]);
    }
    const double inv = 1.0 / m[col * n + col];
    for (std::size_t r = col + 1; r < n; ++r) {
      const double f = m[r * n + col] * inv;
      if (f == 0.0) continue;
      for (std::size_t c = col; c < n; ++c) m[r * n + c] -= f * m[col * n + c];
      rhs[r] -= f * rhs[col];
    }
  }
  for (std::size_t i = n; i-- > 0;) {
    double s = rhs[i];
    for (std::size_t c = i + 1; c < n; ++c) s -= m[i * n + c] * rhs[c];
    rhs[i] = s / m[i * n + i];
  }
  return true;
}

EquilibriumResult minimize_gibbs_ideal_gas(const EquilibriumInput& in, Diagnostics& diag,
                                           const EquilibriumOptions& opt) {
  EquilibriumResult r;
  const std::size_t S = in.g_RT.size();
  if (in.composition.size() != S || in.n0.size() != S) {
    diag.error(DiagCode::InvalidInput, "composition, g_RT and n0 must have one entry per species");
    return r;
  }
  if (!(in.pressure_Pa > 0.0) || !std::isfinite(in.pressure_Pa)) {
    diag.error(DiagCode::InvalidInput, "pressure must be positive and finite");
    return r;
  }
  bool bad = false;
  for (std::size_t i = 0; i < S; ++i) {
    if (!(in.n0[i] >= 0.0) || !std::isfinite(in.n0[i])) {
      diag.error(DiagCode::NegativeAmount, "initial amount of species " + std::to_string(i) + " must be finite and >= 0")
          .on_species(static_cast<std::int64_t>(i));
      bad = true;
    }
    if (!std::isfinite(in.g_RT[i])) {
      diag.error(DiagCode::NoThermoData, "species " + std::to_string(i) + " has no standard chemical potential")
          .on_species(static_cast<std::int64_t>(i));
      bad = true;
    }
  }
  if (bad) return r;
  r.n.assign(S, 0.0);
  r.mole_fraction.assign(S, 0.0);

  // Element totals b = A n0.
  std::array<double, kElementSlots> btot{};
  for (std::size_t i = 0; i < S; ++i)
    for (std::size_t z = 1; z < kElementSlots; ++z) btot[z] += in.composition[i][z] * in.n0[i];

  // A species can only exist if every element it contains is present.
  std::vector<std::size_t> active;
  for (std::size_t i = 0; i < S; ++i) {
    bool ok = false;
    for (std::size_t z = 1; z < kElementSlots; ++z) {
      if (in.composition[i][z] == 0) continue;
      ok = true;
      if (btot[z] <= 0.0) { ok = false; break; }
    }
    if (ok) active.push_back(i);
  }
  std::vector<std::uint8_t> elems;
  for (std::size_t z = 1; z < kElementSlots; ++z)
    if (btot[z] > 0.0) elems.push_back(static_cast<std::uint8_t>(z));
  if (active.empty() || elems.empty()) {
    diag.error(DiagCode::EmptySystem, "no species with a positive amount of every constituent element");
    return r;
  }

  // Keep a linearly independent subset of element rows (modified Gram-Schmidt).
  // Example: a system containing only H2O has H = 2*O; one row must be dropped.
  const std::size_t Sa = active.size();
  std::vector<std::vector<double>> basis;
  std::vector<std::uint8_t> rows;
  for (std::uint8_t z : elems) {
    std::vector<double> v(Sa);
    for (std::size_t k = 0; k < Sa; ++k) v[k] = in.composition[active[k]][z];
    double norm0 = 0.0;
    for (double x : v) norm0 += x * x;
    for (const auto& q : basis) {
      double dot = 0.0;
      for (std::size_t k = 0; k < Sa; ++k) dot += v[k] * q[k];
      for (std::size_t k = 0; k < Sa; ++k) v[k] -= dot * q[k];
    }
    double norm = 0.0;
    for (double x : v) norm += x * x;
    if (norm > 1e-20 * std::max(norm0, 1.0)) {
      norm = std::sqrt(norm);
      for (double& x : v) x /= norm;
      basis.push_back(std::move(v));
      rows.push_back(z);
    }
  }
  const std::size_t E = rows.size();
  std::vector<double> A(E * Sa);  // row-major E x Sa
  std::vector<double> b(E);
  for (std::size_t j = 0; j < E; ++j) {
    b[j] = btot[rows[j]];
    for (std::size_t k = 0; k < Sa; ++k) A[j * Sa + k] = in.composition[active[k]][rows[j]];
  }
  std::vector<double> c(Sa);
  const double lnP = std::log(in.pressure_Pa / kStandardPressure);
  for (std::size_t k = 0; k < Sa; ++k) c[k] = in.g_RT[active[k]] + lnP;

  // Unknowns: pi (E element potentials) and nu = ln N.
  //   x_k = exp(a_k . pi - c_k),  n_k = N x_k
  //   F_j = (sum_k a_jk n_k - b_j) / b_j     (j < E)
  //   F_E = sum_k x_k - 1
  const std::size_t U = E + 1;
  std::vector<double> u(U, 0.0);
  {
    // Initial guess: least-squares fit of ln x_k to the starting composition.
    double N0 = 0.0;
    for (std::size_t k = 0; k < Sa; ++k) N0 += in.n0[active[k]];
    std::vector<double> M(E * E, 0.0), rhs(E, 0.0);
    for (std::size_t k = 0; k < Sa; ++k) {
      const double xg = std::max(in.n0[active[k]] / N0, 1e-6);
      const double target = c[k] + std::log(xg);
      for (std::size_t j = 0; j < E; ++j) {
        rhs[j] += A[j * Sa + k] * target;
        for (std::size_t l = 0; l < E; ++l) M[j * E + l] += A[j * Sa + k] * A[l * Sa + k];
      }
    }
    if (!solve_dense(M, rhs, E)) {
      diag.error(DiagCode::SingularSystem, "element matrix is singular");
      return r;
    }
    for (std::size_t j = 0; j < E; ++j) u[j] = rhs[j];
    u[E] = std::log(N0);
  }

  std::vector<double> x(Sa), F(U), J(U * U);
  auto evaluate = [&](const std::vector<double>& uu, bool jac) {
    const double N = std::exp(uu[E]);
    for (std::size_t k = 0; k < Sa; ++k) {
      double e = -c[k];
      for (std::size_t j = 0; j < E; ++j) e += A[j * Sa + k] * uu[j];
      x[k] = std::exp(std::min(e, 700.0));
    }
    double sx = 0.0;
    for (double v : x) sx += v;
    for (std::size_t j = 0; j < E; ++j) {
      double s = 0.0;
      for (std::size_t k = 0; k < Sa; ++k) s += A[j * Sa + k] * x[k];
      F[j] = (N * s - b[j]) / b[j];
    }
    F[E] = sx - 1.0;
    if (jac) {
      std::fill(J.begin(), J.end(), 0.0);
      for (std::size_t j = 0; j < E; ++j) {
        for (std::size_t l = 0; l < E; ++l) {
          double s = 0.0;
          for (std::size_t k = 0; k < Sa; ++k) s += A[j * Sa + k] * A[l * Sa + k] * x[k];
          J[j * U + l] = N * s / b[j];
        }
        double s = 0.0;
        for (std::size_t k = 0; k < Sa; ++k) s += A[j * Sa + k] * x[k];
        J[j * U + E] = N * s / b[j];
        J[E * U + j] = s;
      }
      J[E * U + E] = 0.0;
    }
    double norm2 = 0.0;
    for (double f : F) norm2 += f * f;
    return norm2;
  };

  double norm2 = evaluate(u, true);
  for (r.iterations = 0; r.iterations < opt.max_iterations; ++r.iterations) {
    double inf = 0.0;
    for (double f : F) inf = std::max(inf, std::fabs(f));
    if (inf < opt.tolerance) { r.converged = true; break; }

    std::vector<double> M = J, step(U);
    for (std::size_t i = 0; i < U; ++i) step[i] = -F[i];
    if (!solve_dense(M, step, U)) {
      diag.error(DiagCode::SingularSystem, "singular Newton matrix in Gibbs minimisation");
      return r;
    }
    double biggest = 0.0;
    for (double s : step) biggest = std::max(biggest, std::fabs(s));
    if (biggest > opt.max_step)
      for (double& s : step) s *= opt.max_step / biggest;

    // Backtracking line search on the residual norm.
    double lambda = 1.0;
    std::vector<double> trial(U);
    double trial_norm2 = norm2;
    for (int bt = 0; bt < 40; ++bt) {
      for (std::size_t i = 0; i < U; ++i) trial[i] = u[i] + lambda * step[i];
      trial_norm2 = evaluate(trial, false);
      if (std::isfinite(trial_norm2) && trial_norm2 <= (1.0 - 1e-4 * lambda) * norm2) break;
      lambda *= 0.5;
    }
    u = trial;
    norm2 = evaluate(u, true);
  }

  const double N = std::exp(u[E]);
  r.total_moles = N;
  r.elements = rows;
  r.element_potentials.assign(u.begin(), u.begin() + static_cast<std::ptrdiff_t>(E));
  for (std::size_t k = 0; k < Sa; ++k) {
    r.n[active[k]] = N * x[k];
    r.mole_fraction[active[k]] = x[k];
  }
  r.G_RT = 0.0;
  for (std::size_t k = 0; k < Sa; ++k)
    if (x[k] > 0.0) r.G_RT += N * x[k] * (c[k] + std::log(x[k]));
  for (std::size_t z = 1; z < kElementSlots; ++z) {
    if (btot[z] <= 0.0) continue;
    double s = 0.0;
    for (std::size_t i = 0; i < S; ++i) s += in.composition[i][z] * r.n[i];
    r.max_relative_element_error = std::max(r.max_relative_element_error, std::fabs(s - btot[z]) / btot[z]);
  }
  if (!r.converged) {
    diag.error(DiagCode::NotConverged, "Gibbs minimisation did not converge in " + std::to_string(opt.max_iterations) +
                                           " iterations (residual " + std::to_string(std::sqrt(norm2)) + ")");
  }
  return r;
}

EquilibriumResult equilibrate(ChemicalSystemState& state, Diagnostics& diag, const EquilibriumOptions& options) {
  bool ok = true;
  for (std::uint32_t i = 0; i < state.species.size(); ++i) {
    if (state.species[i].net_charge != 0) {
      diag.error(DiagCode::ChargedSpeciesUnsupported,
                 "species " + state.species[i].formula + " is charged; ionic equilibria need an electroneutrality "
                 "constraint and activity model (not implemented)")
          .on_species(i);
      ok = false;
    }
  }
  if (!ok) return {};
  if (!bind_thermo(state, diag)) return {};
  EquilibriumInput in;
  in.pressure_Pa = state.conditions.pressure_Pa;
  for (std::uint32_t i = 0; i < state.species.size(); ++i) {
    in.composition.push_back(state.species[i].composition);
    in.g_RT.push_back(state.g_RT[i]);
    in.n0.push_back(state.moles[i]);
  }
  EquilibriumResult r = minimize_gibbs_ideal_gas(in, diag, options);
  if (r.converged) {
    for (std::size_t i = 0; i < r.n.size(); ++i) state.moles[i] = r.n[i];
    ++state.revision;
  }
  return r;
}

}  // namespace omnichem
