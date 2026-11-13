#include "chemlarp/md.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <random>

#include "chemlarp/conformer.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/json.hpp"
#include "chemlarp/system_state.hpp"
#include "chemlarp/uff.hpp"

namespace chemlarp {
namespace {

constexpr double kAcc = 4.184e-4;        // (kcal/mol/Å) / amu -> Å/fs²
constexpr double kBoltz = 0.0019872043;  // kcal/(mol K)
constexpr double kLightCmFs = 2.99792458e-5;  // speed of light, cm/fs

// Reproducible random numbers: mt19937_64 is fully specified; the deviates are computed here because
// std::normal_distribution and std::gamma_distribution differ between standard libraries.
class Rng {
 public:
  explicit Rng(std::uint64_t seed) : g_(seed) {}
  double uniform() { return (static_cast<double>(g_() >> 11) + 0.5) * 0x1.0p-53; }  // (0, 1)
  double normal() {
    if (have_) {
      have_ = false;
      return spare_;
    }
    const double u = uniform(), v = uniform();
    const double r = std::sqrt(-2.0 * std::log(u));
    spare_ = r * std::sin(2.0 * std::numbers::pi * v);
    have_ = true;
    return r * std::cos(2.0 * std::numbers::pi * v);
  }
  // Gamma(shape, 1): Marsaglia & Tsang, ACM TOMS 26 (2000) 363.
  double gamma(double shape) {
    if (shape < 1.0) return gamma(shape + 1.0) * std::pow(uniform(), 1.0 / shape);
    const double d = shape - 1.0 / 3.0, c = 1.0 / std::sqrt(9.0 * d);
    for (;;) {
      double x, v;
      do {
        x = normal();
        v = 1.0 + c * x;
      } while (v <= 0.0);
      v = v * v * v;
      const double u = uniform();
      if (u < 1.0 - 0.0331 * x * x * x * x) return d * v;
      if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return d * v;
    }
  }

 private:
  std::mt19937_64 g_;
  bool have_ = false;
  double spare_ = 0.0;
};

// Bussi stochastic velocity rescaling: a new kinetic energy for a canonical distribution with target mean
// Kbar over nf degrees of freedom, relaxing with c = exp(−dt/τ).
double bussi_kinetic(double K, double Kbar, double nf, double c, Rng& rng) {
  if (!(K > 0.0) || nf < 1.0) return K;
  const double r1 = rng.normal();
  const double sum = nf > 1.0 ? 2.0 * rng.gamma(0.5 * (nf - 1.0)) : 0.0;  // χ² with nf − 1 degrees
  return K + (1.0 - c) * (Kbar * (sum + r1 * r1) / nf - K) + 2.0 * r1 * std::sqrt(K * Kbar / nf * (1.0 - c) * c);
}

// Cyclic Jacobi eigen-decomposition of a symmetric n × n matrix (row-major); eigenvectors in columns of V.
void jacobi_eigen(std::vector<double> a, std::size_t n, std::vector<double>& w, std::vector<double>& V) {
  V.assign(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) V[i * n + i] = 1.0;
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0, diag = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      diag += a[i * n + i] * a[i * n + i];
      for (std::size_t j = i + 1; j < n; ++j) off += a[i * n + j] * a[i * n + j];
    }
    if (off <= 1e-30 * std::max(diag, 1e-300)) break;
    for (std::size_t p = 0; p < n; ++p)
      for (std::size_t q = p + 1; q < n; ++q) {
        const double apq = a[p * n + q];
        if (std::fabs(apq) < 1e-300) continue;
        const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (std::size_t k = 0; k < n; ++k) {
          const double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double vkp = V[k * n + p], vkq = V[k * n + q];
          V[k * n + p] = c * vkp - s * vkq;
          V[k * n + q] = s * vkp + c * vkq;
        }
      }
  }
  w.resize(n);
  for (std::size_t i = 0; i < n; ++i) w[i] = a[i * n + i];
}

// Orthonormal mass-weighted translation and rotation vectors (5 for a linear molecule, 6 otherwise).
std::vector<std::vector<double>> rigid_body_vectors(const std::vector<double>& m, const std::vector<double>& x) {
  const std::size_t N = m.size();
  double M = 0.0, c[3] = {0, 0, 0};
  for (std::size_t i = 0; i < N; ++i) {
    M += m[i];
    for (int a = 0; a < 3; ++a) c[a] += m[i] * x[3 * i + a];
  }
  for (double& v : c) v /= M;
  std::vector<std::vector<double>> raw;
  for (int a = 0; a < 3; ++a) {
    std::vector<double> t(3 * N, 0.0);
    for (std::size_t i = 0; i < N; ++i) t[3 * i + a] = std::sqrt(m[i]);
    raw.push_back(t);
  }
  // Rotations about the principal axes; a moment of inertia below 1e-3 of the largest (the axis of a linear
  // molecule, straight to within numerical precision) has no rotation.
  std::vector<double> I(9, 0.0);
  for (std::size_t i = 0; i < N; ++i) {
    const double p[3] = {x[3 * i] - c[0], x[3 * i + 1] - c[1], x[3 * i + 2] - c[2]};
    const double r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) I[a * 3 + b] += m[i] * ((a == b ? r2 : 0.0) - p[a] * p[b]);
  }
  std::vector<double> moments, axes;
  jacobi_eigen(I, 3, moments, axes);
  const double Imax = std::max({moments[0], moments[1], moments[2]});
  for (int ax = 0; ax < 3; ++ax) {
    if (!(moments[ax] > 1e-3 * Imax)) continue;
    const double e[3] = {axes[0 * 3 + ax], axes[1 * 3 + ax], axes[2 * 3 + ax]};
    std::vector<double> r(3 * N, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
      const double p[3] = {x[3 * i] - c[0], x[3 * i + 1] - c[1], x[3 * i + 2] - c[2]};
      const double cr[3] = {e[1] * p[2] - e[2] * p[1], e[2] * p[0] - e[0] * p[2], e[0] * p[1] - e[1] * p[0]};
      for (int k = 0; k < 3; ++k) r[3 * i + k] = std::sqrt(m[i]) * cr[k];
    }
    raw.push_back(r);
  }
  std::vector<std::vector<double>> out;
  for (auto& v : raw) {
    for (const auto& u : out) {
      double d = 0.0;
      for (std::size_t k = 0; k < v.size(); ++k) d += u[k] * v[k];
      for (std::size_t k = 0; k < v.size(); ++k) v[k] -= d * u[k];
    }
    double nrm = 0.0;
    for (double q : v) nrm += q * q;
    nrm = std::sqrt(nrm);
    if (nrm < 1e-6 * std::sqrt(M)) continue;  // a linear molecule has no rotation about its axis
    for (double& q : v) q /= nrm;
    out.push_back(v);
  }
  return out;
}

// Removes centre-of-mass velocity and overall angular momentum: projects the mass-weighted velocities
// off the rigid-body vectors.
void remove_rigid_motion(const std::vector<double>& m, const std::vector<double>& x, std::vector<double>& v) {
  const auto rb = rigid_body_vectors(m, x);
  std::vector<double> u(v.size());
  for (std::size_t i = 0; i < m.size(); ++i)
    for (int a = 0; a < 3; ++a) u[3 * i + a] = std::sqrt(m[i]) * v[3 * i + a];
  for (const auto& r : rb) {
    double d = 0.0;
    for (std::size_t k = 0; k < u.size(); ++k) d += r[k] * u[k];
    for (std::size_t k = 0; k < u.size(); ++k) u[k] -= d * r[k];
  }
  for (std::size_t i = 0; i < m.size(); ++i)
    for (int a = 0; a < 3; ++a) v[3 * i + a] = u[3 * i + a] / std::sqrt(m[i]);
}

double kinetic(const std::vector<double>& m, const std::vector<double>& v) {
  double K = 0.0;
  for (std::size_t i = 0; i < m.size(); ++i)
    K += m[i] * (v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] + v[3 * i + 2] * v[3 * i + 2]);
  return 0.5 * K / kAcc;
}

void linear_fit(const std::vector<double>& t, const std::vector<double>& y, double& slope, double& sd) {
  const double n = double(t.size());
  slope = 0.0;
  sd = 0.0;
  if (t.size() < 2) return;
  double mt = 0, my = 0;
  for (std::size_t i = 0; i < t.size(); ++i) mt += t[i], my += y[i];
  mt /= n;
  my /= n;
  double stt = 0, sty = 0, syy = 0;
  for (std::size_t i = 0; i < t.size(); ++i) {
    stt += (t[i] - mt) * (t[i] - mt);
    sty += (t[i] - mt) * (y[i] - my);
    syy += (y[i] - my) * (y[i] - my);
  }
  slope = stt > 0 ? sty / stt : 0.0;
  sd = std::sqrt(syy / n);
}

std::string jnum(double v) {
  char b[40];
  std::snprintf(b, sizeof b, "%.9g", v);
  return std::isfinite(v) ? b : "null";
}

}  // namespace

std::optional<MdResult> run_md(const std::vector<double>& masses, const std::vector<double>& x0,
                               const PotentialFunction& potential, const MdOptions& o, Diagnostics& diag,
                               const std::vector<double>* v0) {
  const std::size_t N = masses.size();
  auto bad = [&](const std::string& m) {
    diag.error(DiagCode::InvalidInput, m);
    return std::nullopt;
  };
  if (N == 0 || x0.size() != 3 * N) return bad("molecular dynamics needs one position (x, y, z) per mass");
  for (double m : masses)
    if (!(m > 0.0)) return bad("masses must be positive");
  if (!(o.dt_fs > 0.0) || o.dt_fs > 10.0) return bad("the time step must be between 0 and 10 fs");
  if (!(o.temperature_K >= 0.0)) return bad("the temperature must not be negative");
  if (o.thermostat != Thermostat::None && !(o.tau_fs > 0.0)) return bad("the thermostat time constant must be positive");
  if (o.sample_every == 0) return bad("sample_every must be at least 1");

  Rng rng(o.seed);
  MdResult res;
  std::vector<double> x = x0, v(3 * N, 0.0), g(3 * N, 0.0);
  const bool langevin = o.thermostat == Thermostat::Langevin;
  const std::size_t rigid = (N > 1 && !langevin) ? rigid_body_vectors(masses, x).size() : (langevin ? 0 : 3);
  res.degrees_of_freedom = 3 * N > rigid ? 3 * N - rigid : 0;
  const double nf = double(res.degrees_of_freedom);
  const double kT = kBoltz * o.temperature_K;
  if (v0) {
    if (v0->size() != 3 * N) return bad("initial velocities need one entry per coordinate");
    v = *v0;
  } else {
    for (std::size_t i = 0; i < N; ++i) {
      const double s = std::sqrt(kT * kAcc / masses[i]);
      for (int a = 0; a < 3; ++a) v[3 * i + a] = s * rng.normal();
    }
  }
  if (!langevin && N > 1) remove_rigid_motion(masses, x, v);
  if (!v0 && nf > 0) {
    const double K = kinetic(masses, v);
    if (K > 0) {
      const double f = std::sqrt(0.5 * nf * kT / K);
      for (double& q : v) q *= f;
    }
  }
  double U = potential(x, g);
  if (!std::isfinite(U)) return bad("the potential energy is not finite at the starting coordinates");
  const double dt = o.dt_fs;
  auto accel = [&](std::vector<double>& a) {
    a.resize(3 * N);
    for (std::size_t i = 0; i < N; ++i)
      for (int k = 0; k < 3; ++k) a[3 * i + k] = -g[3 * i + k] * kAcc / masses[i];
  };
  std::vector<double> acc;
  accel(acc);
  double heat = 0.0;  // energy added by the thermostat
  const double c_bussi = std::exp(-dt / o.tau_fs);
  const double c1 = langevin ? std::exp(-dt / o.tau_fs) : 1.0, c2 = std::sqrt(1.0 - c1 * c1);
  auto record = [&](std::size_t step) {
    MdFrame f;
    f.time_fs = double(step) * dt;
    f.kinetic = kinetic(masses, v);
    f.potential = U;
    f.total = f.kinetic + U;
    f.conserved = f.total - heat;
    f.temperature_K = nf > 0 ? 2.0 * f.kinetic / (nf * kBoltz) : 0.0;
    if (o.keep_positions) f.x = x;
    if (o.keep_velocities) f.v = v;
    res.frames.push_back(std::move(f));
  };
  record(0);
  for (std::size_t step = 1; step <= o.steps; ++step) {
    if (langevin) {
      // BAOAB: half kick, half drift, friction and noise, half drift, force, half kick.
      for (std::size_t k = 0; k < 3 * N; ++k) v[k] += 0.5 * dt * acc[k];
      for (std::size_t k = 0; k < 3 * N; ++k) x[k] += 0.5 * dt * v[k];
      const double K0 = kinetic(masses, v);
      for (std::size_t i = 0; i < N; ++i) {
        const double s = std::sqrt(kT * kAcc / masses[i]);
        for (int a = 0; a < 3; ++a) v[3 * i + a] = c1 * v[3 * i + a] + c2 * s * rng.normal();
      }
      heat += kinetic(masses, v) - K0;
      for (std::size_t k = 0; k < 3 * N; ++k) x[k] += 0.5 * dt * v[k];
      U = potential(x, g);
      accel(acc);
      for (std::size_t k = 0; k < 3 * N; ++k) v[k] += 0.5 * dt * acc[k];
    } else {
      for (std::size_t k = 0; k < 3 * N; ++k) {
        v[k] += 0.5 * dt * acc[k];
        x[k] += dt * v[k];
      }
      U = potential(x, g);
      accel(acc);
      for (std::size_t k = 0; k < 3 * N; ++k) v[k] += 0.5 * dt * acc[k];
      if (o.thermostat == Thermostat::Bussi && nf > 0) {
        const double K = kinetic(masses, v);
        const double Kn = bussi_kinetic(K, 0.5 * nf * kT, nf, c_bussi, rng);
        if (K > 0 && Kn > 0) {
          const double f = std::sqrt(Kn / K);
          for (double& q : v) q *= f;
          heat += Kn - K;
        }
      }
    }
    if (!std::isfinite(U) || !(std::fabs(U) < 1e12)) {
      diag.error(DiagCode::NotConverged, "the dynamics became unstable at step " + std::to_string(step) +
                                             "; use a shorter time step");
      return std::nullopt;
    }
    if (step % o.sample_every == 0) record(step);
  }
  res.x = x;
  res.v = v;
  std::vector<double> t, e;
  double sT = 0, sU = 0;
  for (const auto& f : res.frames) {
    t.push_back(f.time_fs * 1e-3);  // ps
    e.push_back(f.conserved);
    sT += f.temperature_K;
    sU += f.potential;
  }
  res.mean_temperature_K = sT / double(res.frames.size());
  res.mean_potential = sU / double(res.frames.size());
  linear_fit(t, e, res.drift_kcal_mol_ps, res.fluctuation);
  return res;
}

std::optional<NormalModes> normal_modes(const std::vector<double>& m, const std::vector<double>& x,
                                        const PotentialFunction& potential, Diagnostics& diag) {
  const std::size_t N = m.size(), n = 3 * N;
  if (N < 2 || x.size() != n) {
    diag.error(DiagCode::InvalidInput, "normal modes need at least two atoms");
    return std::nullopt;
  }
  if (N > 150) {
    diag.error(DiagCode::InvalidInput, "normal modes are limited to 150 atoms");
    return std::nullopt;
  }
  const double h = 1e-4;
  std::vector<double> H(n * n), gp(n), gm(n), xx = x;
  for (std::size_t k = 0; k < n; ++k) {
    xx[k] = x[k] + h;
    potential(xx, gp);
    xx[k] = x[k] - h;
    potential(xx, gm);
    xx[k] = x[k];
    for (std::size_t j = 0; j < n; ++j) H[j * n + k] = (gp[j] - gm[j]) / (2 * h);
  }
  std::vector<double> sm(n);
  for (std::size_t i = 0; i < N; ++i)
    for (int a = 0; a < 3; ++a) sm[3 * i + a] = std::sqrt(m[i]);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i; j < n; ++j) {
      const double s = 0.5 * (H[i * n + j] + H[j * n + i]) / (sm[i] * sm[j]);
      H[i * n + j] = H[j * n + i] = s;
    }
  // Project out translation and rotation: H ← P H P with P = 1 − Σ u uᵀ.
  const auto rb = rigid_body_vectors(m, x);
  std::vector<double> P(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) P[i * n + i] = 1.0;
  for (const auto& u : rb)
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j) P[i * n + j] -= u[i] * u[j];
  std::vector<double> T(n * n, 0.0), Hp(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t k = 0; k < n; ++k) {
      const double pik = P[i * n + k];
      if (pik == 0.0) continue;
      for (std::size_t j = 0; j < n; ++j) T[i * n + j] += pik * H[k * n + j];
    }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t k = 0; k < n; ++k) {
      const double tik = T[i * n + k];
      if (tik == 0.0) continue;
      for (std::size_t j = 0; j < n; ++j) Hp[i * n + j] += tik * P[k * n + j];
    }
  std::vector<double> w, V;
  jacobi_eigen(Hp, n, w, V);
  // Drop the eigenvectors that lie in the rigid-body space.
  struct Mode {
    double lambda;
    std::size_t col;
  };
  std::vector<Mode> modes;
  for (std::size_t c = 0; c < n; ++c) {
    double in = 0.0;
    for (const auto& u : rb) {
      double d = 0.0;
      for (std::size_t k = 0; k < n; ++k) d += u[k] * V[k * n + c];
      in += d * d;
    }
    if (in < 0.5) modes.push_back({w[c], c});
  }
  std::sort(modes.begin(), modes.end(), [](const Mode& a, const Mode& b) { return a.lambda < b.lambda; });
  NormalModes out;
  for (const auto& md : modes) {
    const double nu = std::sqrt(std::fabs(md.lambda) * kAcc) / (2.0 * std::numbers::pi * kLightCmFs);
    out.frequencies_cm1.push_back(md.lambda < 0 ? -nu : nu);
    std::vector<double> d(n);
    double nrm = 0.0, inv_mu = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
      d[k] = V[k * n + md.col] / sm[k];
      nrm += d[k] * d[k];
      inv_mu += V[k * n + md.col] * V[k * n + md.col] / (sm[k] * sm[k]);
    }
    for (double& q : d) q /= std::sqrt(nrm);
    out.displacements.push_back(d);
    out.reduced_masses.push_back(inv_mu > 0 ? 1.0 / inv_mu : 0.0);
  }
  return out;
}

std::vector<double> velocity_spectrum(const MdResult& r, const std::vector<double>& m, double sample_dt_fs,
                                      double max_cm1, std::size_t points) {
  std::vector<double> out(points, 0.0);
  const std::size_t T = r.frames.size();
  if (T < 4 || points < 2) return out;
  const std::size_t n = r.frames[0].v.size();
  if (n == 0 || n != 3 * m.size()) return out;
  std::vector<double> win(T);
  for (std::size_t t = 0; t < T; ++t) win[t] = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * double(t) / double(T - 1));
  double peak = 0.0;
  for (std::size_t p = 0; p < points; ++p) {
    const double nu = max_cm1 * double(p) / double(points - 1);
    const double om = 2.0 * std::numbers::pi * kLightCmFs * nu * sample_dt_fs;  // rad per sample
    double s = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
      const double sq = std::sqrt(m[k / 3]);
      double re = 0.0, im = 0.0;
      for (std::size_t t = 0; t < T; ++t) {
        const double val = win[t] * sq * r.frames[t].v[k];
        re += val * std::cos(om * double(t));
        im -= val * std::sin(om * double(t));
      }
      s += re * re + im * im;
    }
    out[p] = s;
    peak = std::max(peak, s);
  }
  if (peak > 0)
    for (double& v : out) v /= peak;
  return out;
}

std::optional<MoleculeMdResult> molecule_md(std::string_view smiles, const MdOptions& options, bool compute_modes,
                                            Diagnostics& diag) {
  ChemicalSystemState st;
  EmbedResult er;
  const std::int32_t sp = embed_smiles(st, smiles, diag, &er);
  if (sp < 0) return std::nullopt;
  const auto s = static_cast<std::uint32_t>(sp);
  auto ff = UffForceField::build(st, s, diag);
  if (!ff) return std::nullopt;
  MoleculeMdResult res;
  const auto& spc = st.species[s];
  const std::uint32_t b0 = spc.atom_begin;
  std::vector<double> x;
  for (std::uint32_t a = spc.atom_begin; a < spc.atom_begin + spc.atom_count; ++a) {
    const ElementInfo* e = find_element_by_z(st.atoms.element[a]);
    res.elements.push_back(e ? std::string(e->symbol) : "X");
    res.masses.push_back(e ? e->mass : 1.0);
    x.push_back(st.atoms.x[a]);
    x.push_back(st.atoms.y[a]);
    x.push_back(st.atoms.z[a]);
  }
  for (std::size_t b = 0; b < st.bonds.size(); ++b)
    if (st.bonds.species[b] == s) res.bonds.push_back({st.bonds.begin[b] - b0, st.bonds.end[b] - b0, int(st.bonds.order[b])});
  if (res.masses.size() > 400) {
    diag.error(DiagCode::InvalidInput, "molecular dynamics of a single molecule is limited to 400 atoms");
    return std::nullopt;
  }
  const UffForceField& f = *ff;
  PotentialFunction pot = [&f](const std::vector<double>& xx, std::vector<double>& g) { return f.energy_gradient(xx, g); };
  res.minimized_energy = f.energy(x);
  if (compute_modes && res.masses.size() <= 150) res.modes = normal_modes(res.masses, x, pot, diag);
  auto md = run_md(res.masses, x, pot, options, diag);
  if (!md) return std::nullopt;
  res.md = std::move(*md);
  return res;
}

std::string molecule_md_json(std::string_view smiles, const MdOptions& options, bool compute_modes) {
  Diagnostics d;
  MdOptions o = options;
  o.keep_positions = true;
  o.keep_velocities = true;
  auto r = molecule_md(smiles, o, compute_modes, d);
  std::string out = std::string("{\"ok\":") + (r && !d.has_errors() ? "true" : "false");
  if (r) {
    out += ",\"elements\":[";
    for (std::size_t i = 0; i < r->elements.size(); ++i) out += (i ? "," : "") + json_escape(r->elements[i]);
    out += "],\"bonds\":[";
    for (std::size_t i = 0; i < r->bonds.size(); ++i)
      out += (i ? ",[" : "[") + std::to_string(r->bonds[i].a) + "," + std::to_string(r->bonds[i].b) + "," +
             std::to_string(r->bonds[i].order) + "]";
    const auto& m = r->md;
    out += "],\"minimized_energy\":" + jnum(r->minimized_energy) + ",\"dof\":" + std::to_string(m.degrees_of_freedom) +
           ",\"mean_temperature_K\":" + jnum(m.mean_temperature_K) + ",\"mean_potential\":" + jnum(m.mean_potential) +
           ",\"drift_kcal_mol_ps\":" + jnum(m.drift_kcal_mol_ps) + ",\"fluctuation\":" + jnum(m.fluctuation) + ",\"frames\":[";
    for (std::size_t i = 0; i < m.frames.size(); ++i) {
      const auto& f = m.frames[i];
      out += (i ? ",{" : "{") + std::string("\"t\":") + jnum(f.time_fs) + ",\"K\":" + jnum(f.kinetic) + ",\"U\":" +
             jnum(f.potential) + ",\"E\":" + jnum(f.total) + ",\"conserved\":" + jnum(f.conserved) + ",\"T\":" +
             jnum(f.temperature_K) + "}";
    }
    // Positions, at most 400 frames (evenly thinned), to 4 decimals.
    out += "],\"positions\":[";
    const std::size_t nf = m.frames.size(), stride = std::max<std::size_t>(1, (nf + 399) / 400);
    bool first = true;
    for (std::size_t i = 0; i < nf; i += stride) {
      out += first ? "[" : ",[";
      first = false;
      const auto& x = m.frames[i].x;
      for (std::size_t k = 0; k < x.size(); ++k) {
        char b[24];
        std::snprintf(b, sizeof b, "%s%.4f", k ? "," : "", x[k]);
        out += b;
      }
      out += "]";
    }
    out += "],\"position_stride\":" + std::to_string(stride) + ",\"modes\":[";
    if (r->modes)
      for (std::size_t i = 0; i < r->modes->frequencies_cm1.size(); ++i) out += (i ? "," : "") + jnum(r->modes->frequencies_cm1[i]);
    const double maxcm = 4000.0;
    const auto spec = velocity_spectrum(m, r->masses, o.dt_fs * double(o.sample_every), maxcm, 801);
    out += "],\"spectrum\":{\"max_cm1\":" + jnum(maxcm) + ",\"intensity\":[";
    for (std::size_t i = 0; i < spec.size(); ++i) {
      char b[24];
      std::snprintf(b, sizeof b, "%s%.5g", i ? "," : "", spec[i]);
      out += b;
    }
    // Multi-frame XYZ of the stored positions.
    std::string xyz;
    for (std::size_t i = 0; i < nf; i += stride) {
      xyz += std::to_string(r->elements.size()) + "\nt = " + jnum(m.frames[i].time_fs) + " fs\n";
      const auto& x = m.frames[i].x;
      for (std::size_t a = 0; a < r->elements.size(); ++a) {
        char b[96];
        std::snprintf(b, sizeof b, "%-2s %12.6f %12.6f %12.6f\n", r->elements[a].c_str(), x[3 * a], x[3 * a + 1], x[3 * a + 2]);
        xyz += b;
      }
    }
    out += "]},\"xyz\":" + json_escape(xyz);
  }
  return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

std::optional<LjResult> lennard_jones_fluid(const LjOptions& o, Diagnostics& diag) {
  if (o.cells < 1 || o.cells > 8 || !(o.density > 0 && o.density < 1.5) || !(o.temperature > 0) || !(o.cutoff > 1.0) ||
      !(o.dt > 0 && o.dt < 0.02)) {
    diag.error(DiagCode::InvalidInput, "Lennard-Jones fluid: 1-8 cells, 0 < density < 1.5, T > 0, cutoff > 1, 0 < dt < 0.02");
    return std::nullopt;
  }
  const std::size_t N = std::size_t(4 * o.cells * o.cells * o.cells);
  const double L = std::cbrt(double(N) / o.density);
  if (o.cutoff > 0.5 * L) {
    diag.error(DiagCode::InvalidInput, "the cutoff must not exceed half the box length (use more cells)");
    return std::nullopt;
  }
  LjResult res;
  res.particles = N;
  res.box = L;
  std::vector<double> x(3 * N), v(3 * N), f(3 * N);
  const double a = L / o.cells;
  const double basis[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  std::size_t p = 0;
  for (int i = 0; i < o.cells; ++i)
    for (int j = 0; j < o.cells; ++j)
      for (int k = 0; k < o.cells; ++k)
        for (const auto& b : basis) {
          x[3 * p] = (i + b[0] + 0.25) * a;
          x[3 * p + 1] = (j + b[1] + 0.25) * a;
          x[3 * p + 2] = (k + b[2] + 0.25) * a;
          ++p;
        }
  Rng rng(o.seed);
  for (double& q : v) q = std::sqrt(o.temperature) * rng.normal();
  double vc[3] = {0, 0, 0};
  for (std::size_t i = 0; i < N; ++i)
    for (int d = 0; d < 3; ++d) vc[d] += v[3 * i + d] / double(N);
  for (std::size_t i = 0; i < N; ++i)
    for (int d = 0; d < 3; ++d) v[3 * i + d] -= vc[d];
  const double rc2 = o.cutoff * o.cutoff;
  const double ir6c = 1.0 / (rc2 * rc2 * rc2), ushift = 4.0 * ir6c * (ir6c - 1.0);
  double U = 0.0, W = 0.0;
  auto forces = [&]() {
    std::fill(f.begin(), f.end(), 0.0);
    U = 0.0;
    W = 0.0;
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = i + 1; j < N; ++j) {
        double d[3], r2 = 0.0;
        for (int k = 0; k < 3; ++k) {
          d[k] = x[3 * i + k] - x[3 * j + k];
          d[k] -= L * std::round(d[k] / L);
          r2 += d[k] * d[k];
        }
        if (r2 >= rc2) continue;
        const double ir2 = 1.0 / r2, ir6 = ir2 * ir2 * ir2;
        U += 4.0 * ir6 * (ir6 - 1.0) - ushift;
        const double ff = 48.0 * ir6 * (ir6 - 0.5) * ir2;  // F·r / r²
        W += ff * r2;
        for (int k = 0; k < 3; ++k) {
          f[3 * i + k] += ff * d[k];
          f[3 * j + k] -= ff * d[k];
        }
      }
  };
  auto kin = [&]() {
    double K = 0.0;
    for (double q : v) K += q * q;
    return 0.5 * K;
  };
  forces();
  const double nf = 3.0 * double(N) - 3.0, dt = o.dt;
  const bool langevin = o.thermostat == Thermostat::Langevin;
  const double cb = std::exp(-dt / o.tau), c1 = std::exp(-dt / o.tau), c2 = std::sqrt(1 - c1 * c1);
  double heat = 0.0;
  std::vector<double> us, ps, ts, tt, ec;
  const double V = L * L * L;
  for (std::size_t step = 1; step <= o.equilibration + o.production; ++step) {
    if (langevin) {
      for (std::size_t k = 0; k < 3 * N; ++k) v[k] += 0.5 * dt * f[k], x[k] += 0.5 * dt * v[k];
      const double K0 = kin();
      for (double& q : v) q = c1 * q + c2 * std::sqrt(o.temperature) * rng.normal();
      heat += kin() - K0;
      for (std::size_t k = 0; k < 3 * N; ++k) x[k] += 0.5 * dt * v[k];
      forces();
      for (std::size_t k = 0; k < 3 * N; ++k) v[k] += 0.5 * dt * f[k];
    } else {
      for (std::size_t k = 0; k < 3 * N; ++k) v[k] += 0.5 * dt * f[k], x[k] += dt * v[k];
      forces();
      for (std::size_t k = 0; k < 3 * N; ++k) v[k] += 0.5 * dt * f[k];
      if (o.thermostat == Thermostat::Bussi) {
        const double K = kin(), Kn = bussi_kinetic(K, 0.5 * nf * o.temperature, nf, cb, rng);
        if (K > 0 && Kn > 0) {
          const double s = std::sqrt(Kn / K);
          for (double& q : v) q *= s;
          heat += Kn - K;
        }
      }
    }
    if (!std::isfinite(U)) {
      diag.error(DiagCode::NotConverged, "the Lennard-Jones dynamics became unstable; use a shorter time step");
      return std::nullopt;
    }
    if (step > o.equilibration) {
      const double K = kin(), Tk = 2.0 * K / nf;
      us.push_back(U / double(N));
      // Canonical ideal-gas term from the set temperature (MD conserves momentum: 3N − 3 kinetic degrees).
      ps.push_back((double(N) * (o.thermostat == Thermostat::None ? Tk : o.temperature) + W / 3.0) / V);
      ts.push_back(Tk);
      tt.push_back(double(step) * dt);
      ec.push_back((K + U - heat) / double(N));
    }
  }
  auto block = [](const std::vector<double>& y, double& mean, double& err) {
    const std::size_t nb = 20, bs = y.size() / nb;
    mean = 0.0;
    for (double q : y) mean += q;
    mean /= double(y.size());
    err = 0.0;
    if (bs == 0) return;
    std::vector<double> bm(nb, 0.0);
    for (std::size_t b = 0; b < nb; ++b) {
      for (std::size_t i = 0; i < bs; ++i) bm[b] += y[b * bs + i];
      bm[b] /= double(bs);
    }
    double bmean = 0.0;
    for (double q : bm) bmean += q / double(nb);
    double var = 0.0;
    for (double q : bm) var += (q - bmean) * (q - bmean);
    err = std::sqrt(var / double(nb - 1) / double(nb));
  };
  double terr = 0.0;
  block(us, res.energy_per_particle, res.energy_error);
  block(ps, res.pressure, res.pressure_error);
  block(ts, res.temperature, terr);
  double sd = 0.0;
  linear_fit(tt, ec, res.drift, sd);
  return res;
}

}  // namespace chemlarp
