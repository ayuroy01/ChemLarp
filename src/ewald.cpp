// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#include "chemlarp/ewald.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "chemlarp/parallel.hpp"

namespace chemlarp {

namespace {

constexpr double kPi = std::numbers::pi;

// Cardinal B-spline weights of order n at fractional offset w in [0, 1): data[j] applies to grid point
// floor(u) − n + 1 + j; ddata[j] = d data[j] / dw (Essmann et al. 1995, appendix).
void bspline(double w, int n, double* data, double* ddata) {
  data[n - 1] = 0.0;
  data[1] = w;
  data[0] = 1.0 - w;
  for (int j = 3; j < n; ++j) {
    const double div = 1.0 / (j - 1);
    data[j - 1] = div * w * data[j - 2];
    for (int k = 1; k < j - 1; ++k) data[j - k - 1] = div * ((w + k) * data[j - k - 2] + (j - k - w) * data[j - k - 1]);
    data[0] = div * (1.0 - w) * data[0];
  }
  ddata[0] = -data[0];
  for (int k = 1; k < n; ++k) ddata[k] = data[k - 1] - data[k];
  const double div = 1.0 / (n - 1);
  data[n - 1] = div * w * data[n - 2];
  for (int k = 1; k < n - 1; ++k) data[n - k - 1] = div * ((w + k) * data[n - k - 2] + (n - k - w) * data[n - k - 1]);
  data[0] = div * (1.0 - w) * data[0];
}

// Recursive mixed-radix FFT on a strided view. `tw` holds exp(sign 2πi j / root) for j < root, where root is
// the length of the whole transform (n divides it); `scratch` has room for root values.
void fft_rec(std::complex<double>* out, const std::complex<double>* in, std::size_t n, std::size_t stride,
             const std::complex<double>* tw, std::size_t root, std::complex<double>* scratch) {
  if (n == 1) {
    out[0] = in[0];
    return;
  }
  const std::size_t p = n % 2 == 0 ? 2 : n % 3 == 0 ? 3 : 5;
  const std::size_t m = n / p;
  for (std::size_t r = 0; r < p; ++r) fft_rec(out + r * m, in + r * stride, m, stride * p, tw, root, scratch);
  const std::size_t step = root / n;
  for (std::size_t k = 0; k < m; ++k)
    for (std::size_t q = 0; q < p; ++q) {
      std::complex<double> sum = 0.0;
      const std::size_t kk = k + q * m;
      for (std::size_t r = 0; r < p; ++r) sum += out[r * m + k] * tw[((r * kk) % n) * step];
      scratch[kk] = sum;
    }
  std::copy(scratch, scratch + n, out);
}

struct FftCache {
  std::size_t n = 0;
  std::vector<std::complex<double>> fwd, inv, out, scratch;
};

FftCache& fft_cache(std::size_t n) {
  thread_local std::vector<FftCache> caches;
  for (auto& c : caches)
    if (c.n == n) return c;
  FftCache c;
  c.n = n;
  c.fwd.resize(n);
  c.inv.resize(n);
  for (std::size_t j = 0; j < n; ++j) {
    const double a = 2.0 * kPi * double(j) / double(n);
    c.fwd[j] = {std::cos(a), -std::sin(a)};
    c.inv[j] = {std::cos(a), std::sin(a)};
  }
  c.out.resize(n);
  c.scratch.resize(n);
  caches.push_back(std::move(c));
  return caches.back();
}

}  // namespace

double ewald_beta(double cutoff, double tolerance) {
  // erfc(β rc) = tol: bisection on β.
  double lo = 0.0, hi = 10.0 / cutoff;
  for (int i = 0; i < 100; ++i) {
    const double mid = 0.5 * (lo + hi);
    (std::erfc(mid * cutoff) > tolerance ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

std::size_t fft_size(std::size_t n) {
  for (std::size_t k = std::max<std::size_t>(n, 1);; ++k) {
    std::size_t x = k;
    for (std::size_t p : {2, 3, 5})
      while (x % p == 0) x /= p;
    if (x == 1) return k;
  }
}

void fft_inplace(std::vector<std::complex<double>>& a, int sign) {
  const std::size_t n = a.size();
  if (n <= 1) return;
  FftCache& c = fft_cache(n);
  fft_rec(c.out.data(), a.data(), n, 1, sign < 0 ? c.fwd.data() : c.inv.data(), n, c.scratch.data());
  std::copy(c.out.begin(), c.out.end(), a.begin());
}

Pme::Pme(const std::array<double, 3>& box, double beta, int order, double spacing) : L_(box), beta_(beta), n_(order) {
  for (int d = 0; d < 3; ++d) K_[d] = fft_size(static_cast<std::size_t>(std::ceil(box[d] / spacing)));
  moduli();
}

void Pme::set_box(const std::array<double, 3>& box) { L_ = box; }

void Pme::moduli() {
  // |b(m)|^-2 terms: the spline values at integer points, Fourier-transformed.
  std::vector<double> data(n_), ddata(n_), bsp(static_cast<std::size_t>(n_) + 1, 0.0);
  bspline(0.0, n_, data.data(), ddata.data());
  for (int i = 0; i < n_; ++i) bsp[static_cast<std::size_t>(i) + 1] = data[static_cast<std::size_t>(i)];
  for (int d = 0; d < 3; ++d) {
    const std::size_t K = K_[d];
    bmod_[d].assign(K, 0.0);
    for (std::size_t m = 0; m < K; ++m) {
      double sc = 0.0, ss = 0.0;
      for (std::size_t j = 0; j < K && j <= static_cast<std::size_t>(n_); ++j) {
        const double arg = 2.0 * kPi * double(m * j) / double(K);
        sc += bsp[j] * std::cos(arg);
        ss += bsp[j] * std::sin(arg);
      }
      bmod_[d][m] = sc * sc + ss * ss;
    }
    // Zeros (odd order at the Nyquist frequency): average the neighbours.
    for (std::size_t m = 0; m < K; ++m)
      if (bmod_[d][m] < 1e-7) bmod_[d][m] = 0.5 * (bmod_[d][(m + K - 1) % K] + bmod_[d][(m + 1) % K]);
  }
}

double Pme::energy(const std::vector<double>& xyz, const std::vector<double>& q, std::vector<double>* forces) {
  const std::size_t N = q.size();
  const std::size_t K0 = K_[0], K1 = K_[1], K2 = K_[2];
  const std::size_t G = K0 * K1 * K2;
  grid_.assign(G, 0.0);
  const int n = n_;
  std::vector<double> th(3 * N * n), dth(3 * N * n);
  std::vector<int> base(3 * N);
  // Spread.
  for (std::size_t i = 0; i < N; ++i) {
    for (int d = 0; d < 3; ++d) {
      double s = xyz[3 * i + d] / L_[d];
      s -= std::floor(s);
      const double u = s * double(K_[d]);
      const double fl = std::floor(u);
      bspline(u - fl, n, &th[(3 * i + d) * n], &dth[(3 * i + d) * n]);
      base[3 * i + d] = static_cast<int>(fl) - n + 1;
    }
    for (int a = 0; a < n; ++a) {
      const std::size_t ia = static_cast<std::size_t>((base[3 * i] + a + 4 * int(K0)) % int(K0));
      const double wa = q[i] * th[(3 * i) * n + a];
      for (int b = 0; b < n; ++b) {
        const std::size_t ib = static_cast<std::size_t>((base[3 * i + 1] + b + 4 * int(K1)) % int(K1));
        const double wab = wa * th[(3 * i + 1) * n + b];
        for (int c = 0; c < n; ++c) {
          const std::size_t ic = static_cast<std::size_t>((base[3 * i + 2] + c + 4 * int(K2)) % int(K2));
          grid_[(ia * K1 + ib) * K2 + ic] += wab * th[(3 * i + 2) * n + c];
        }
      }
    }
  }
  // 3-D FFT along each axis.
  // 3-D FFT along each axis; the lines of one axis are independent, so they run in parallel.
  auto fft3 = [&](int sign) {
    parallel_for(K0 * K1, [&](std::size_t ab) {
      std::vector<std::complex<double>> line(K2);
      for (std::size_t c = 0; c < K2; ++c) line[c] = grid_[ab * K2 + c];
      fft_inplace(line, sign);
      for (std::size_t c = 0; c < K2; ++c) grid_[ab * K2 + c] = line[c];
    });
    parallel_for(K0 * K2, [&](std::size_t ac) {
      const std::size_t a = ac / K2, c = ac % K2;
      std::vector<std::complex<double>> line(K1);
      for (std::size_t b = 0; b < K1; ++b) line[b] = grid_[(a * K1 + b) * K2 + c];
      fft_inplace(line, sign);
      for (std::size_t b = 0; b < K1; ++b) grid_[(a * K1 + b) * K2 + c] = line[b];
    });
    parallel_for(K1 * K2, [&](std::size_t bc) {
      const std::size_t b = bc / K2, c = bc % K2;
      std::vector<std::complex<double>> line(K0);
      for (std::size_t a = 0; a < K0; ++a) line[a] = grid_[(a * K1 + b) * K2 + c];
      fft_inplace(line, sign);
      for (std::size_t a = 0; a < K0; ++a) grid_[(a * K1 + b) * K2 + c] = line[a];
    });
  };
  fft3(+1);
  // Multiply by B·C and accumulate the energy.
  const double V = L_[0] * L_[1] * L_[2];
  const double fac = kPi * kPi / (beta_ * beta_);
  double E = 0.0;
  for (std::size_t a = 0; a < K0; ++a) {
    const double mx = double(a <= K0 / 2 ? int(a) : int(a) - int(K0)) / L_[0];
    for (std::size_t b = 0; b < K1; ++b) {
      const double my = double(b <= K1 / 2 ? int(b) : int(b) - int(K1)) / L_[1];
      for (std::size_t c = 0; c < K2; ++c) {
        const std::size_t idx = (a * K1 + b) * K2 + c;
        if (a == 0 && b == 0 && c == 0) {
          grid_[idx] = 0.0;
          continue;
        }
        const double mz = double(c <= K2 / 2 ? int(c) : int(c) - int(K2)) / L_[2];
        const double m2 = mx * mx + my * my + mz * mz;
        const double denom = kPi * V * bmod_[0][a] * bmod_[1][b] * bmod_[2][c] * m2;
        const double eterm = kCoulomb * std::exp(-fac * m2) / denom;
        E += 0.5 * eterm * std::norm(grid_[idx]);
        grid_[idx] *= eterm;
      }
    }
  }
  if (!forces) return E;
  fft3(-1);  // back to real space: the convolution θ_rec ⋆ Q
  parallel_for(N, [&](std::size_t i) {
    double fx = 0.0, fy = 0.0, fz = 0.0;
    for (int a = 0; a < n; ++a) {
      const std::size_t ia = static_cast<std::size_t>((base[3 * i] + a + 4 * int(K0)) % int(K0));
      for (int b = 0; b < n; ++b) {
        const std::size_t ib = static_cast<std::size_t>((base[3 * i + 1] + b + 4 * int(K1)) % int(K1));
        for (int c = 0; c < n; ++c) {
          const std::size_t ic = static_cast<std::size_t>((base[3 * i + 2] + c + 4 * int(K2)) % int(K2));
          const double g = grid_[(ia * K1 + ib) * K2 + ic].real();
          const double ta = th[(3 * i) * n + a], tb = th[(3 * i + 1) * n + b], tc = th[(3 * i + 2) * n + c];
          fx += dth[(3 * i) * n + a] * tb * tc * g;
          fy += ta * dth[(3 * i + 1) * n + b] * tc * g;
          fz += ta * tb * dth[(3 * i + 2) * n + c] * g;
        }
      }
    }
    (*forces)[3 * i] -= q[i] * fx * double(K0) / L_[0];
    (*forces)[3 * i + 1] -= q[i] * fy * double(K1) / L_[1];
    (*forces)[3 * i + 2] -= q[i] * fz * double(K2) / L_[2];
  });
  return E;
}

double ewald_reciprocal_exact(const std::array<double, 3>& L, double beta, const std::vector<double>& xyz,
                              const std::vector<double>& q, int kmax, std::vector<double>* forces) {
  const std::size_t N = q.size();
  const double V = L[0] * L[1] * L[2];
  double E = 0.0;
  for (int nx = -kmax; nx <= kmax; ++nx)
    for (int ny = -kmax; ny <= kmax; ++ny)
      for (int nz = -kmax; nz <= kmax; ++nz) {
        if (nx == 0 && ny == 0 && nz == 0) continue;
        const double kx = 2 * kPi * nx / L[0], ky = 2 * kPi * ny / L[1], kz = 2 * kPi * nz / L[2];
        const double k2 = kx * kx + ky * ky + kz * kz;
        const double A = std::exp(-k2 / (4 * beta * beta)) / k2;
        double sr = 0.0, si = 0.0;
        for (std::size_t j = 0; j < N; ++j) {
          const double ph = kx * xyz[3 * j] + ky * xyz[3 * j + 1] + kz * xyz[3 * j + 2];
          sr += q[j] * std::cos(ph);
          si += q[j] * std::sin(ph);
        }
        E += 2 * kPi / V * kCoulomb * A * (sr * sr + si * si);
        if (forces)
          for (std::size_t j = 0; j < N; ++j) {
            const double ph = kx * xyz[3 * j] + ky * xyz[3 * j + 1] + kz * xyz[3 * j + 2];
            // −∂/∂r_j of A|S|²: 2 A q_j k (sin(ph) S_r − cos(ph) S_i)
            const double s = 4 * kPi / V * kCoulomb * A * q[j] * (std::sin(ph) * sr - std::cos(ph) * si);
            (*forces)[3 * j] += s * kx;
            (*forces)[3 * j + 1] += s * ky;
            (*forces)[3 * j + 2] += s * kz;
          }
      }
  return E;
}

double ewald_self(double beta, const std::vector<double>& q) {
  double s = 0.0;
  for (double x : q) s += x * x;
  return -kCoulomb * beta / std::sqrt(kPi) * s;
}

}  // namespace chemlarp
