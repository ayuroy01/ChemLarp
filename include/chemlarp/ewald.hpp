// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#pragma once
// Ewald summation for point charges in an orthorhombic periodic box, and its smooth particle-mesh form.
//
// Units: Å, elementary charges, kcal/mol (Coulomb constant kCoulomb). The Ewald split with parameter β:
//   real space    ke Σ_{i<j, r<rc} q_i q_j erfc(β r)/r                 (minimum image, not done here)
//   reciprocal    ke (2π/V) Σ_{k≠0} exp(−k²/4β²)/k² |S(k)|²,  S(k) = Σ q_j exp(i k·r_j)
//   self          −ke β/√π Σ q_i²
//   excluded      −ke q_i q_j erf(β r)/r for pairs whose direct interaction is left out (bonded atoms of a
//                 molecule): the reciprocal sum includes them, this removes them again.
// A neutral box is assumed (no background term).
//
// Smooth PME (Essmann, Perera, Berkowitz, Darden, Lee & Pedersen, J. Chem. Phys. 103 (1995) 8577): charges
// spread on a grid with cardinal B-splines of order n, the reciprocal sum by FFT, forces from the gradient of
// the splines. The FFT is a mixed-radix (2, 3, 5) Cooley–Tukey transform written here; grid sizes are rounded
// up to such numbers. ewald_reciprocal_exact() is the direct k-space sum, used to check PME.

#include <array>
#include <complex>
#include <cstddef>
#include <vector>

namespace chemlarp {

inline constexpr double kCoulomb = 332.0637;  // kcal Å / (mol e²)

// β such that erfc(β rc) = tolerance.
double ewald_beta(double cutoff, double tolerance);

// Smallest n >= min_n whose prime factors are 2, 3 and 5.
std::size_t fft_size(std::size_t min_n);

// In-place forward (sign -1) or inverse (sign +1, unnormalized) complex FFT of length n (factors 2, 3, 5).
void fft_inplace(std::vector<std::complex<double>>& a, int sign);

class Pme {
 public:
  // box: edge lengths (Å); grid sizes chosen so that the spacing is at most `spacing`.
  Pme(const std::array<double, 3>& box, double beta, int order = 6, double spacing = 1.0);

  void set_box(const std::array<double, 3>& box);  // keeps the grid sizes
  const std::array<std::size_t, 3>& grid() const noexcept { return K_; }

  // Reciprocal energy (kcal/mol); adds −∂E/∂r to forces (3N) when not null.
  double energy(const std::vector<double>& xyz, const std::vector<double>& q, std::vector<double>* forces);

 private:
  void moduli();
  std::array<double, 3> L_{};
  double beta_;
  int n_;
  std::array<std::size_t, 3> K_{};
  std::array<std::vector<double>, 3> bmod_;  // |b(m)|² per dimension
  std::vector<std::complex<double>> grid_;
};

// Direct reciprocal-space Ewald sum over |n_i| <= kmax; adds forces when not null.
double ewald_reciprocal_exact(const std::array<double, 3>& box, double beta, const std::vector<double>& xyz,
                              const std::vector<double>& q, int kmax, std::vector<double>* forces);

double ewald_self(double beta, const std::vector<double>& q);

}  // namespace chemlarp
