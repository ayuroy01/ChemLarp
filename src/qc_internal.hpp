#pragma once
// Internal pieces shared by qc.cpp (integrals, SCF) and qc_dft.cpp (grid, functionals).

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace chemlarp::qc {

inline constexpr double kAngstromToBohr = 1.0 / 0.529177210903;  // CODATA 2018 Bohr radius
inline constexpr double kAuToDebye = 2.541746473;

struct Shell {
  std::array<double, 3> center{};  // bohr
  int L = 0;
  int atom = 0;                    // index in Molecule::atoms
  int first = 0;                   // first basis function
  std::vector<double> exps, coefs; // coefficients include the primitive normalisation for x^L
};

struct Molecule {
  std::vector<int> z;                         // atomic numbers
  std::vector<std::array<double, 3>> xyz;     // bohr
  int charge = 0;
};

struct Basis {
  std::vector<Shell> shells;
  std::vector<double> scale;                  // per function: makes each Cartesian component unit-normalised
  std::vector<int> atom_of;                   // per function
  std::vector<std::array<int, 3>> powers;     // per function: (lx, ly, lz)
  int n = 0;
  int max_l = 0;
};

// Cartesian components of angular momentum L in the usual order (xx, xy, xz, yy, yz, zz).
std::vector<std::array<int, 3>> cartesian_components(int L);

// Builds the basis for a molecule; false (with `error`) if an element is missing from the set.
bool build_basis(const Molecule& mol, const std::string& name, Basis& out, std::string& error);

// Symmetric N x N matrices, row-major.
using Mat = std::vector<double>;

// Eigenvalues (ascending) and eigenvectors (columns of `vectors`) of symmetric `a` (cyclic Jacobi).
void symmetric_eigen(const Mat& a, int n, std::vector<double>& values, Mat& vectors);

// Exchange-correlation on a molecular grid (qc_dft.cpp).
enum class Functional { SVWN5, PBE, B3LYP };

struct XcResult {
  double energy = 0.0;
  double electrons = 0.0;
  Mat va, vb;  // Kohn-Sham potential matrices for alpha and beta
};

class XcIntegrator {
 public:
  XcIntegrator(const Molecule& mol, const Basis& basis, Functional f, int grid_level);
  ~XcIntegrator();
  XcIntegrator(const XcIntegrator&) = delete;
  XcIntegrator& operator=(const XcIntegrator&) = delete;
  // Densities in the AO basis; with `restricted`, pb must equal pa.
  XcResult evaluate(const Mat& pa, const Mat& pb, bool restricted) const;
  std::size_t points() const;

 private:
  struct Impl;
  Impl* impl_;
};

// Fraction of exact exchange for a functional (0 for pure DFT).
double exact_exchange_fraction(Functional f);

}  // namespace chemlarp::qc
