#pragma once
// Rigorous equilibrium-stage distillation column.
//
// Stages are numbered from the top: stage 1 is the condenser, stage N the reboiler (a kettle reboiler,
// itself an equilibrium stage). On every stage the MESH equations hold with the stream model of
// process.hpp: component Material balances, Equilibrium y_i = K_i x_i (modified Raoult, ideal or UNIFAC
// liquid), Summation (implicit in the component-flow variables), and Heat balances without heat loss.
// The condenser is total (liquid distillate at its bubble point, the reflux of the same composition) or
// partial (vapour distillate in equilibrium with the reflux). The two degrees of freedom are the reflux
// ratio L/D and the distillate or the bottoms flow; the condenser and reboiler duties follow.
//
// Method: Naphtali–Sandholm. The variables of each stage are its temperature and the component flows of
// the liquid and vapour leaving it (for a total condenser: the reflux and the distillate); all N(2C+1)
// equations are solved together by Newton's method with a block-tridiagonal Jacobian (finite differences)
// and a backtracking line search. The starting profile comes from constant molar overflow and a few
// bubble-point (Wang–Henke) iterations. References: Naphtali & Sandholm, AIChE J. 17 (1971) 148;
// Wang & Henke, Hydrocarbon Process. 45(8) (1966) 155; Seader, Henley & Roper, Separation Process
// Principles, ch. 10.
//
// Every component needs a vapour-pressure equation (no non-condensables), and the liquid is one phase.

#include <optional>
#include <vector>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/process.hpp"

namespace chemlarp {

enum class Condenser { Total, Partial };

struct ColumnFeed {
  std::size_t stage = 2;  // 1-based: 1 is the condenser, N the reboiler
  Stream stream;
};

struct ColumnSpec {
  std::size_t stages = 10;  // N, including the condenser and the reboiler
  std::vector<ColumnFeed> feeds;
  double P_top_Pa = 0.0, P_bottom_Pa = 0.0;  // linear profile; P_bottom <= 0 means the top pressure throughout
  Condenser condenser = Condenser::Total;
  double reflux_ratio = 0.0;    // L/D, > 0
  double distillate_rate = 0.0;  // mol/s; give this or bottoms_rate
  double bottoms_rate = 0.0;
  int max_iterations = 60;
};

struct ColumnStage {
  double T_K = 0.0, P_Pa = 0.0;
  double L = 0.0, V = 0.0;   // liquid and vapour flows leaving the stage (for a total condenser, L is the reflux)
  std::vector<double> x, y;  // mole fractions
};

struct ColumnResult {
  Stream distillate, bottoms;
  double condenser_duty_W = 0.0;  // heat added (negative: removed)
  double reboiler_duty_W = 0.0;
  std::vector<ColumnStage> stages;
  int iterations = 0;
  double residual = 0.0;  // largest scaled MESH residual at the solution
};

std::optional<ColumnResult> distillation_column(const ProcessComponents& c, const ColumnSpec& spec, Diagnostics& diag);

}  // namespace chemlarp
