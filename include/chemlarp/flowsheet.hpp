#pragma once
// Steady-state flowsheets: unit operations (process.hpp) connected by named streams.
//
// Solution method: sequential modular. Recycle loops are broken at tear streams, chosen from the back edges
// of a depth-first search over the unit graph (or named by the user); the units are then computed in a
// topological order, and the tear streams are iterated to a fixed point by Wegstein's method (bounded
// acceleration factor, q in [-5, 0], after two direct-substitution steps) or by direct substitution.
// The iteration variables of a tear stream are its component flows and enthalpy flow; its pressure is
// taken from the computed value, and its temperature and phase split come from a PH flash.
// References: Wegstein, Commun. ACM 1 (1958) 9; Westerberg, Hutchison, Motard & Winter, Process
// Flowsheeting (1979); Biegler, Grossmann & Westerberg, Systematic Methods of Chemical Process Design (1997).
//
// Text format, one statement per line; '#' at the start of a line or after a space starts a comment.
// Quantities are SI (K, Pa, W, mol/s) unless a unit is attached: T=80C, P=1.5bar|kPa|MPa|atm|Pa, Q=-20kW|MW|W.
//   components <SMILES> <SMILES> ...
//   model ideal | unifac                       (liquid activity model; default ideal)
//   feed  <stream> T=<T> P=<P> flows <SMILES>=<mol/s> ...
//   guess <stream> T=<T> P=<P> flows <SMILES>=<mol/s> ...   (initial value; makes the stream a tear)
//   tear  <stream> ...                         (optional: choose tear streams)
//   method wegstein | direct
//   tolerance <relative>                       (default 1e-7)
//   max-iterations <n>                         (default 100)
//   mixer      <name> <in> <in> ... -> <out> [P=]
//   splitter   <name> <in> -> <out> <out> ... fractions=<f1>,<f2>,...
//   heater     <name> <in> -> <out> T= | Q= [P=]           (also: cooler)
//   valve      <name> <in> -> <out> P=
//   flash      <name> <in> -> <vapour> <liquid> T= P= | P=  (P only: adiabatic)
//   separator  <name> <in> -> <first> <second> [T=] [P=] split <SMILES>=<fraction to first> ...
//   pump       <name> <in> -> <out> P= [eff=]               (efficiency, default 1)
//   compressor <name> <in> -> <out> P= [eff=]               (also: expander, turbine)
//   exchanger  <name> <hot in> <cold in> -> <hot out> <cold out> Thot= | Tcold= | Q=
//   reactor    <name> <in> -> <out> key=<SMILES> X=<conversion> [T=] [P=] nu <SMILES>=<coefficient> ...
//              (no T: adiabatic)
// A stream has one source (a feed or a unit outlet) and at most one destination; outlets without a
// destination are products.

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/process.hpp"

namespace chemlarp {

enum class UnitKind { Mixer, Splitter, Heater, Valve, Flash, Separator, Pump, Compressor, Exchanger, Reactor };
const char* to_string(UnitKind kind) noexcept;

struct FeedSpec {
  std::string name;
  double T_K = 0.0, P_Pa = 0.0;
  std::vector<double> n;  // mol/s, one per component
};

struct UnitSpec {
  std::string name;
  UnitKind kind = UnitKind::Mixer;
  std::vector<std::string> inlets, outlets;
  // Parameters; 0 (or empty) when not given.
  double T_K = 0.0, P_Pa = 0.0, efficiency = 0.0, conversion = 0.0, hot_T_K = 0.0, cold_T_K = 0.0;
  double duty_W = 0.0;
  bool has_duty = false;
  std::vector<double> fractions;      // splitter
  std::vector<double> per_component;  // separator split fractions, or reactor coefficients
  std::size_t key = 0;                // reactor key component
  int line = 0;
};

enum class TearMethod { Wegstein, Direct };

struct FlowsheetSpec {
  std::vector<std::string> components;
  LiquidModel model = LiquidModel::Ideal;
  std::vector<FeedSpec> feeds;
  std::vector<UnitSpec> units;
  std::vector<std::string> tears;  // user-chosen tear streams
  std::vector<FeedSpec> guesses;   // initial values for tear streams
  TearMethod method = TearMethod::Wegstein;
  int max_iterations = 100;
  double tolerance = 1e-7;
};

std::optional<FlowsheetSpec> parse_flowsheet(std::string_view text, Diagnostics& diag);

struct UnitReport {
  std::string name;
  UnitKind kind = UnitKind::Mixer;
  double duty_W = 0.0;  // heat added from outside (zero for an exchanger, whose duty is internal)
  double work_W = 0.0;  // shaft work added
  std::vector<std::pair<std::string, double>> values;  // other results (exchanger duty, UA, extent, ...)
};

struct NamedStream {
  std::string name, source, destination;  // source: unit name or "feed"; destination: unit name or "" (product)
  Stream stream;
};

struct FlowsheetResult {
  ProcessComponents components;
  bool converged = false;
  int iterations = 0;
  std::vector<std::string> tears, order;  // tear streams; units in calculation order
  std::vector<double> residuals;          // largest relative tear-stream change, per iteration
  std::vector<NamedStream> streams;       // feeds first, then unit outlets in unit order
  std::vector<UnitReport> units;
  double mass_in_kgs = 0.0, mass_out_kgs = 0.0;
  double energy_imbalance_W = 0.0;  // feed enthalpy + duties + work − product enthalpy
};

// Solves the flowsheet. A result without convergence (RecycleNotConverged) is still returned, with the
// last iterate; errors in a unit or in the topology return nullopt.
std::optional<FlowsheetResult> solve_flowsheet(const FlowsheetSpec& spec, Diagnostics& diag);

// Parses and solves a flowsheet text and returns the result as JSON:
// {"ok","converged","iterations","method","tears","order","residuals","components":[{smiles,name,molar_mass}],
//  "streams":[{name,from,to,T_K,P_Pa,vapor_fraction,H_W,mass_kgs,flows,x,y}],
//  "units":[{name,type,duty_W,work_W,values:{...}}],"balance":{mass_in_kgs,mass_out_kgs,energy_imbalance_W},
//  "diagnostics"}. "ok" is true when the flowsheet solved and converged without errors.
std::string flowsheet_json(std::string_view text);

}  // namespace chemlarp
