#pragma once
// The error vector. Every stage (parser, valence validator, stoichiometry checker,
// thermo binder, equilibrium solver) appends structured diagnostics here instead
// of throwing, so a single call reports *all* problems with exact locations.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace chemlarp {

enum class DiagCode : std::uint16_t {
  // 1xx — SMILES syntax
  EmptyInput = 100,
  UnexpectedCharacter,
  UnknownElement,
  UnclosedBracketAtom,
  InvalidBracketAtom,
  UnmatchedOpenParen,
  UnmatchedCloseParen,
  EmptyBranch,
  BondWithoutAtom,
  UnclosedRing,
  RingBondConflict,
  RingSelfBond,
  DuplicateBond,
  UnsupportedFeature,

  // 2xx — chemistry / valence
  ValenceExceeded = 200,
  InvalidAromaticAtom,
  KekulizationFailed,
  KekulizationLimit,
  ConflictingBondStereo,  // two '/' '\' marks on one double-bond end put both neighbours on the same side
  StereoIgnored,          // warning: stereo mark that cannot be honoured (not stereogenic, small ring, ...)
  CipUndetermined,        // warning: CIP rules 1-5 cannot rank the substituents (not stereogenic, or digraph too large)

  // 3xx — stoichiometry
  MalformedReaction = 300,
  EmptyReactionSide,
  InvalidCoefficient,
  ElementImbalance,
  ChargeImbalance,

  // 4xx — thermodynamics / equilibrium
  NoThermoData = 400,
  ThermoFormulaMatch,  // warning: database entry matched by formula only, not by canonical structure
  NoJobackGroup,
  ThermoApproximation,
  TemperatureOutOfRange,
  InvalidInput,
  NegativeAmount,
  EmptySystem,
  ChargedSpeciesUnsupported,
  SingularSystem,
  NotConverged,

  // 5xx — kinetics / ODE integration
  MechanismSyntax = 500,
  UnknownSpecies,
  DuplicateSpecies,
  InvalidRateParameter,
  UnsupportedReactionType,
  MissingHeatCapacity,
  InvalidReactorSetup,
  StepSizeUnderflow,
  TooManySteps,
  DuplicateReaction,    // warning: same reaction twice without DUPLICATE (or DUPLICATE without a partner)
  CompositionMismatch,  // SMILES structure and thermo-card composition disagree

  // 6xx — 3D structure / force fields
  NoForceFieldType = 600,  // atom has no UFF type (e.g. metals)
  ImplicitHydrogens,       // 3D work needs explicit hydrogens (add_explicit_hydrogens)
  EmbeddingFailed,         // distance geometry found no acceptable coordinates
  InconsistentBounds,      // triangle smoothing produced lower > upper
  MinimizationNotConverged,  // warning
  StereoLost,              // 3D structure does not reproduce the requested stereo

  // 7xx — phase equilibrium / process
  NoCriticalData = 700,    // component has no critical constants (Tc, Pc, omega)
  SupercriticalTemperature,  // saturation pressure asked above the critical temperature
  FlashNotConverged,       // warning: successive substitution did not converge
  InvalidMixtureInput,     // mole fractions, temperature or model parameters out of range
  NoInteractionParameter,  // activity model: no published parameter for a pair of groups or components
  NoUnifacGroups,          // a compound cannot be split into supported UNIFAC subgroups
  NoVaporPressureData,     // component has no vapour-pressure correlation
  VleNotConverged,         // bubble/dew point iteration did not converge
  NoPropertyData,          // the pure-component table lacks the property (or the compound)
  PropertyExtrapolated,    // warning: a correlation was used outside its stated range
  TemperatureCross,        // heat exchanger: the hot stream would be colder than the cold stream
  UnbalancedReaction,      // reactor stoichiometry does not conserve elements or charge
  FlowsheetSyntax,         // a flowsheet line that cannot be read
  FlowsheetTopology,       // streams not connected correctly (missing source, two sources, unused unit)
  RecycleNotConverged,     // the tear streams did not converge within the iteration limit

  // 8xx — identifiers (InChI)
  InchiMetalDisconnected = 800,  // warning: bonds to metals were cut, as standard InChI does
  InchiSymmetryLimit,            // warning: symmetry search truncated; the numbering may not be canonical
  InchiUnsupported,              // warning: a feature the generator does not encode (e.g. isotopic stereo)

  // 9xx — structure files (MOL/SDF, PDB, CML)
  UnknownFileFormat = 900,  // content is not a recognised structure file
  MalformedRecord,          // a line or element that cannot be read
  UnsupportedFileFeature,   // query atoms/bonds, atom lists, R-groups, ...
  BondOrdersGuessed,        // warning: the file has no bond orders for some atoms (PDB)
  FileContentIgnored,       // warning: content read but not used (S-groups, extra models, alternate locations)

  // 10xx — queries (SMARTS)
  SmartsSyntax = 1000,      // malformed SMARTS
  SmartsUnsupported,        // valid SMARTS feature this engine does not implement

  // 11xx — nomenclature
  NamingUnsupported = 1100,  // structure needs a nomenclature rule the namer does not implement
};

enum class Severity : std::uint8_t { Warning, Error };

const char* to_string(DiagCode code) noexcept;
const char* to_string(Severity severity) noexcept;

struct Diagnostic {
  static constexpr std::size_t npos = static_cast<std::size_t>(-1);

  DiagCode code{};
  Severity severity = Severity::Error;
  std::size_t position = npos;  // character offset in the input string
  std::int32_t atom = -1;       // global atom index in ChemicalSystemState
  std::int32_t species = -1;    // species index
  std::uint8_t element = 0;     // atomic number for stoichiometry diagnostics
  double delta = 0.0;           // imbalance (products - reactants) for ElementImbalance / ChargeImbalance
  std::string message;

  Diagnostic& at(std::size_t pos) { position = pos; return *this; }
  Diagnostic& on_atom(std::int64_t a) { atom = static_cast<std::int32_t>(a); return *this; }
  Diagnostic& on_species(std::int64_t s) { species = static_cast<std::int32_t>(s); return *this; }
  Diagnostic& with_element(std::uint8_t z, double d) { element = z; delta = d; return *this; }
};

class Diagnostics {
 public:
  Diagnostic& error(DiagCode code, std::string message) { return push(code, Severity::Error, std::move(message)); }
  Diagnostic& warn(DiagCode code, std::string message) { return push(code, Severity::Warning, std::move(message)); }

  bool has_errors() const noexcept { return error_count_ > 0; }
  std::size_t error_count() const noexcept { return error_count_; }
  std::size_t size() const noexcept { return items_.size(); }
  bool empty() const noexcept { return items_.empty(); }
  bool contains(DiagCode code) const noexcept {
    for (const auto& d : items_) if (d.code == code) return true;
    return false;
  }
  std::size_t count(DiagCode code) const noexcept {
    std::size_t n = 0;
    for (const auto& d : items_) n += (d.code == code);
    return n;
  }
  const std::vector<Diagnostic>& items() const noexcept { return items_; }
  void clear() noexcept { items_.clear(); error_count_ = 0; }

 private:
  Diagnostic& push(DiagCode code, Severity sev, std::string message) {
    items_.push_back(Diagnostic{});
    Diagnostic& d = items_.back();
    d.code = code;
    d.severity = sev;
    d.message = std::move(message);
    if (sev == Severity::Error) ++error_count_;
    return d;
  }
  std::vector<Diagnostic> items_;
  std::size_t error_count_ = 0;
};

}  // namespace chemlarp
