#pragma once
// The error vector. Every stage (parser, valence validator, stoichiometry checker,
// thermo binder, equilibrium solver) appends structured diagnostics here instead
// of throwing, so a single call reports *all* problems with exact locations.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace omnichem {

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
  CipUndetermined,        // warning: CIP rules 1-2 cannot rank the substituents (rules 3-5 not implemented)

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

}  // namespace omnichem
