#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

const char* to_string(DiagCode code) noexcept {
  switch (code) {
    case DiagCode::EmptyInput: return "EmptyInput";
    case DiagCode::UnexpectedCharacter: return "UnexpectedCharacter";
    case DiagCode::UnknownElement: return "UnknownElement";
    case DiagCode::UnclosedBracketAtom: return "UnclosedBracketAtom";
    case DiagCode::InvalidBracketAtom: return "InvalidBracketAtom";
    case DiagCode::UnmatchedOpenParen: return "UnmatchedOpenParen";
    case DiagCode::UnmatchedCloseParen: return "UnmatchedCloseParen";
    case DiagCode::EmptyBranch: return "EmptyBranch";
    case DiagCode::BondWithoutAtom: return "BondWithoutAtom";
    case DiagCode::UnclosedRing: return "UnclosedRing";
    case DiagCode::RingBondConflict: return "RingBondConflict";
    case DiagCode::RingSelfBond: return "RingSelfBond";
    case DiagCode::DuplicateBond: return "DuplicateBond";
    case DiagCode::UnsupportedFeature: return "UnsupportedFeature";
    case DiagCode::ValenceExceeded: return "ValenceExceeded";
    case DiagCode::InvalidAromaticAtom: return "InvalidAromaticAtom";
    case DiagCode::KekulizationFailed: return "KekulizationFailed";
    case DiagCode::KekulizationLimit: return "KekulizationLimit";
    case DiagCode::ConflictingBondStereo: return "ConflictingBondStereo";
    case DiagCode::StereoIgnored: return "StereoIgnored";
    case DiagCode::CipUndetermined: return "CipUndetermined";
    case DiagCode::ThermoFormulaMatch: return "ThermoFormulaMatch";
    case DiagCode::MalformedReaction: return "MalformedReaction";
    case DiagCode::EmptyReactionSide: return "EmptyReactionSide";
    case DiagCode::InvalidCoefficient: return "InvalidCoefficient";
    case DiagCode::ElementImbalance: return "ElementImbalance";
    case DiagCode::ChargeImbalance: return "ChargeImbalance";
    case DiagCode::NoThermoData: return "NoThermoData";
    case DiagCode::NoJobackGroup: return "NoJobackGroup";
    case DiagCode::ThermoApproximation: return "ThermoApproximation";
    case DiagCode::TemperatureOutOfRange: return "TemperatureOutOfRange";
    case DiagCode::InvalidInput: return "InvalidInput";
    case DiagCode::NegativeAmount: return "NegativeAmount";
    case DiagCode::EmptySystem: return "EmptySystem";
    case DiagCode::ChargedSpeciesUnsupported: return "ChargedSpeciesUnsupported";
    case DiagCode::SingularSystem: return "SingularSystem";
    case DiagCode::NotConverged: return "NotConverged";
    case DiagCode::MechanismSyntax: return "MechanismSyntax";
    case DiagCode::UnknownSpecies: return "UnknownSpecies";
    case DiagCode::DuplicateSpecies: return "DuplicateSpecies";
    case DiagCode::InvalidRateParameter: return "InvalidRateParameter";
    case DiagCode::UnsupportedReactionType: return "UnsupportedReactionType";
    case DiagCode::MissingHeatCapacity: return "MissingHeatCapacity";
    case DiagCode::InvalidReactorSetup: return "InvalidReactorSetup";
    case DiagCode::StepSizeUnderflow: return "StepSizeUnderflow";
    case DiagCode::TooManySteps: return "TooManySteps";
    case DiagCode::DuplicateReaction: return "DuplicateReaction";
    case DiagCode::CompositionMismatch: return "CompositionMismatch";
    case DiagCode::NoForceFieldType: return "NoForceFieldType";
    case DiagCode::ImplicitHydrogens: return "ImplicitHydrogens";
    case DiagCode::EmbeddingFailed: return "EmbeddingFailed";
    case DiagCode::InconsistentBounds: return "InconsistentBounds";
    case DiagCode::MinimizationNotConverged: return "MinimizationNotConverged";
    case DiagCode::StereoLost: return "StereoLost";
    case DiagCode::NoCriticalData: return "NoCriticalData";
    case DiagCode::SupercriticalTemperature: return "SupercriticalTemperature";
    case DiagCode::FlashNotConverged: return "FlashNotConverged";
    case DiagCode::InvalidMixtureInput: return "InvalidMixtureInput";
    case DiagCode::NoInteractionParameter: return "NoInteractionParameter";
    case DiagCode::NoUnifacGroups: return "NoUnifacGroups";
    case DiagCode::NoVaporPressureData: return "NoVaporPressureData";
    case DiagCode::VleNotConverged: return "VleNotConverged";
    case DiagCode::NoPropertyData: return "NoPropertyData";
    case DiagCode::PropertyExtrapolated: return "PropertyExtrapolated";
    case DiagCode::TemperatureCross: return "TemperatureCross";
    case DiagCode::UnbalancedReaction: return "UnbalancedReaction";
    case DiagCode::FlowsheetSyntax: return "FlowsheetSyntax";
    case DiagCode::FlowsheetTopology: return "FlowsheetTopology";
    case DiagCode::RecycleNotConverged: return "RecycleNotConverged";
    case DiagCode::InchiMetalDisconnected: return "InchiMetalDisconnected";
    case DiagCode::InchiSymmetryLimit: return "InchiSymmetryLimit";
    case DiagCode::InchiUnsupported: return "InchiUnsupported";
    case DiagCode::UnknownFileFormat: return "UnknownFileFormat";
    case DiagCode::MalformedRecord: return "MalformedRecord";
    case DiagCode::UnsupportedFileFeature: return "UnsupportedFileFeature";
    case DiagCode::BondOrdersGuessed: return "BondOrdersGuessed";
    case DiagCode::FileContentIgnored: return "FileContentIgnored";
    case DiagCode::SmartsSyntax: return "SmartsSyntax";
    case DiagCode::SmartsUnsupported: return "SmartsUnsupported";
    case DiagCode::NamingUnsupported: return "NamingUnsupported";
  }
  return "Unknown";
}

const char* to_string(Severity severity) noexcept {
  return severity == Severity::Error ? "error" : "warning";
}

}  // namespace chemlarp
