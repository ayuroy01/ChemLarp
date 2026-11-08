#include "omnichem/reaction.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "omnichem/elements.hpp"
#include "omnichem/smiles.hpp"
#include "omnichem/thermo.hpp"

namespace omnichem {
namespace {

struct Token {
  std::string_view text;
  std::size_t pos;
};

std::vector<Token> tokenize(std::string_view s) {
  std::vector<Token> out;
  std::size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
    if (i >= s.size()) break;
    const std::size_t start = i;
    while (i < s.size() && s[i] != ' ' && s[i] != '\t') ++i;
    out.push_back({s.substr(start, i - start), start});
  }
  return out;
}

bool looks_numeric(std::string_view t) {
  if (t.empty()) return false;
  bool digit = false;
  for (char c : t) {
    if (c >= '0' && c <= '9') digit = true;
    else if (c != '.' && c != 'e' && c != 'E' && c != '-' && c != '+') return false;
  }
  return digit && ((t[0] >= '0' && t[0] <= '9') || t[0] == '.');
}

}  // namespace

std::optional<Reaction> parse_reaction(ChemicalSystemState& state, std::string_view eq, Diagnostics& diag) {
  const auto tokens = tokenize(eq);
  Reaction rxn;
  bool on_products = false;
  bool arrow_seen = false;
  bool ok = true;
  bool expect_term = true;  // after start, '+' or arrow
  double coef = 1.0;
  bool have_coef = false;
  std::size_t coef_pos = 0;

  for (std::size_t t = 0; t < tokens.size(); ++t) {
    const Token& tk = tokens[t];
    if (tk.text == "->" || tk.text == ">>") {
      if (arrow_seen) {
        diag.error(DiagCode::MalformedReaction, "more than one reaction arrow").at(tk.pos);
        ok = false;
      }
      if (have_coef) {
        diag.error(DiagCode::MalformedReaction, "coefficient is not followed by a species").at(coef_pos);
        ok = false;
        have_coef = false;
      }
      if (t > 0 && tokens[t - 1].text == "+") {
        diag.error(DiagCode::MalformedReaction, "'+' directly before the arrow").at(tokens[t - 1].pos);
        ok = false;
      }
      arrow_seen = true;
      on_products = true;
      expect_term = true;
      continue;
    }
    if (tk.text == "+") {
      if (expect_term) {
        diag.error(DiagCode::MalformedReaction, "'+' without a preceding species").at(tk.pos);
        ok = false;
      }
      if (have_coef) {
        diag.error(DiagCode::MalformedReaction, "coefficient is not followed by a species").at(coef_pos);
        ok = false;
        have_coef = false;
      }
      expect_term = true;
      continue;
    }
    if (!expect_term) {
      diag.error(DiagCode::MalformedReaction, "missing '+' between species").at(tk.pos);
      ok = false;
    }
    if (looks_numeric(tk.text) && !have_coef) {
      char* endp = nullptr;
      std::string num(tk.text);
      const double v = std::strtod(num.c_str(), &endp);
      if (endp != num.c_str() + num.size() || !std::isfinite(v) || v <= 0.0) {
        diag.error(DiagCode::InvalidCoefficient,
                   "stoichiometric coefficient '" + num + "' must be a positive finite number")
            .at(tk.pos);
        ok = false;
        coef = 1.0;
      } else {
        coef = v;
      }
      have_coef = true;
      coef_pos = tk.pos;
      continue;
    }
    const std::int32_t sp = parse_smiles(state, tk.text, diag, tk.pos);
    if (sp < 0) {
      ok = false;
    } else {
      ReactionTerm term{have_coef ? coef : 1.0, static_cast<std::uint32_t>(sp), tk.pos};
      (on_products ? rxn.products : rxn.reactants).push_back(term);
    }
    coef = 1.0;
    have_coef = false;
    expect_term = false;
  }

  if (have_coef) {
    diag.error(DiagCode::MalformedReaction, "coefficient is not followed by a species").at(coef_pos);
    ok = false;
  }
  if (!arrow_seen) {
    diag.error(DiagCode::MalformedReaction, "reaction has no '->' or '>>' arrow").at(0);
    return std::nullopt;
  }
  if (expect_term && on_products && !tokens.empty() && tokens.back().text == "+") {
    diag.error(DiagCode::MalformedReaction, "trailing '+'").at(tokens.back().pos);
    ok = false;
  }
  if (rxn.reactants.empty() && ok) {
    diag.error(DiagCode::EmptyReactionSide, "reaction has no reactants").at(0);
    ok = false;
  }
  if (rxn.products.empty() && ok) {
    diag.error(DiagCode::EmptyReactionSide, "reaction has no products").at(eq.size());
    ok = false;
  }
  if (!ok) return std::nullopt;
  return rxn;
}

bool check_balance(const ChemicalSystemState& state, const Reaction& rxn, Diagnostics& diag, double tolerance) {
  std::array<double, kElementSlots> delta{};
  double charge = 0.0;
  auto accumulate = [&](const std::vector<ReactionTerm>& side, double sign) {
    for (const auto& t : side) {
      const Species& s = state.species[t.species];
      for (std::size_t z = 1; z < kElementSlots; ++z) delta[z] += sign * t.coefficient * s.composition[z];
      charge += sign * t.coefficient * s.net_charge;
    }
  };
  accumulate(rxn.reactants, -1.0);
  accumulate(rxn.products, +1.0);

  bool balanced = true;
  for (std::size_t z = 1; z < kElementSlots; ++z) {
    if (std::fabs(delta[z]) <= tolerance) continue;
    balanced = false;
    const std::string sym(element_symbol(static_cast<std::uint8_t>(z)));
    char buf[64];
    std::snprintf(buf, sizeof buf, "%+g", delta[z]);
    diag.error(DiagCode::ElementImbalance,
               sym + " is not conserved: products - reactants = " + buf + " atoms")
        .with_element(static_cast<std::uint8_t>(z), delta[z]);
  }
  if (std::fabs(charge) > tolerance) {
    balanced = false;
    char buf[64];
    std::snprintf(buf, sizeof buf, "%+g", charge);
    diag.error(DiagCode::ChargeImbalance, std::string("charge is not conserved: products - reactants = ") + buf)
        .with_element(0, charge);
  }
  return balanced;
}

std::optional<ReactionThermo> reaction_thermo(const ChemicalSystemState& state, const Reaction& rxn,
                                              double T, Diagnostics& diag) {
  ReactionThermo out;
  double dh_RT = 0.0, ds_R = 0.0;
  bool ok = true;
  auto accumulate = [&](const std::vector<ReactionTerm>& side, double sign) {
    for (const auto& t : side) {
      auto th = resolve_thermo(state, t.species, T, diag);
      if (!th) { ok = false; continue; }
      dh_RT += sign * t.coefficient * th->h_RT(T);
      ds_R += sign * t.coefficient * th->s_R(T);
    }
  };
  accumulate(rxn.reactants, -1.0);
  accumulate(rxn.products, +1.0);
  if (!ok) return std::nullopt;
  out.dG_RT = dh_RT - ds_R;
  out.dH_J = dh_RT * kGasConstant * T;
  out.dS_J = ds_R * kGasConstant;
  out.K = std::exp(-out.dG_RT);
  return out;
}

}  // namespace omnichem
