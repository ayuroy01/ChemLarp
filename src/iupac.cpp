// Structure -> IUPAC name (2013 Recommendations, P-rules cited inline).
//
// Outline: a heavy-atom view of the species is perceived into characteristic
// groups; the namer then builds names recursively. One routine names a
// "fragment" -- the whole molecule, or a substituent reached through one bond --
// by enumerating parent candidates (ring systems with fixed numbering, carbon
// chains), choosing a numbering for each by the lowest-locant rules of
// P-31.1.4, ranking the candidates by P-44/P-45 and assembling prefixes,
// parent hydride, unsaturation and suffix. Substituent prefixes are named by
// the same routine, so nesting and enclosing marks follow automatically.

#include "chemlarp/iupac.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "chemlarp/cip.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/smarts.hpp"
#include "chemlarp/smiles.hpp"

namespace chemlarp {
namespace {

// ---- small helpers -----------------------------------------------------------------------

using Key = std::int64_t;  // stereo element: atom (a) or bond (flag | a << 21 | b)
Key atom_key(int a) { return a; }
Key bond_key(int a, int b) {
  if (a > b) std::swap(a, b);
  return (Key{1} << 42) | (Key(a) << 21) | Key(b);
}

// Sort key for a locant: 1 < 2 < 3a < 4 < 10 < N.
double locant_key(const std::string& l) {
  if (l.empty()) return 1e9;
  if (l[0] == 'N') return 1e6 + static_cast<double>(l.size());
  double v = 0;
  std::size_t i = 0;
  while (i < l.size() && l[i] >= '0' && l[i] <= '9') v = v * 10 + (l[i++] - '0');
  if (i < l.size() && l[i] >= 'a' && l[i] <= 'z') v += 0.01 * (l[i] - 'a' + 1);
  return v;
}

std::vector<double> sorted_keys(const std::vector<std::string>& l) {
  std::vector<double> k;
  k.reserve(l.size());
  for (const auto& s : l) k.push_back(locant_key(s));
  std::sort(k.begin(), k.end());
  return k;
}

// First point of difference (P-31.1.4).
int compare_keys(const std::vector<double>& a, const std::vector<double>& b) {
  for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i) {
    if (a[i] < b[i]) return -1;
    if (a[i] > b[i]) return 1;
  }
  return 0;
}

std::string join_locants(std::vector<std::string> l) {
  std::sort(l.begin(), l.end(), [](const std::string& x, const std::string& y) { return locant_key(x) < locant_key(y); });
  std::stable_partition(l.begin(), l.end(), [](const std::string& x) { return !x.empty() && x[0] == 'N'; });
  std::string s;
  for (std::size_t i = 0; i < l.size(); ++i) s += (i ? "," : "") + l[i];
  return s;
}

// P-21.2.1: stems of the saturated unbranched hydrocarbons (1-99 carbons).
std::string alkane_stem(int n) {
  static const char* small[] = {"",         "meth",     "eth",      "prop",     "but",     "pent",    "hex",
                                "hept",     "oct",      "non",      "dec",      "undec",   "dodec",   "tridec",
                                "tetradec", "pentadec", "hexadec",  "heptadec", "octadec", "nonadec", "icos"};
  if (n <= 0 || n >= 100) return "";
  if (n <= 20) return small[n];
  static const char* units[] = {"", "hen", "do", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona"};
  static const char* tens[] = {"", "", "cos", "triacont", "tetracont", "pentacont", "hexacont", "heptacont", "octacont", "nonacont"};
  std::string u = units[n % 10];
  std::string t = tens[n / 10];
  if (!u.empty() && u.back() == 'a' && !t.empty() && (t[0] == 'a' || t[0] == 'o')) u.pop_back();
  return u + t;
}

std::string multiplier(int n, bool compound) {
  static const char* simple[] = {"",          "",          "di",        "tri",        "tetra",     "penta",    "hexa",
                                 "hepta",     "octa",      "nona",      "deca",       "undeca",    "dodeca",   "trideca",
                                 "tetradeca", "pentadeca", "hexadeca",  "heptadeca",  "octadeca",  "nonadeca", "icosa"};
  static const char* comp[] = {"", "", "bis", "tris", "tetrakis", "pentakis", "hexakis", "heptakis", "octakis", "nonakis", "decakis"};
  if (compound) return n <= 10 ? comp[n] : "";
  return n <= 20 ? simple[n] : "";
}

bool starts_with_vowel(const std::string& s) { return !s.empty() && std::string("aeiouy").find(s[0]) != std::string::npos; }

// Hydride + suffix, eliding a final 'e' of the hydride before a vowel (P-16.7.1).
std::string attach(std::string base, const std::string& locants, const std::string& word) {
  if (!base.empty() && base.back() == 'e' && starts_with_vowel(word)) base.pop_back();
  return base + (locants.empty() ? "" : "-" + locants + "-") + word;
}

// Enclosing marks nested ( [ { (P-16.5.4).
// A bracket holding only locants and single letters ("[3.2.1]", "[1,5-a]", "[b,d]", "[1,2,4]") is a
// fusion or von Baeyer descriptor, not an enclosing mark.
bool descriptor_bracket(const std::string& s, std::size_t open) {
  const std::size_t close = s.find(']', open);
  if (close == std::string::npos) return false;
  int run = 0;
  for (std::size_t i = open + 1; i < close; ++i) {
    const char c = s[i];
    if (c >= 'a' && c <= 'z') {
      if (++run > 1) return false;
    } else if ((c >= '0' && c <= '9') || c == ',' || c == '.' || c == '-') {
      run = 0;
    } else {
      return false;
    }
  }
  return true;
}

std::string enclose(const std::string& s) {
  int depth = 0, cur = 0;
  for (std::size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (c == '[' && descriptor_bracket(s, i)) {
      i = s.find(']', i);
      continue;
    }
    if (c == '(' || c == '[' || c == '{') depth = std::max(depth, ++cur);
    if (c == ')' || c == ']' || c == '}') --cur;
  }
  static const char* open = "([{";
  static const char* close = ")]}";
  return open[depth % 3] + s + close[depth % 3];
}

bool descriptor_bracket(const std::string& s, std::size_t open);
bool is_enclosed(const std::string& s) { return !s.empty() && (s[0] == '(' || s[0] == '{' || (s[0] == '[' && !descriptor_bracket(s, 0))); }

// Alphanumerical key (P-14.5.2): italic letters, locants, stereodescriptors and enclosing marks ignored.
std::string alpha_key(const std::string& s) {
  std::string k;
  for (char c : s)
    if (c >= 'a' && c <= 'z') k += c;
  return k;
}

bool ends_with(const std::string& s, const std::string& t) { return s.size() >= t.size() && s.compare(s.size() - t.size(), t.size(), t) == 0; }

bool has_digit(const std::string& s) {
  return std::any_of(s.begin(), s.end(), [](char c) { return c >= '0' && c <= '9'; });
}

// Joins name pieces: a hyphen goes between a letter or closing mark and a following locant.
void append_piece(std::string& out, const std::string& piece) {
  if (!out.empty() && !piece.empty()) {
    const char a = out.back(), b = piece[0];
    const bool b_locant = (b >= '0' && b <= '9') || b == 'N';
    const bool a_end = (a >= 'a' && a <= 'z') || a == ')' || a == ']' || a == '}';
    if (b_locant && a_end) out += '-';
  }
  out += piece;
}

// ---- ring templates with fixed IUPAC numbering ----------------------------------------------
//
// Two kinds. "Exact" templates are retained names of fully saturated rings
// (piperidine, oxolane, ...) and match only that form. "Mancude" templates are
// ring skeletons (elements and connectivity only) of parents with the maximum
// number of non-cumulative double bonds; the actual saturation is expressed by
// indicated hydrogen and hydro prefixes (P-31.1.4.2.4).

struct Template {
  const char* name;     // parent hydride (mancude: without indicated hydrogen)
  const char* pattern;  // atoms in locant order
  std::vector<std::string> locants;
  bool mancude;
};

const std::vector<std::string> kL3 = {"1", "2", "3"};
const std::vector<std::string> kL4 = {"1", "2", "3", "4"};
const std::vector<std::string> kL5 = {"1", "2", "3", "4", "5"};
const std::vector<std::string> kL6 = {"1", "2", "3", "4", "5", "6"};
const std::vector<std::string> kL7 = {"1", "2", "3", "4", "5", "6", "7"};
const std::vector<std::string> kL8 = {"1", "2", "3", "4", "5", "6", "7", "8"};
const std::vector<std::string> kNaph = {"1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a"};
const std::vector<std::string> kInd = {"1", "2", "3", "3a", "4", "5", "6", "7", "7a"};
const std::vector<std::string> kAzul = {"1", "2", "3", "3a", "4", "5", "6", "7", "8", "8a"};
const std::vector<std::string> kFluor = {"1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "9a"};
const std::vector<std::string> kAnth = {"1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"};
const std::vector<std::string> kPhen = {"1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "10", "10a"};
const std::vector<std::string> kPyr = {"1", "2", "3", "4", "4a", "5", "6", "7", "7a"};
const std::vector<std::string> kBzep = {"1", "2", "3", "4", "5", "5a", "6", "7", "8", "9", "9a"};
const std::vector<std::string> kDbf = {"1", "2", "3", "4", "4a", "5", "5a", "6", "7", "8", "9", "9a", "9b"};
const std::vector<std::string> kPhz = {"1", "2", "3", "4", "4a", "5", "5a", "6", "7", "8", "9", "9a", "10", "10a"};

const std::vector<Template>& templates() {
  static const std::vector<Template> t = {
      // retained / Hantzsch-Widman names of saturated heteromonocycles (exact form)
      {"piperidine", "[NX3]1-[C]-[C]-[C]-[C]-[C]-1", kL6, false},
      {"piperazine", "[NX3]1-[C]-[C]-[NX3]-[C]-[C]-1", kL6, false},
      {"morpholine", "[OX2]1-[C]-[C]-[NX3]-[C]-[C]-1", kL6, false},
      {"thiomorpholine", "[SX2]1-[C]-[C]-[NX3]-[C]-[C]-1", kL6, false},
      {"pyrrolidine", "[NX3]1-[C]-[C]-[C]-[C]-1", kL5, false},
      {"pyrazolidine", "[NX3]1-[NX3]-[C]-[C]-[C]-1", kL5, false},
      {"imidazolidine", "[NX3]1-[C]-[NX3]-[C]-[C]-1", kL5, false},
      {"1,3-oxazolidine", "[OX2]1-[C]-[NX3]-[C]-[C]-1", kL5, false},
      {"1,2-oxazolidine", "[OX2]1-[NX3]-[C]-[C]-[C]-1", kL5, false},
      {"1,3-thiazolidine", "[SX2]1-[C]-[NX3]-[C]-[C]-1", kL5, false},
      {"azetidine", "[NX3]1-[C]-[C]-[C]-1", kL4, false},
      {"aziridine", "[NX3]1-[C]-[C]-1", kL3, false},
      {"oxolane", "[OX2]1-[C]-[C]-[C]-[C]-1", kL5, false},
      {"oxane", "[OX2]1-[C]-[C]-[C]-[C]-[C]-1", kL6, false},
      {"oxetane", "[OX2]1-[C]-[C]-[C]-1", kL4, false},
      {"oxirane", "[OX2]1-[C]-[C]-1", kL3, false},
      {"oxepane", "[OX2]1-[C]-[C]-[C]-[C]-[C]-[C]-1", kL7, false},
      {"azepane", "[NX3]1-[C]-[C]-[C]-[C]-[C]-[C]-1", kL7, false},
      {"azocane", "[NX3]1-[C]-[C]-[C]-[C]-[C]-[C]-[C]-1", kL8, false},
      {"1,4-diazepane", "[NX3]1-[C]-[C]-[NX3]-[C]-[C]-[C]-1", kL7, false},
      {"1,3-diazinane", "[NX3]1-[C]-[NX3]-[C]-[C]-[C]-1", kL6, false},
      {"thiolane", "[SX2]1-[C]-[C]-[C]-[C]-1", kL5, false},
      {"thiane", "[SX2]1-[C]-[C]-[C]-[C]-[C]-1", kL6, false},
      {"thietane", "[SX2]1-[C]-[C]-[C]-1", kL4, false},
      {"thiirane", "[SX2]1-[C]-[C]-1", kL3, false},
      {"1,3-dioxolane", "[OX2]1-[C]-[OX2]-[C]-[C]-1", kL5, false},
      {"1,3-dioxane", "[OX2]1-[C]-[OX2]-[C]-[C]-[C]-1", kL6, false},
      {"1,4-dioxane", "[OX2]1-[C]-[C]-[OX2]-[C]-[C]-1", kL6, false},
      {"1,3-dithiolane", "[SX2]1-[C]-[SX2]-[C]-[C]-1", kL5, false},
      {"1,3-oxathiolane", "[OX2]1-[C]-[SX2]-[C]-[C]-1", kL5, false},
      {"adamantane", "[CX4]12-[CX4]-[CX4]3-[CX4]-[CX4](-[CX4]-2)-[CX4]-[CX4](-[CX4]-3)-[CX4]-1", {"1", "2", "3", "4", "5", "9", "6", "7", "10", "8"}, false},
      // mancude heteromonocycles (P-22.2)
      {"pyridine", "n1ccccc1", kL6, true},
      {"pyrimidine", "n1cnccc1", kL6, true},
      {"pyrazine", "n1ccncc1", kL6, true},
      {"pyridazine", "n1ncccc1", kL6, true},
      {"1,3,5-triazine", "n1cncnc1", kL6, true},
      {"1,2,4-triazine", "n1ncncc1", kL6, true},
      {"1,2,3-triazine", "n1nnccc1", kL6, true},
      {"pyran", "o1ccccc1", kL6, true},
      {"thiopyran", "s1ccccc1", kL6, true},
      {"1,4-dioxine", "o1ccocc1", kL6, true},
      {"1,4-oxazine", "o1ccncc1", kL6, true},
      {"1,3-oxazine", "o1cnccc1", kL6, true},
      {"1,4-thiazine", "s1ccncc1", kL6, true},
      {"1,3-thiazine", "s1cnccc1", kL6, true},
      {"pyrrole", "n1cccc1", kL5, true},
      {"furan", "o1cccc1", kL5, true},
      {"thiophene", "s1cccc1", kL5, true},
      {"imidazole", "n1cncc1", kL5, true},
      {"pyrazole", "n1nccc1", kL5, true},
      {"1,3-oxazole", "o1cncc1", kL5, true},
      {"1,2-oxazole", "o1nccc1", kL5, true},
      {"1,3-thiazole", "s1cncc1", kL5, true},
      {"1,2-thiazole", "s1nccc1", kL5, true},
      {"1,2,4-triazole", "n1ncnc1", kL5, true},
      {"1,2,3-triazole", "n1nncc1", kL5, true},
      {"tetrazole", "n1nnnc1", kL5, true},
      {"1,3,4-oxadiazole", "o1cnnc1", kL5, true},
      {"1,2,4-oxadiazole", "o1ncnc1", kL5, true},
      {"1,2,5-oxadiazole", "o1nccn1", kL5, true},
      {"1,2,3-oxadiazole", "o1nncc1", kL5, true},
      {"1,3,4-thiadiazole", "s1cnnc1", kL5, true},
      {"1,2,4-thiadiazole", "s1ncnc1", kL5, true},
      {"1,2,5-thiadiazole", "s1nccn1", kL5, true},
      {"1,2,3-thiadiazole", "s1nncc1", kL5, true},
      {"1,3-dioxole", "o1cocc1", kL5, true},
      {"1,3-oxathiole", "o1cscc1", kL5, true},
      {"1,3-dithiole", "s1cscc1", kL5, true},
      {"azepine", "n1cccccc1", kL7, true},
      {"oxepine", "o1cccccc1", kL7, true},
      {"thiepine", "s1cccccc1", kL7, true},
      {"1,4-diazepine", "n1ccnccc1", kL7, true},
      {"1,3-diazepine", "n1cncccc1", kL7, true},
      {"azocine", "n1ccccccc1", kL8, true},
      // fused carbocycles (P-25.1.1)
      {"naphthalene", "c1cccc2ccccc12", kNaph, true},
      {"azulene", "c1ccc2cccccc12", kAzul, true},
      {"indene", "c1ccc2ccccc12", kInd, true},
      {"anthracene", "c1cccc2cc3ccccc3cc12", kAnth, true},
      {"phenanthrene", "c1cccc2c3ccccc3ccc12", kPhen, true},
      {"fluorene", "c1cccc2c3ccccc3cc12", kFluor, true},
      {"cyclopenta[a]phenanthrene", "C1CCCC2CCC3C(C21)CCC1C3CCC1",
       {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"}, true},
      // fused heterocycles (P-25.2)
      {"quinoline", "n1cccc2ccccc12", kNaph, true},
      {"isoquinoline", "c1nccc2ccccc12", kNaph, true},
      {"quinazoline", "n1cncc2ccccc12", kNaph, true},
      {"quinoxaline", "n1ccnc2ccccc12", kNaph, true},
      {"cinnoline", "n1nccc2ccccc12", kNaph, true},
      {"phthalazine", "c1nncc2ccccc12", kNaph, true},
      {"pteridine", "n1cncc2nccnc12", kNaph, true},
      {"1,8-naphthyridine", "n1cccc2cccnc12", kNaph, true},
      {"1-benzopyran", "o1cccc2ccccc12", kNaph, true},
      {"2-benzopyran", "c1occc2ccccc12", kNaph, true},
      {"1,4-benzodioxine", "o1ccoc2ccccc12", kNaph, true},
      {"1,4-benzoxazine", "o1ccnc2ccccc12", kNaph, true},
      {"indole", "n1ccc2ccccc12", kInd, true},
      {"isoindole", "c1ncc2ccccc12", kInd, true},
      {"indazole", "n1ncc2ccccc12", kInd, true},
      {"benzimidazole", "n1cnc2ccccc12", kInd, true},
      {"1-benzofuran", "o1ccc2ccccc12", kInd, true},
      {"2-benzofuran", "c1occ2ccccc12", kInd, true},
      {"1-benzothiophene", "s1ccc2ccccc12", kInd, true},
      {"2-benzothiophene", "c1scc2ccccc12", kInd, true},
      {"1,3-benzoxazole", "o1cnc2ccccc12", kInd, true},
      {"1,2-benzoxazole", "o1ncc2ccccc12", kInd, true},
      {"1,3-benzothiazole", "s1cnc2ccccc12", kInd, true},
      {"1,2-benzothiazole", "s1ncc2ccccc12", kInd, true},
      {"1,3-benzodioxole", "o1coc2ccccc12", kInd, true},
      {"indolizine", "c1ccn2ccccc12", {"1", "2", "3", "4", "5", "6", "7", "8", "8a"}, true},
      {"purine", "n1cnc2c(c1)ncn2", {"1", "2", "3", "4", "5", "6", "7", "8", "9"}, true},
      {"carbazole", "c1cccc2c3ccccc3nc12", kFluor, true},
      {"dibenzo[b,d]furan", "c1cccc2oc3ccccc3c12", kDbf, true},
      {"dibenzo[b,d]thiophene", "c1cccc2sc3ccccc3c12", kDbf, true},
      {"xanthene", "c1cccc2oc3ccccc3cc12", kAnth, true},
      {"thioxanthene", "c1cccc2sc3ccccc3cc12", kAnth, true},
      {"acridine", "c1cccc2nc3ccccc3cc12", kAnth, true},
      {"phenothiazine", "c1cccc2sc3ccccc3nc12", kPhz, true},
      {"phenoxazine", "c1cccc2oc3ccccc3nc12", kPhz, true},
      {"phenazine", "c1cccc2nc3ccccc3nc12", kPhz, true},
      {"1,2,3-benzotriazole", "n1nnc2ccccc12", kInd, true},
      {"1,3-benzodioxine", "o1cocc2ccccc12", kNaph, true},
      {"1-benzazepine", "n1ccccc2ccccc12", kBzep, true},
      {"2-benzazepine", "c1ncccc2ccccc12", kBzep, true},
      {"3-benzazepine", "c1cnccc2ccccc12", kBzep, true},
      {"1,4-benzodiazepine", "n1ccncc2ccccc12", kBzep, true},
      {"1,5-benzodiazepine", "n1cccnc2ccccc12", kBzep, true},
      {"dibenzo[b,f]azepine", "c1cccc2nc3ccccc3ccc12", {"1", "2", "3", "4", "4a", "5", "5a", "6", "7", "8", "9", "9a", "10", "11", "11a"}, true},
      {"thieno[2,3-d]pyrimidine", "n1cncc2ccsc12", kPyr, true},
      {"thieno[3,2-d]pyrimidine", "n1cncc2sccc12", kPyr, true},
      {"furo[2,3-d]pyrimidine", "n1cncc2ccoc12", kPyr, true},
      {"pyrrolo[2,3-d]pyrimidine", "n1cncc2ccnc12", kPyr, true},
      {"pyrazolo[3,4-d]pyrimidine", "n1ncc2cncnc12", kInd, true},
      {"pyrazolo[1,5-a]pyrimidine", "n1ccc2ncccn12", {"1", "2", "3", "3a", "4", "5", "6", "7", "8"}, true},
      {"imidazo[1,2-a]pyridine", "n1ccn2ccccc12", {"1", "2", "3", "4", "5", "6", "7", "8", "8a"}, true},
      {"imidazo[1,2-a]pyrimidine", "n1ccn2cccnc12", {"1", "2", "3", "4", "5", "6", "7", "8", "8a"}, true},
      {"[1,2,4]triazolo[1,5-a]pyrimidine", "n1cnc2ncccn12", {"1", "2", "3", "3a", "4", "5", "6", "7", "8"}, true},
      {"[1,3]thiazolo[3,2-a]pyrimidine", "s1ccn2cccnc12", {"1", "2", "3", "4", "5", "6", "7", "8", "8a"}, true},
      {"pyrido[1,2-a]pyrimidine", "n1cccn2ccccc12", {"1", "2", "3", "4", "5", "6", "7", "8", "9", "9a"}, true},
      {"pyrido[2,3-d]pyrimidine", "n1cncc2cccnc12", kNaph, true},
      {"1,6-naphthyridine", "n1cccc2cnccc12", kNaph, true},
      {"1,5-naphthyridine", "n1cccc2ncccc12", kNaph, true},
  };
  return t;
}

// Element-only SMARTS of a template pattern: "n1cc..." -> "[#7]~1~[#6]~[#6]...".
std::string generic_smarts(const char* p) {
  std::string out;
  bool bond_next = false;
  for (const char* c = p; *c; ++c) {
    const char ch = *c;
    int z = 0;
    switch (ch) {
      case 'c': case 'C': z = 6; break;
      case 'n': case 'N': z = 7; break;
      case 'o': case 'O': z = 8; break;
      case 's': case 'S': z = 16; break;
      default: break;
    }
    if (z) {
      if (bond_next) out += '~';
      out += "[#" + std::to_string(z) + "]";
      bond_next = true;
    } else if (ch >= '0' && ch <= '9') {
      out += '~';
      out += ch;
    } else if (ch == '(') {
      out += '(';
      bond_next = true;
    } else if (ch == ')') {
      out += ')';
    }
  }
  return out;
}

int pattern_rings(const char* p) {
  int d = 0;
  bool bracket = false;
  for (const char* c = p; *c; ++c) {
    if (*c == '[') bracket = true;
    if (*c == ']') bracket = false;
    if (!bracket) d += (*c >= '0' && *c <= '9');
  }
  return d / 2;
}

const std::vector<SmartsQuery>& template_queries() {
  static const std::vector<SmartsQuery> q = [] {
    std::vector<SmartsQuery> v;
    for (const auto& t : templates()) {
      Diagnostics d;
      v.push_back(parse_smarts(t.mancude ? generic_smarts(t.pattern) : std::string(t.pattern), d));
    }
    return v;
  }();
  return q;
}

// ---- data ------------------------------------------------------------------------------------

// Seniority of classes (P-41), lower = more senior.
// Sulfonyl and Sulfinyl mark sulfones/sulfoxides (never suffixes).
enum class Cls : int {
  Carboxylate = -6, Sulfonate = -5, Olate = -4, Ammonium = -3, RingCation = -2,  // anions, then cations (P-41: ions first)
  Acid = 1, SulfonicAcid = 2, Phosphonic = 3, Ester = 4, Carbonate = 5, Carbamate = 6, PhosMono = 7, Amide = 8, Thioamide = 9,
  Hydrazide = 10, Sulfonamide = 11, Urea = 12, Thiourea = 13, Nitrile = 14, Aldehyde = 15, Ketone = 16, Thione = 17, Alcohol = 18,
  Thiol = 19, Amine = 20, Imine = 21,
  Sulfonyl = 50, Sulfinyl = 51, None = 99
};

bool is_acyl(Cls c) {
  return c == Cls::Acid || c == Cls::Carboxylate || c == Cls::Ester || c == Cls::Amide || c == Cls::Thioamide || c == Cls::Hydrazide ||
         c == Cls::Nitrile ||
         c == Cls::Aldehyde;
}
bool is_sulfur_suffix(Cls c) { return c == Cls::SulfonicAcid || c == Cls::Sulfonate || c == Cls::Sulfonamide; }

struct NAtom {
  int z = 0, h = 0, q = 0;
  bool arom = false, ring = false;
  std::vector<int> nb, bo;
  std::vector<char> ab;  // aromatic bond
};

struct Named {
  std::string text;  // without its own enclosing marks
  bool compound = false;
  std::vector<Key> stereo;
};

enum class Mode { Molecule, Yl, Ylidene, Oyl, Carbonyl };
enum class PKind { Chain, Carbocycle, Benzene, Template, VonBaeyer };

struct RingSys {
  PKind kind = PKind::Carbocycle;
  int templ = -1;
  std::vector<int> atoms;
  std::vector<std::vector<int>> orders;
  int nrings = 1, nhetero = 0, nN = 0, nO = 0, nS = 0, variety = 0, unsat = 0;
  bool mancude = false;
  int ih = 0;                   // mancude parent needs indicated hydrogen
  std::vector<int> X;           // atoms saturated relative to the mancude parent
  std::vector<int> matchable;   // atoms that carry a double bond in the mancude parent
  int vb[3] = {0, 0, 0};        // von Baeyer bridge lengths
};

struct Spec {
  Mode mode = Mode::Molecule;
  int root = -1, from = -1, link = -1;
  std::string word;  // Carbonyl mode: "carbonyl", "sulfonyl" or "sulfinyl"
  std::vector<char> region;
  Cls principal = Cls::None;
  std::vector<char> pok;  // acyl atoms that may count as principal groups
};

struct Cand {
  PKind kind = PKind::Chain;
  int rsys = -1;
  std::vector<int> atoms;
  std::vector<std::vector<int>> orders;
  std::vector<char> in;
  std::vector<int> suffix_at;
  bool inchain = false;
  std::vector<std::pair<int, Named>> pre;       // parent atom (-1: N locant) and prefix
  std::vector<std::pair<int, int>> exo_double;  // parent atom = substituent root
  std::vector<std::pair<int, int>> ester_o;     // ester oxygen, parent atom carrying that ester suffix
  int nmult = 0, ndouble = 0;
};

struct Score {
  std::vector<double> het, ind, suf, unsat, ene, pre, alpha;
};

int compare_score(const Score& a, const Score& b) {
  if (int c = compare_keys(a.het, b.het)) return c;
  if (int c = compare_keys(a.ind, b.ind)) return c;
  if (int c = compare_keys(a.suf, b.suf)) return c;
  if (int c = compare_keys(a.unsat, b.unsat)) return c;
  if (int c = compare_keys(a.ene, b.ene)) return c;
  if (int c = compare_keys(a.pre, b.pre)) return c;
  for (std::size_t i = 0; i < std::min(a.alpha.size(), b.alpha.size()); ++i) {
    if (a.alpha[i] < b.alpha[i]) return -1;
    if (a.alpha[i] > b.alpha[i]) return 1;
  }
  return 0;
}

struct Built {
  std::string name;
  bool compound = false;
  std::vector<Key> stereo;
  std::vector<std::pair<int, std::string>> ester_o;  // ester oxygen, locant of its suffix
};

// ---- the namer -----------------------------------------------------------------------------

class Namer {
 public:
  Namer(const ChemicalSystemState& st, std::uint32_t sp, std::vector<std::uint32_t> only = {}) : st_(st), sp_(sp), only_(std::move(only)) {}
  std::optional<std::string> run();
  const std::string& error() const { return err_; }

 private:
  const ChemicalSystemState& st_;
  std::uint32_t sp_;
  std::vector<NAtom> at_;
  std::vector<std::uint32_t> glob_;  // local -> global atom
  std::map<std::uint32_t, int> local_;
  std::vector<Cls> acyl_;
  std::vector<char> amineN_, amideN_, nitro_, imineC_, imineN_, cationN_, ringCat_, anionO_;
  std::vector<char> hydrN_, hydrN2_, hydroxamO_;  // hydrazide N (on the acyl) and N'; N-hydroxy amide O
  std::vector<std::uint32_t> only_;  // restrict to these global atoms (one component of a salt)
  std::vector<int> ringsys_;
  std::vector<RingSys> rs_;
  std::map<Key, CipLabel> stereo_;
  std::map<std::tuple<int, int, int>, std::optional<Named>> memo_;
  std::vector<std::vector<char>> ringbond_;
  std::string err_;

  int n() const { return static_cast<int>(at_.size()); }
  int deg(int a) const { return static_cast<int>(at_[a].nb.size()); }
  int order(int a, int b) const {
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
      if (at_[a].nb[k] == b) return at_[a].bo[k];
    return 0;
  }
  bool fail(const std::string& why) {
    if (err_.empty()) err_ = why;
    return false;
  }
  std::optional<Named> fail_n(const std::string& why) {
    fail(why);
    return std::nullopt;
  }
  std::string sym(int a) const {
    return std::string(element_symbol(static_cast<std::uint8_t>(at_[a].z))) + std::to_string(glob_[a] - st_.species[sp_].atom_begin + 1);
  }
  int dO(int c) const {  // =O neighbour of a carbon
    for (std::size_t k = 0; k < at_[c].nb.size(); ++k) {
      const int o = at_[c].nb[k];
      if (at_[c].bo[k] == 2 && at_[o].z == 8 && deg(o) == 1) return o;
    }
    return -1;
  }
  int dS(int c) const {  // =S neighbour of a carbon (thiocarbonyl)
    for (std::size_t k = 0; k < at_[c].nb.size(); ++k) {
      const int x = at_[c].nb[k];
      if (at_[c].bo[k] == 2 && at_[x].z == 16 && deg(x) == 1) return x;
    }
    return -1;
  }
  bool is_Ominus(int o) const { return at_[o].z == 8 && at_[o].q == -1 && deg(o) == 1 && at_[o].h == 0; }
  bool is_OH(int o) const { return at_[o].z == 8 && at_[o].h == 1 && deg(o) == 1 && at_[o].q == 0; }
  bool is_SH(int s) const { return at_[s].z == 16 && at_[s].h == 1 && deg(s) == 1; }
  int ester_O(int c) const {
    for (int x : at_[c].nb)
      if (at_[x].z == 8 && deg(x) == 2) return x;
    return -1;
  }
  int amide_N(int c) const {
    int r = -1;
    for (int x : at_[c].nb)
      if (at_[x].z == 7 && (r < 0 || amideN_[x])) r = x;
    return r;
  }
  int imine_N(int c) const {
    for (std::size_t k = 0; k < at_[c].nb.size(); ++k)
      if (at_[c].bo[k] == 2 && imineN_[at_[c].nb[k]]) return at_[c].nb[k];
    return -1;
  }
  std::string format_prefixes(const std::vector<std::pair<std::string, Named>>& entries, bool omit) const;

  bool load();
  bool perceive();
  bool rings();
  std::vector<char> region_from(int root, const std::vector<char>& blocked) const;
  std::optional<Built> build(const Spec& sp);
  bool analyze(const Spec& sp, Cand& c);
  Score score(const Spec& sp, const Cand& c, const std::vector<int>& o, bool& valid) const;
  std::vector<std::string> locants(const Cand& c, const std::vector<int>& o) const;
  bool perfect_matching(std::vector<int> atoms) const;
  bool von_baeyer(RingSys& R);
  void hydro(const Cand& c, const std::vector<std::string>& L, std::string& ind, std::vector<std::string>& hyd) const;
  std::optional<Built> assemble(const Spec& sp, const Cand& c, const std::vector<int>& o);
  std::string hydride(const Cand& c, const std::vector<int>& o, bool omit_unsat, bool& saturated) const;
  std::optional<Named> sub(int q, int from, int bo);
  std::optional<Named> sub_impl(int q, int from, int bo);
  std::optional<Named> build_sub(Mode mode, int root, int from, int link = -1, const std::string& word = "carbonyl");
  std::optional<Named> sulfonyl_name(int S, int from);
  std::optional<Named> acyl_name(int c, int from);
  std::optional<Named> ligands(const std::vector<int>& xs, int n_atom, const std::string& tail);
};

bool Namer::load() {
  if (sp_ >= st_.species.size()) return fail("no such species");
  const Species& s = st_.species[sp_];
  std::vector<int> g2l(s.atom_count, -1);
  std::vector<char> inside(s.atom_count, only_.empty() ? 1 : 0);
  for (std::uint32_t g : only_) inside[g - s.atom_begin] = 1;
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    const std::uint32_t g = s.atom_begin + i;
    if (!inside[i]) continue;
    const int z = st_.atoms.element[g];
    if (z == 1) {
      const auto nbs = st_.adjacency.neighbors(g);
      if (nbs.size() == 1 && st_.atoms.isotope[g] == 0 && st_.atoms.formal_charge[g] == 0 && st_.atoms.element[nbs[0]] != 1) continue;
      return fail("hydrogen that is not a plain substituent (H2, hydride or isotope)");
    }
    if (st_.atoms.isotope[g] != 0) return fail("isotopically modified atom");
    g2l[i] = n();
    NAtom a;
    a.z = z;
    a.h = st_.atoms.implicit_h[g];
    a.q = st_.atoms.formal_charge[g];
    a.arom = st_.atoms.flags[g] & kAtomAromatic;
    a.ring = st_.atoms.flags[g] & kAtomInRing;
    at_.push_back(a);
    glob_.push_back(g);
    local_[g] = g2l[i];
  }
  for (int a = 0; a < n(); ++a) {
    const std::uint32_t g = glob_[a];
    const auto nbs = st_.adjacency.neighbors(g);
    const auto bds = st_.adjacency.bonds_of(g);
    for (std::size_t k = 0; k < nbs.size(); ++k) {
      const int l = g2l[nbs[k] - s.atom_begin];
      if (l < 0) {
        ++at_[a].h;
        continue;
      }
      at_[a].nb.push_back(l);
      at_[a].bo.push_back(st_.bonds.order[bds[k]]);
      at_[a].ab.push_back((st_.bonds.flags[bds[k]] & kBondAromatic) ? 1 : 0);
    }
  }
  if (n() == 0) return fail("no heavy atoms");
  ringbond_.assign(n(), {});
  for (int a = 0; a < n(); ++a) {
    ringbond_[a].assign(at_[a].nb.size(), 0);
    const auto nbs = st_.adjacency.neighbors(glob_[a]);
    const auto bds = st_.adjacency.bonds_of(glob_[a]);
    for (std::size_t k = 0; k < nbs.size(); ++k) {
      if (!(st_.bonds.flags[bds[k]] & kBondInRing)) continue;
      auto it = local_.find(nbs[k]);
      if (it == local_.end()) continue;
      for (std::size_t j = 0; j < at_[a].nb.size(); ++j)
        if (at_[a].nb[j] == it->second) ringbond_[a][j] = 1;
    }
  }
  std::vector<char> seen(n(), 0);
  std::vector<int> stack = {0};
  seen[0] = 1;
  int count = 1;
  while (!stack.empty()) {
    const int a = stack.back();
    stack.pop_back();
    for (int b : at_[a].nb)
      if (!seen[b]) seen[b] = 1, ++count, stack.push_back(b);
  }
  if (count != n()) return fail("several disconnected components");
  auto local = [&](std::uint32_t g) {
    auto it = local_.find(g);
    return it == local_.end() ? -1 : it->second;
  };
  auto mine = [&](std::uint32_t g) { return g >= s.atom_begin && g < s.atom_begin + s.atom_count && inside[g - s.atom_begin]; };
  for (const auto& c : cip_centres(st_, sp_)) {
    if (c.label == CipLabel::None || !mine(c.index)) continue;
    const int a = local(c.index);
    if (a < 0) return fail("stereocentre on hydrogen");
    stereo_[atom_key(a)] = c.label;
  }
  for (const auto& c : cip_double_bonds(st_, sp_)) {
    if (c.label == CipLabel::None || !mine(st_.bonds.begin[c.index])) continue;
    const int a = local(st_.bonds.begin[c.index]), b = local(st_.bonds.end[c.index]);
    if (a < 0 || b < 0) return fail("stereo double bond to hydrogen");
    stereo_[bond_key(a, b)] = c.label;
  }
  return true;
}

bool Namer::perceive() {
  acyl_.assign(n(), Cls::None);
  amineN_.assign(n(), 0);
  amideN_.assign(n(), 0);
  nitro_.assign(n(), 0);
  imineC_.assign(n(), 0);
  imineN_.assign(n(), 0);
  std::vector<char> nitroO(n(), 0);
  for (int a = 0; a < n(); ++a) {  // nitro groups, charge-separated or pentavalent
    const NAtom& A = at_[a];
    if (A.z != 7) continue;
    int nc = 0, nO = 0, ndO = 0, nmO = 0;
    for (std::size_t k = 0; k < A.nb.size(); ++k) {
      const int b = A.nb[k];
      if (at_[b].z == 6 && A.bo[k] == 1) ++nc;
      else if (at_[b].z == 8 && deg(b) == 1) {
        ++nO;
        if (A.bo[k] == 2 && at_[b].q == 0) ++ndO;
        if (A.bo[k] == 1 && at_[b].q == -1) ++nmO;
      }
    }
    const bool ionic = A.q == 1 && nc == 1 && nO == 2 && ndO == 1 && nmO == 1 && A.h == 0;
    const bool neutral = A.q == 0 && nc == 1 && nO == 2 && ndO == 2 && A.h == 0;
    if (ionic || neutral) {
      nitro_[a] = 1;
      for (int b : A.nb)
        if (at_[b].z == 8) nitroO[b] = 1;
    }
  }
  // ions: ammonium and aromatic N+ cations, O- anions (carboxylate, sulfonate, olate)
  cationN_.assign(n(), 0);
  ringCat_.assign(n(), 0);
  anionO_.assign(n(), 0);
  for (int a = 0; a < n(); ++a) {
    const NAtom& A = at_[a];
    int val = A.h;
    for (int o : A.bo) val += o;
    if (A.z == 7 && A.q == 1 && !nitro_[a] && val == 4) {
      bool ok = true;
      for (int b : A.nb) ok &= at_[b].z == 6 || (A.arom && at_[b].z == 7 && at_[b].arom);
      if (ok && !A.ring && std::all_of(A.bo.begin(), A.bo.end(), [](int o) { return o == 1; })) cationN_[a] = 1;
      else if (ok && A.arom) ringCat_[a] = 1;
    }
    if (is_Ominus(a) && A.bo[0] == 1 && (at_[A.nb[0]].z == 6 || at_[A.nb[0]].z == 16)) anionO_[a] = 1;
  }
  // hydrazides C(=O)-NR-NR2 / C(=O)-NR-N=CR2 and N-hydroxy amides C(=O)-NR-OH
  hydrN_.assign(n(), 0);
  hydrN2_.assign(n(), 0);
  hydroxamO_.assign(n(), 0);
  for (int N = 0; N < n(); ++N) {
    const NAtom& A = at_[N];
    if (A.z != 7 || A.ring || A.arom || A.q != 0) continue;
    int acylC = -1, other = -1, nother = 0;
    bool singles = std::all_of(A.bo.begin(), A.bo.end(), [](int o) { return o == 1; });
    for (int x : A.nb) {
      if (at_[x].z == 6 && !at_[x].ring && dO(x) >= 0 && acylC < 0) acylC = x;
      else if (at_[x].z != 6) other = x, ++nother;
    }
    if (!singles || acylC < 0 || nother != 1) continue;
    int cheteros = 0;
    for (int y : at_[acylC].nb) cheteros += at_[y].z != 6 && y != dO(acylC);
    if (cheteros != 1) continue;  // only plain acyl groups
    if (is_OH(other)) {
      hydroxamO_[other] = 1;
      amideN_[N] = 1;
      continue;
    }
    if (at_[other].z != 7 || at_[other].ring || at_[other].arom || at_[other].q != 0) continue;
    bool ok = true;
    for (std::size_t k = 0; k < at_[other].nb.size(); ++k) {
      const int y = at_[other].nb[k];
      if (y == N) continue;
      if (at_[y].z != 6 || dO(y) >= 0) ok = false;
      if (at_[other].bo[k] == 2)
        for (int w : at_[y].nb) ok &= w == other || at_[w].z == 6;
    }
    if (!ok) continue;
    hydrN_[N] = 1;
    hydrN2_[other] = 1;
  }
  // sulfur and phosphorus oxo groups
  std::vector<char> okO(n(), 0);
  for (int a = 0; a < n(); ++a) {
    const NAtom& A = at_[a];
    if ((A.z != 16 && A.z != 15) || A.arom || A.q != 0) continue;
    if (A.z == 16 && deg(a) == 1 && A.bo[0] == 2 && at_[A.nb[0]].z == 6) continue;  // thiocarbonyl, classified on its carbon
    std::vector<int> dO_, single;
    for (std::size_t k = 0; k < A.nb.size(); ++k) {
      const int b = A.nb[k];
      if (A.bo[k] == 2 && at_[b].z == 8 && deg(b) == 1) dO_.push_back(b);
      else if (A.bo[k] == 1) single.push_back(b);
      else return fail("unsupported multiple bond at " + sym(a));
    }
    if (A.z == 16) {
      if (dO_.empty()) continue;
      if (A.ring || A.h) return fail("cyclic or hydrido sulfur oxide at " + sym(a));
      int nC = 0, nOH = 0, nOm = 0;
      int N = -1;
      for (int b : single) {
        if (at_[b].z == 6) ++nC;
        else if (is_OH(b)) ++nOH;
        else if (is_Ominus(b)) ++nOm;
        else if (at_[b].z == 7 && !nitro_[b]) N = b;
        else return fail("sulfur ester or halide at " + sym(a));
      }
      if (dO_.size() == 1 && nC == 2 && single.size() == 2) acyl_[a] = Cls::Sulfinyl;
      else if (dO_.size() == 2 && single.size() == 2 && nC == 2) acyl_[a] = Cls::Sulfonyl;
      else if (dO_.size() == 2 && single.size() == 2 && nC == 1 && nOH == 1) acyl_[a] = Cls::SulfonicAcid;
      else if (dO_.size() == 2 && single.size() == 2 && nC == 1 && nOm == 1) acyl_[a] = Cls::Sulfonate;
      else if (dO_.size() == 2 && single.size() == 2 && nC == 1 && N >= 0) {
        if (at_[N].ring) acyl_[a] = Cls::Sulfonyl;
        else {
          for (int y : at_[N].nb) {
            if (y == a) continue;
            if (at_[y].z != 6 || dO(y) >= 0) return fail("N-acyl or N-heteroatom sulfonamide at " + sym(N));
          }
          acyl_[a] = Cls::Sulfonamide;
          amideN_[N] = 1;
        }
      } else return fail("unsupported sulfur oxide at " + sym(a));
      for (int b : dO_) okO[b] = 1;
      for (int b : single)
        if (at_[b].z == 8) okO[b] = 1;
    } else {
      int nOH = 0, nC = 0, bridge = -1;
      for (int b : single) {
        if (is_OH(b)) ++nOH;
        else if (at_[b].z == 6) ++nC;
        else if (at_[b].z == 8 && deg(b) == 2) {
          const int y = at_[b].nb[0] == a ? at_[b].nb[1] : at_[b].nb[0];
          if (at_[y].z != 6 || at_[b].ring) return fail("phosphoric anhydride or cyclic phosphate at " + sym(a));
          bridge = b;
        } else return fail("unsupported phosphorus group at " + sym(a));
      }
      if (dO_.size() != 1 || nOH != 2 || single.size() != 3) return fail("phosphorus group other than a phosphate monoester or phosphonic acid at " + sym(a));
      if (nC == 1) acyl_[a] = Cls::Phosphonic;
      else if (bridge >= 0) acyl_[a] = Cls::PhosMono;
      else return fail("unsupported phosphorus group at " + sym(a));
      okO[dO_[0]] = 1;
      for (int b : single)
        if (at_[b].z == 8) okO[b] = 1;
    }
  }
  for (int a = 0; a < n(); ++a) {
    const NAtom& A = at_[a];
    int val = A.h;
    for (int o : A.bo) val += o;
    if (A.q != 0 && !nitro_[a] && !nitroO[a] && !cationN_[a] && !ringCat_[a] && !anionO_[a])
      return fail("charged atom " + sym(a) + " (ion type not covered)");
    switch (A.z) {
      case 6: {
        if (val != 4) return fail("carbon " + sym(a) + " with unusual valence (radical or carbene)");
        int nd = 0, nt = 0, imN = -1, nS = 0;
        for (std::size_t k = 0; k < A.nb.size(); ++k) {
          const int b = A.nb[k];
          if (at_[b].z == 6 || A.ab[k] || (ringbond_[a][k] && at_[b].z != 8)) continue;
          if (A.bo[k] == 2) {
            if (at_[b].z == 8 && deg(b) == 1) ++nd;
            else if (at_[b].z == 16 && deg(b) == 1) ++nS;
            else if (at_[b].z == 7 && imN < 0 && at_[b].q == 0) imN = b;
            else return fail("C=" + std::string(element_symbol(static_cast<std::uint8_t>(at_[b].z))) + " double bond at " + sym(a) + " (imines, thiones, ...)");
          } else if (A.bo[k] == 3) {
            if (at_[b].z == 7 && deg(b) == 1) ++nt;
            else return fail("unsupported triple bond at " + sym(a));
          }
        }
        if (nS) {  // thiocarbonyls: thiones, thioamides, thioureas
          if (nd || nt || imN >= 0 || nS > 1) return fail("unsupported thiocarbonyl at " + sym(a));
          std::vector<int> hs;
          for (std::size_t k = 0; k < A.nb.size(); ++k)
            if (at_[A.nb[k]].z != 6 && A.bo[k] == 1) hs.push_back(A.nb[k]);
          auto plainN = [&](int x) {
            if (at_[x].z != 7 || at_[x].ring || nitro_[x]) return false;
            for (int y : at_[x].nb)
              if (y != a && (at_[y].z != 6 || dO(y) >= 0 || dS(y) >= 0)) return false;
            return true;
          };
          if (A.ring || hs.empty()) acyl_[a] = Cls::Thione;
          else if (hs.size() == 1 && at_[hs[0]].z == 7 && at_[hs[0]].ring) acyl_[a] = Cls::Thione;  // N-thioacyl heterocycle
          else if (hs.size() == 1 && plainN(hs[0])) acyl_[a] = Cls::Thioamide, amideN_[hs[0]] = 1;
          else if (hs.size() == 2 && plainN(hs[0]) && plainN(hs[1])) acyl_[a] = Cls::Thiourea, amideN_[hs[0]] = amideN_[hs[1]] = 1;
          else if (hs.size() == 2 && ((plainN(hs[0]) && at_[hs[1]].z == 7 && at_[hs[1]].ring) || (plainN(hs[1]) && at_[hs[0]].z == 7 && at_[hs[0]].ring)))
            acyl_[a] = Cls::Thioamide, amideN_[at_[hs[0]].ring ? hs[1] : hs[0]] = 1;
          else return fail("thiocarbonyl ester, acid or halide at " + sym(a));
          if (acyl_[a] == Cls::Thioamide && A.h > 1) return fail("unsupported thioformamide at " + sym(a));
          break;
        }
        if (imN >= 0 && hydrN2_[imN]) {  // hydrazone carbon of a hydrazide: an ylidene on N'
          if (nd || nt) return fail("unsupported carbon " + sym(a));
          break;
        }
        if (imN >= 0) {  // imines and oximes (C=N-H, C=N-R, C=N-OH, C=N-OR)
          if (nd || nt) return fail("unsupported carbon " + sym(a));
          for (int b : A.nb)
            if (b != imN && at_[b].z != 6) return fail("amidine, imidate or related C(=N)X group at " + sym(a));
          for (int y : at_[imN].nb) {
            if (y == a) continue;
            if (at_[y].z == 6 && order(imN, y) == 1 && (dO(y) < 0)) continue;
            if (at_[y].z == 8 && order(imN, y) == 1 && (is_OH(y) || (deg(y) == 2 && at_[at_[y].nb[0] == imN ? at_[y].nb[1] : at_[y].nb[0]].z == 6))) continue;
            return fail("hydrazone, azine or N-acyl imine at " + sym(imN));
          }
          if (deg(imN) > 2) return fail("unsupported imine nitrogen at " + sym(imN));
          imineC_[a] = 1;
          imineN_[imN] = 1;
          break;
        }
        if (nt) {
          if (nd || deg(a) > 2) return fail("unsupported nitrile at " + sym(a));
          for (int b : A.nb)
            if (at_[b].z != 6 && !(at_[b].z == 7 && order(a, b) == 3)) return fail("cyanate or cyanamide at " + sym(a));
          acyl_[a] = Cls::Nitrile;
        } else if (nd == 1) {
          if (A.ring) {
            acyl_[a] = Cls::Ketone;  // ketones and heterocyclic pseudoketones (lactams, lactones, cyclic imides)
            break;
          }
          std::vector<int> hs;
          int ncarb = 0;
          for (std::size_t k = 0; k < A.nb.size(); ++k) {
            const int b = A.nb[k];
            if (at_[b].z == 6) ++ncarb;
            else if (A.bo[k] == 1) hs.push_back(b);
          }
          if (hs.empty()) {
            acyl_[a] = ncarb == 2 ? Cls::Ketone : Cls::Aldehyde;
          } else if (hs.size() == 1) {
            const int x = hs[0];
            if (at_[x].z == 8) {
              if (is_OH(x)) acyl_[a] = Cls::Acid;
              else if (is_Ominus(x)) acyl_[a] = Cls::Carboxylate;
              else if (deg(x) == 2) {
                const int y = at_[x].nb[0] == a ? at_[x].nb[1] : at_[x].nb[0];
                if (at_[y].z != 6 || dO(y) >= 0) return fail("anhydride at " + sym(a));
                acyl_[a] = Cls::Ester;
              } else return fail("unsupported acyl group at " + sym(a));
            } else if (at_[x].z == 7) {
              if (at_[x].ring) {
                acyl_[a] = Cls::Ketone;  // N-acyl heterocycle: named as a ketone
              } else if (nitro_[x]) {
                return fail("unsupported acyl group at " + sym(a));
              } else {
                for (int y : at_[x].nb) {
                  if (y == a || hydroxamO_[y] || hydrN2_[y]) continue;
                  if (at_[y].z != 6) return fail("N-heteroatom amide at " + sym(x));
                  if (dO(y) >= 0) return fail("imide or diacylamine at " + sym(x));
                }
                acyl_[a] = hydrN_[x] ? Cls::Hydrazide : Cls::Amide;
                amideN_[x] = 1;
              }
            } else return fail("acyl halide or thioester at " + sym(a));
          } else if (hs.size() == 2) {
            // carbonic acid derivatives: carbonates, carbamates, ureas and ring-N carboxylic acids/esters/amides
            auto esterO = [&](int x) {
              if (at_[x].z != 8 || deg(x) != 2) return false;
              const int y = at_[x].nb[0] == a ? at_[x].nb[1] : at_[x].nb[0];
              return at_[y].z == 6 && dO(y) < 0;
            };
            auto plainN = [&](int x) {
              if (at_[x].z != 7 || at_[x].ring || nitro_[x]) return false;
              for (int y : at_[x].nb)
                if (y != a && (at_[y].z != 6 || dO(y) >= 0)) return false;
              return true;
            };
            const int x = hs[0], y = hs[1];
            const bool ringN_x = at_[x].z == 7 && at_[x].ring, ringN_y = at_[y].z == 7 && at_[y].ring;
            if (ringN_x || ringN_y) {
              const int o = ringN_x ? y : x;
              if (ringN_x && ringN_y) return fail("carbonyl between two ring nitrogens at " + sym(a));
              if (is_OH(o)) acyl_[a] = Cls::Acid;
              else if (esterO(o)) acyl_[a] = Cls::Ester;
              else if (plainN(o)) acyl_[a] = Cls::Amide, amideN_[o] = 1;
              else return fail("unsupported carbonyl group at " + sym(a));
            } else if (esterO(x) && esterO(y)) {
              acyl_[a] = Cls::Carbonate;
            } else if ((esterO(x) && plainN(y)) || (esterO(y) && plainN(x))) {
              acyl_[a] = Cls::Carbamate;
              amideN_[at_[x].z == 7 ? x : y] = 1;
            } else if (plainN(x) && plainN(y)) {
              acyl_[a] = Cls::Urea;
              amideN_[x] = amideN_[y] = 1;
            } else return fail("carbonic acid derivative (carbonic or carbamic acid, halide, ...) at " + sym(a));
          } else return fail("unsupported carbonyl group at " + sym(a));
        } else if (nd > 1) return fail("unsupported carbon " + sym(a));
        break;
      }
      case 7: {
        if (nitro_[a] || A.arom || imineN_[a] || cationN_[a] || ringCat_[a] || hydrN_[a] || hydrN2_[a]) break;
        if (amideN_[a]) {  // N-hydroxy amide nitrogen, checked with its acyl carbon
          bool hx = false;
          for (int b : A.nb) hx |= hydroxamO_[b] != 0;
          if (hx) break;
        }
        if (deg(a) == 1 && A.bo[0] == 3) break;  // nitrile N
        if (val != 3) return fail("nitrogen " + sym(a) + " with unusual valence");
        if (deg(a) == 0) return fail("ammonia");
        for (std::size_t k = 0; k < A.nb.size(); ++k) {
          const int b = A.nb[k];
          if (A.bo[k] != 1 && !A.ab[k] && !ringbond_[a][k]) return fail("N=X double bond at " + sym(a) + " (imines, azo compounds, oximes, ...)");
          const bool sulfonyl_N = at_[b].z == 16 && (acyl_[b] == Cls::Sulfonamide || acyl_[b] == Cls::Sulfonyl);
          if (at_[b].z != 6 && !ringbond_[a][k] && !sulfonyl_N) return fail("N-heteroatom bond at " + sym(a) + " (hydrazines, hydroxylamines, ...)");
        }
        if (!A.ring) amineN_[a] = 1;
        break;
      }
      case 8: {
        if (nitroO[a] || okO[a] || anionO_[a] || hydroxamO_[a]) break;
        {
          bool oxime = false;
          for (int b : A.nb) oxime |= imineN_[b] != 0;
          if (oxime) break;
        }
        if (val != 2) return fail("oxygen " + sym(a) + " with unusual valence");
        if (deg(a) == 0) return fail("water");
        for (std::size_t k = 0; k < A.nb.size(); ++k)
          if (at_[A.nb[k]].z != 6 && !ringbond_[a][k])
            return fail("O-" + std::string(element_symbol(static_cast<std::uint8_t>(at_[A.nb[k]].z))) + " bond at " + sym(a) + " (peroxides, N-oxides, sulfur or phosphorus esters, ...)");
        break;
      }
      case 15:
        if (acyl_[a] == Cls::None) return fail("unsupported phosphorus at " + sym(a));
        break;
      case 16: {
        if (A.arom || acyl_[a] != Cls::None) break;
        if (deg(a) == 1 && A.bo[0] == 2 && at_[A.nb[0]].z == 6) break;  // thiocarbonyl
        if (val != 2) return fail("sulfur " + sym(a) + " in a higher oxidation state (sulfoxides, sulfones, sulfonic acids)");
        if (deg(a) == 0) return fail("hydrogen sulfide");
        for (std::size_t k = 0; k < A.nb.size(); ++k)
          if (at_[A.nb[k]].z != 6 || A.bo[k] != 1) return fail("unsupported sulfur bond at " + sym(a));
        break;
      }
      case 9: case 17: case 35: case 53:
        if (deg(a) != 1 || A.bo[0] != 1 || at_[A.nb[0]].z != 6) return fail("halogen " + sym(a) + " not bonded to carbon");
        break;
      default:
        return fail("element " + std::string(element_symbol(static_cast<std::uint8_t>(A.z))) + " is outside the naming scope (C, H, N, O, S, halogens)");
    }
  }
  for (int a = 0; a < n(); ++a)
    if (amideN_[a] || nitro_[a] || imineN_[a] || cationN_[a] || hydrN_[a] || hydrN2_[a] || at_[a].ring || at_[a].arom) amineN_[a] = 0;
  return true;
}

bool Namer::rings() {
  ringsys_.assign(n(), -1);
  std::vector<int> par(n());
  for (int i = 0; i < n(); ++i) par[i] = i;
  auto find = [&](int x) {
    while (par[x] != x) x = par[x] = par[par[x]];
    return x;
  };
  const auto& ringbond = ringbond_;
  for (int a = 0; a < n(); ++a)
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
      if (ringbond[a][k]) par[find(a)] = find(at_[a].nb[k]);
  std::map<int, std::vector<int>> groups;
  for (int a = 0; a < n(); ++a)
    if (at_[a].ring) groups[find(a)].push_back(a);
  for (auto& [rep, atoms] : groups) {
    RingSys R;
    R.atoms = atoms;
    std::set<int> in(atoms.begin(), atoms.end());
    int nbonds = 0;
    for (int a : atoms)
      for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
        if (ringbond[a][k] && at_[a].nb[k] > a) ++nbonds;
    R.nrings = nbonds - static_cast<int>(atoms.size()) + 1;
    bool allC = true, allArom = true, anyArom = false;
    for (int a : atoms) {
      allC &= at_[a].z == 6;
      allArom &= at_[a].arom;
      anyArom |= at_[a].arom;
      if (at_[a].z != 6) {
        ++R.nhetero;
        R.nN += at_[a].z == 7;
        R.nO += at_[a].z == 8;
        R.nS += at_[a].z == 16;
      }
    }
    R.variety = (R.nN > 0) + (R.nO > 0) + (R.nS > 0);
    // monocycles: carbocycles, and heteromonocycles of more than ten members by replacement
    // nomenclature (P-22.2.3, "1-oxacyclohexadecane"); smaller heterocycles use the templates
    if (R.nrings == 1 && (allC || (atoms.size() > 10 && !anyArom))) {
      std::vector<int> cyc = {atoms[0]};
      int prev = -1, cur = atoms[0];
      while (cyc.size() <= atoms.size()) {
        int next = -1;
        for (std::size_t k = 0; k < at_[cur].nb.size(); ++k)
          if (ringbond[cur][k] && at_[cur].nb[k] != prev) {
            next = at_[cur].nb[k];
            break;
          }
        if (next < 0 || next == atoms[0]) break;
        cyc.push_back(next);
        prev = cur;
        cur = next;
      }
      if (cyc.size() != atoms.size()) return fail("ring perception");
      const int m = static_cast<int>(cyc.size());
      for (int st = 0; st < m; ++st)
        for (int dir : {1, -1}) {
          std::vector<int> o;
          for (int i = 0; i < m; ++i) o.push_back(cyc[((st + dir * i) % m + m) % m]);
          R.orders.push_back(o);
        }
      if (allArom && m == 6) {
        R.kind = PKind::Benzene;
        R.unsat = 100;
      } else if (anyArom) {
        return fail("aromatic carbocycle other than benzene");
      } else {
        R.kind = PKind::Carbocycle;
        for (int i = 0; i < m; ++i)
          if (order(cyc[i], cyc[(i + 1) % m]) > 1) ++R.unsat;
      }
    } else {
      const auto& qs = template_queries();
      SmartsMatchOptions mo;
      mo.unique = false;
      mo.max_matches = 0;
      for (std::size_t t = 0; t < qs.size() && R.orders.empty(); ++t) {
        if (!qs[t].valid() || templates()[t].locants.size() != atoms.size() || pattern_rings(templates()[t].pattern) != R.nrings) continue;
        for (const auto& mt : match_smarts(st_, sp_, qs[t], mo)) {
          std::vector<int> o;
          bool ok = true;
          for (std::uint32_t g : mt) {
            auto it = local_.find(g);
            if (it == local_.end() || !in.count(it->second)) {
              ok = false;
              break;
            }
            o.push_back(it->second);
          }
          if (ok) {
            R.templ = static_cast<int>(t);
            R.orders.push_back(o);
          }
        }
      }
      if (R.orders.empty() && R.nrings == 2 && von_baeyer(R)) {
        const int id = static_cast<int>(rs_.size());
        for (int a : atoms) ringsys_[a] = id;
        rs_.push_back(std::move(R));
        continue;
      }
      if (R.orders.empty())
        return fail(std::string(R.nrings == 1 ? "monocycle" : "ring system") + " of " + std::to_string(atoms.size()) +
                    " atoms without a supported parent name (fused, bridged, spiro, partially hydrogenated or uncommon heterocycle)");
      R.kind = PKind::Template;
      R.mancude = templates()[R.templ].mancude;
      for (int a : atoms)
        for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
          if (ringbond[a][k] && at_[a].nb[k] > a && at_[a].bo[k] == 2) ++R.unsat;
      if (R.mancude) {
        for (int a : atoms) {
          int rdeg = 0;
          bool dbl = false;
          for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
            if (ringbond[a][k]) ++rdeg, dbl |= at_[a].bo[k] == 2;
          const bool inherent = ((at_[a].z == 8 || at_[a].z == 16) && rdeg == 2) || (at_[a].z == 7 && rdeg == 3);
          if (inherent) continue;
          R.matchable.push_back(a);
          if (!dbl) R.X.push_back(a);
        }
        R.ih = static_cast<int>(R.matchable.size() % 2);
      }
    }
    const int id = static_cast<int>(rs_.size());
    for (int a : atoms) ringsys_[a] = id;
    rs_.push_back(std::move(R));
  }
  return true;
}

std::vector<char> Namer::region_from(int root, const std::vector<char>& blocked) const {
  std::vector<char> reg(n(), 0);
  std::vector<int> stack = {root};
  reg[root] = 1;
  while (!stack.empty()) {
    const int a = stack.back();
    stack.pop_back();
    for (int b : at_[a].nb)
      if (!reg[b] && !blocked[b]) reg[b] = 1, stack.push_back(b);
  }
  return reg;
}

// ---- substituent prefixes ------------------------------------------------------------------

std::optional<Named> Namer::sub(int q, int from, int bo) {
  const auto key = std::make_tuple(q, from, bo);
  if (auto it = memo_.find(key); it != memo_.end()) return it->second;
  auto r = sub_impl(q, from, bo);
  memo_[key] = r;
  return r;
}

std::optional<Named> Namer::build_sub(Mode mode, int root, int from, int link, const std::string& word) {
  Spec sp;
  sp.mode = mode;
  sp.root = root;
  sp.from = from;
  sp.link = link;
  sp.word = word;
  std::vector<char> blocked(n(), 0);
  blocked[from] = 1;
  sp.region = region_from(root, blocked);
  sp.pok.assign(n(), 0);
  auto b = build(sp);
  if (!b) return std::nullopt;
  Named r{b->name, b->compound || (mode == Mode::Carbonyl && b->name != "benzoyl"), b->stereo};
  if (r.text == "phenylmethyl") r = Named{"benzyl", false, {}};
  return r;
}

// Acyl prefix of carbonyl carbon c (its =O implied), reached from `from`.
std::optional<Named> Namer::acyl_name(int c, int from) {
  const int o = dO(c);
  std::vector<int> others;
  for (int x : at_[c].nb)
    if (x != from && x != o) others.push_back(x);
  if (others.empty()) return Named{"formyl", false, {}};
  if (others.size() != 1) return fail_n("unsupported acyl group at " + sym(c));
  const int x = others[0];
  if (at_[x].ring) return build_sub(Mode::Carbonyl, c, from, x);
  if (at_[x].z == 6) return build_sub(Mode::Oyl, c, from);
  if (at_[x].z == 8 && deg(x) == 2) {  // alkoxycarbonyl
    auto alk = sub(x, c, 1);
    if (!alk) return std::nullopt;
    return Named{(alk->compound ? enclose(alk->text) : alk->text) + "carbonyl", true, alk->stereo};
  }
  if (at_[x].z == 7 && amideN_[x]) {  // carbamoyl
    std::vector<int> xs;
    for (int y : at_[x].nb)
      if (y != c) xs.push_back(y);
    if (xs.empty()) return Named{"carbamoyl", false, {}};
    return ligands(xs, x, "carbamoyl");
  }
  return fail_n("unsupported acyl group at " + sym(c));
}

// "methanesulfonyl", "benzenesulfonyl", "methanesulfinyl": sulfur S reached from `from`.
std::optional<Named> Namer::sulfonyl_name(int S, int from) {
  int link = -1;
  for (int x : at_[S].nb)
    if (x != from && at_[x].z != 8) link = x;
  if (link < 0) return fail_n("unsupported sulfonyl group at " + sym(S));
  return build_sub(Mode::Carbonyl, S, from, link, acyl_[S] == Cls::Sulfinyl ? "sulfinyl" : "sulfonyl");
}

// "ethyl(methyl)" + tail for the substituents xs of nitrogen n_atom.
std::optional<Named> Namer::ligands(const std::vector<int>& xs, int n_atom, const std::string& tail) {
  Named out;
  out.compound = true;
  std::map<std::string, std::pair<int, bool>> groups;
  for (int x : xs) {
    const bool acyl = at_[x].z == 6 && !at_[x].ring && dO(x) >= 0;
    const bool sulf = at_[x].z == 16 && (acyl_[x] == Cls::Sulfonamide || acyl_[x] == Cls::Sulfonyl);
    std::optional<Named> s = acyl ? acyl_name(x, n_atom) : sulf ? sulfonyl_name(x, n_atom) : sub(x, n_atom, order(x, n_atom));
    if (!s) return std::nullopt;
    out.stereo.insert(out.stereo.end(), s->stereo.begin(), s->stereo.end());
    auto& g = groups[s->text];
    ++g.first;
    g.second = s->compound;
  }
  std::vector<std::pair<std::string, std::pair<int, bool>>> v(groups.begin(), groups.end());
  std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
    const std::string ka = alpha_key(a.first), kb = alpha_key(b.first);
    return ka != kb ? ka < kb : a.first < b.first;
  });
  std::string s;
  for (std::size_t i = 0; i < v.size(); ++i) {
    const auto& [t, g] = v[i];
    std::string piece = g.second ? (g.first > 1 ? multiplier(g.first, true) : "") + enclose(t) : multiplier(g.first, false) + t;
    if (i > 0 && !is_enclosed(piece) && !(g.second && g.first > 1)) piece = enclose(piece);
    s += piece;
  }
  out.text = s + tail;
  return out;
}

std::optional<Named> Namer::sub_impl(int q, int from, int bo) {
  const NAtom& A = at_[q];
  if (bo == 3) return fail_n("ylidyne substituent at " + sym(q));
  if (bo == 2) {
    if (A.z == 8 && deg(q) == 1) return Named{"oxo", false, {}};
    if (A.z == 16 && deg(q) == 1) return Named{"sulfanylidene", false, {}};
    if (A.z == 7 && imineN_[q]) {  // =N-H, =N-R, =N-OH, =N-OR
      int x = -1;
      for (int y : A.nb)
        if (y != from) x = y;
      if (x < 0) return Named{"imino", false, {}};
      auto R = sub(x, q, 1);
      if (!R) return std::nullopt;
      return Named{(R->compound ? enclose(R->text) : R->text) + "imino", true, R->stereo};
    }
    if (A.z == 6) return build_sub(Mode::Ylidene, q, from);
    return fail_n("double-bonded heteroatom substituent at " + sym(q));
  }
  if (A.ring) return build_sub(Mode::Yl, q, from);
  std::vector<int> others;
  for (int x : A.nb)
    if (x != from) others.push_back(x);
  switch (A.z) {
    case 9: return Named{"fluoro", false, {}};
    case 17: return Named{"chloro", false, {}};
    case 35: return Named{"bromo", false, {}};
    case 53: return Named{"iodo", false, {}};
    case 8: {
      if (others.empty()) return Named{A.q == -1 ? "oxido" : "hydroxy", false, {}};
      if (others.size() == 1 && at_[others[0]].z == 15 && acyl_[others[0]] == Cls::PhosMono) return Named{"phosphonooxy", false, {}};
      if (others.size() == 1 && imineN_[others[0]]) {
        auto R = sub(others[0], q, 1);
        if (!R) return std::nullopt;
        return Named{enclose(R->text) + "oxy", true, R->stereo};
      }
      if (others.size() != 1 || at_[others[0]].z != 6) return fail_n("unsupported oxygen substituent at " + sym(q));
      const int x = others[0];
      if (dO(x) >= 0 && !at_[x].ring) {
        auto ac = acyl_name(x, q);
        if (!ac) return std::nullopt;
        return Named{(ac->compound ? enclose(ac->text) : ac->text) + "oxy", true, ac->stereo};
      }
      auto R = sub(x, q, 1);
      if (!R) return std::nullopt;
      static const std::map<std::string, std::string> contracted = {
          {"methyl", "methoxy"}, {"ethyl", "ethoxy"}, {"propyl", "propoxy"}, {"butyl", "butoxy"}, {"phenyl", "phenoxy"}};
      for (const auto& [yl, oxy] : contracted)
        if (ends_with(R->text, yl) && (R->text.size() == yl.size() || (!ends_with(R->text, "cyclo" + yl) && R->text[R->text.size() - yl.size() - 1] != '-')))
          return Named{R->text.substr(0, R->text.size() - yl.size()) + oxy, R->compound, R->stereo};
      return Named{(R->compound ? enclose(R->text) : R->text) + "oxy", true, R->stereo};
    }
    case 15:
      if (acyl_[q] == Cls::Phosphonic) return Named{"phosphono", false, {}};
      return fail_n("unsupported phosphorus substituent at " + sym(q));
    case 16: {
      switch (acyl_[q]) {
        case Cls::SulfonicAcid: return Named{"sulfo", false, {}};
        case Cls::Sulfonate: return Named{"sulfonato", false, {}};
        case Cls::Sulfonamide: {
          int N = -1;
          for (int x : A.nb)
            if (at_[x].z == 7) N = x;
          std::vector<int> xs;
          for (int x : at_[N].nb)
            if (x != q) xs.push_back(x);
          if (xs.empty()) return Named{"sulfamoyl", false, {}};
          return ligands(xs, N, "sulfamoyl");
        }
        case Cls::Sulfonyl:
        case Cls::Sulfinyl: return sulfonyl_name(q, from);
        default: break;
      }
      if (others.empty()) return Named{"sulfanyl", false, {}};
      if (others.size() != 1 || at_[others[0]].z != 6) return fail_n("unsupported sulfur substituent at " + sym(q));
      auto R = sub(others[0], q, 1);
      if (!R) return std::nullopt;
      return Named{(R->compound ? enclose(R->text) : R->text) + "sulfanyl", true, R->stereo};
    }
    case 7: {
      if (nitro_[q]) return Named{"nitro", false, {}};
      if (cationN_[q]) {
        if (others.empty()) return Named{"azaniumyl", false, {}};
        return ligands(others, q, "azaniumyl");
      }
      if (imineN_[q]) {  // -N=CR2: (ylideneamino)
        int c = -1;
        for (std::size_t k = 0; k < A.nb.size(); ++k)
          if (A.bo[k] == 2) c = A.nb[k];
        if (c < 0 || c == from) return fail_n("unsupported imine at " + sym(q));
        auto R = sub(c, q, 2);
        if (!R) return std::nullopt;
        return Named{(R->compound ? enclose(R->text) : R->text) + "amino", true, R->stereo};
      }
      if (others.empty()) return Named{"amino", false, {}};
      if (others.size() == 1 && hydrN2_[others[0]] && deg(others[0]) == 1) return Named{"hydrazinyl", false, {}};
      if (others.size() == 1 && at_[others[0]].z == 16 && acyl_[others[0]] == Cls::Sulfonamide) {
        auto sn = sulfonyl_name(others[0], q);
        if (!sn) return std::nullopt;
        if (!ends_with(sn->text, "sulfonyl")) return Named{enclose(sn->text) + "amino", true, sn->stereo};
        return Named{sn->text.substr(0, sn->text.size() - 8) + "sulfonamido", sn->compound, sn->stereo};
      }
      if (others.size() == 1 && at_[others[0]].z == 6 && !at_[others[0]].ring && dO(others[0]) >= 0 &&
          acyl_[others[0]] != Cls::Carbamate && acyl_[others[0]] != Cls::Urea) {
        auto ac = acyl_name(others[0], q);
        if (!ac) return std::nullopt;
        std::string t = ac->text;
        if (ends_with(t, "acetyl")) t = t.substr(0, t.size() - 6) + "acetamido";
        else if (ends_with(t, "formyl")) t = t.substr(0, t.size() - 6) + "formamido";
        else if (ends_with(t, "oyl")) t = t.substr(0, t.size() - 3) + "amido";
        else if (ends_with(t, "carbonyl")) t = t.substr(0, t.size() - 8) + "carboxamido";
        else return Named{enclose(t) + "amino", true, ac->stereo};
        return Named{t, ac->compound, ac->stereo};
      }
      return ligands(others, q, "amino");
    }
    case 6: {
      if (acyl_[q] == Cls::Nitrile) return Named{"cyano", false, {}};
      if (acyl_[q] == Cls::Thioamide || acyl_[q] == Cls::Thiourea) {  // carbamothioyl, from the carbon side only
        int N = -1;
        for (int x : A.nb)
          if (x != from && at_[x].z == 7 && amideN_[x]) N = x;
        if (N < 0 || (acyl_[q] == Cls::Thioamide && at_[from].z == 7)) return fail_n("thioacyl substituent at " + sym(q));
        std::vector<int> xs;
        for (int x : at_[N].nb)
          if (x != q) xs.push_back(x);
        if (xs.empty()) return Named{"carbamothioyl", false, {}};
        return ligands(xs, N, "carbamothioyl");
      }
      if (dO(q) >= 0) {
        switch (acyl_[q]) {
          case Cls::Acid: return Named{"carboxy", false, {}};
          case Cls::Carboxylate: return Named{"carboxylato", false, {}};
          case Cls::Ester: {
            auto alk = sub(ester_O(q), q, 1);
            if (!alk) return std::nullopt;
            return Named{(alk->compound ? enclose(alk->text) : alk->text) + "carbonyl", true, alk->stereo};
          }
          case Cls::Amide: {
            const int N = amide_N(q);
            std::vector<int> xs;
            for (int x : at_[N].nb)
              if (x != q) xs.push_back(x);
            if (xs.empty()) return Named{"carbamoyl", false, {}};
            return ligands(xs, N, "carbamoyl");
          }
          case Cls::Aldehyde: return Named{"formyl", false, {}};
          default: return acyl_name(q, from);
        }
      }
      return build_sub(Mode::Yl, q, from);
    }
    default: return fail_n("unsupported substituent at " + sym(q));
  }
}

// ---- parents -------------------------------------------------------------------------------

bool Namer::analyze(const Spec& sp, Cand& c) {
  const bool mol = sp.mode == Mode::Molecule;
  c.in.assign(n(), 0);
  for (int a : c.atoms) c.in[a] = 1;
  std::vector<char> used(c.in);
  const bool acylP = mol && (is_acyl(sp.principal) || is_sulfur_suffix(sp.principal));
  auto principal_acyl = [&](int a) { return acylP && acyl_[a] == sp.principal && sp.pok[a]; };
  if (c.kind == PKind::Chain)
    for (int a : c.atoms)
      if (principal_acyl(a)) c.inchain = true;
  for (int a : c.atoms)
    if (ringCat_[a] && !(mol && sp.principal == Cls::RingCation)) return fail("ring cation that is not the principal characteristic group");
  bool nprefix = false;
  auto consume_acyl = [&](int A, int at) -> bool {
    for (int b : at_[A].nb) {
      if (c.in[b] || b == sp.from) continue;
      if (at_[b].z == 6) return fail("acyl carbon " + sym(A) + " with a second carbon");
      used[b] = 1;
      if (acyl_[A] == Cls::Ester && at_[b].z == 8 && deg(b) == 2) c.ester_o.push_back({b, at});
      if (acyl_[A] == Cls::Hydrazide && at_[b].z == 7)
        for (int x : at_[b].nb) {
          if (x == A) continue;
          if (hydrN2_[x]) {
            used[x] = 1;
            for (std::size_t k = 0; k < at_[x].nb.size(); ++k) {
              const int y = at_[x].nb[k];
              if (y == b) continue;
              auto s = sub(y, x, at_[x].bo[k]);
              if (!s) return false;
              c.pre.push_back({-2, *s});
              nprefix = true;
            }
            continue;
          }
          auto s = sub(x, b, order(x, b));
          if (!s) return false;
          c.pre.push_back({-1, *s});
          nprefix = true;
        }
      if ((acyl_[A] == Cls::Amide || acyl_[A] == Cls::Thioamide || acyl_[A] == Cls::Sulfonamide) && at_[b].z == 7)
        for (int x : at_[b].nb) {
          if (x == A) continue;
          auto s = sub(x, b, order(x, b));
          if (!s) return false;
          c.pre.push_back({-1, *s});
          nprefix = true;
        }
    }
    return true;
  };
  for (int p : c.atoms) {
    if (c.inchain && principal_acyl(p)) {
      c.suffix_at.push_back(p);
      if (!consume_acyl(p, p)) return false;
      continue;
    }
    if (mol && sp.principal == Cls::Ketone && acyl_[p] == Cls::Ketone) {
      used[dO(p)] = 1;
      c.suffix_at.push_back(p);
    }
    if (mol && sp.principal == Cls::RingCation && ringCat_[p]) c.suffix_at.push_back(p);
    if (mol && sp.principal == Cls::Thione && acyl_[p] == Cls::Thione) {
      used[dS(p)] = 1;
      c.suffix_at.push_back(p);
    }
    if (mol && sp.principal == Cls::Imine && imineC_[p]) {
      const int N = imine_N(p);
      used[N] = 1;
      c.suffix_at.push_back(p);
      c.exo_double.push_back({p, N});
      for (int x : at_[N].nb) {
        if (x == p) continue;
        auto s2 = sub(x, N, 1);
        if (!s2) return false;
        c.pre.push_back({-1, *s2});
        nprefix = true;
      }
    }
    if (sp.mode == Mode::Oyl && p == sp.root) used[dO(p)] = 1;
    for (std::size_t k = 0; k < at_[p].nb.size(); ++k) {
      const int q = at_[p].nb[k];
      const int bo = at_[p].bo[k];
      if (used[q] || q == sp.from) continue;
      if (sp.mode == Mode::Carbonyl && q == sp.root) continue;
      if (mol && bo == 1) {
        if ((sp.principal == Cls::Alcohol && is_OH(q)) || (sp.principal == Cls::Thiol && is_SH(q)) ||
            (sp.principal == Cls::Olate && is_Ominus(q) && at_[p].z == 6 && acyl_[p] != Cls::Carboxylate)) {
          used[q] = 1;
          c.suffix_at.push_back(p);
          continue;
        }
        if ((sp.principal == Cls::Amine && amineN_[q]) || (sp.principal == Cls::Ammonium && cationN_[q])) {
          used[q] = 1;
          c.suffix_at.push_back(p);
          for (int x : at_[q].nb) {
            if (x == p) continue;
            if (c.in[x]) return fail("amine nitrogen bridging the parent");
            auto s = sub(x, q, order(x, q));
            if (!s) return false;
            c.pre.push_back({-1, *s});
            nprefix = true;
          }
          continue;
        }
        if (acylP && !c.inchain && principal_acyl(q)) {
          used[q] = 1;
          c.suffix_at.push_back(p);
          if (!consume_acyl(q, p)) return false;
          continue;
        }
      }
      auto s = sub(q, p, bo);
      if (!s) return false;
      used[q] = 1;
      c.pre.push_back({p, *s});
      if (bo == 2 && at_[q].z != 6) c.exo_double.push_back({p, q});  // C=C is cited by the ylidene prefix itself
    }
  }
  if (nprefix && (sp.principal == Cls::Amine || sp.principal == Cls::Ammonium || sp.principal == Cls::Amide || sp.principal == Cls::Hydrazide || sp.principal == Cls::Thioamide || sp.principal == Cls::Sulfonamide ||
                  sp.principal == Cls::Imine) &&
      c.suffix_at.size() > 1)
    return fail("N-substituents on several amine or amide groups (N,N' locants)");
  if (c.kind == PKind::VonBaeyer) c.nmult = rs_[c.rsys].unsat;
  if (c.kind == PKind::Chain || c.kind == PKind::Carbocycle) {
    const auto& o = c.orders[0];
    const std::size_t m = o.size();
    const std::size_t lim = c.kind == PKind::Chain ? (m ? m - 1 : 0) : m;
    for (std::size_t i = 0; i < lim; ++i) {
      const int b = order(o[i], o[(i + 1) % m]);
      if (b > 1) ++c.nmult;
      if (b == 2) ++c.ndouble;
    }
  }
  return true;
}

std::vector<std::string> Namer::locants(const Cand& c, const std::vector<int>& o) const {
  std::vector<std::string> L(n());
  const bool templ = c.kind == PKind::Template;
  for (std::size_t i = 0; i < o.size(); ++i) L[o[i]] = templ ? templates()[rs_[c.rsys].templ].locants[i] : std::to_string(i + 1);
  return L;
}

// Bicyclic von Baeyer systems (P-23.2): bridged bicycles, and ortho-fused bicycles with a ring
// smaller than five members (larger ortho-fused systems need fusion names, which are not built).
bool Namer::von_baeyer(RingSys& R) {
  std::vector<int> heads;
  auto rdeg = [&](int a) {
    int d = 0;
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k) d += ringbond_[a][k];
    return d;
  };
  for (int a : R.atoms) {
    const int d = rdeg(a);
    if (d == 3) heads.push_back(a);
    else if (d != 2) return fail("spiro or polycyclic ring junction at " + sym(a));
    if (at_[a].arom) return fail("aromatic atom in a von Baeyer ring system at " + sym(a));
  }
  if (heads.size() != 2) return fail("unsupported bicyclic ring system");
  const int h1 = heads[0], h2 = heads[1];
  std::vector<std::vector<int>> bridges;  // atoms from h1 towards h2
  for (std::size_t k = 0; k < at_[h1].nb.size(); ++k) {
    if (!ringbond_[h1][k]) continue;
    std::vector<int> path;
    int prev = h1, cur = at_[h1].nb[k];
    while (cur != h2) {
      path.push_back(cur);
      int next = -1;
      for (std::size_t j = 0; j < at_[cur].nb.size(); ++j)
        if (ringbond_[cur][j] && at_[cur].nb[j] != prev) next = at_[cur].nb[j];
      if (next < 0 || path.size() > R.atoms.size()) return fail("ring perception");
      prev = cur;
      cur = next;
    }
    bridges.push_back(path);
  }
  if (bridges.size() != 3) return fail("unsupported bicyclic ring system");
  std::sort(bridges.begin(), bridges.end(), [](const auto& x, const auto& y) { return x.size() > y.size(); });
  const int a = static_cast<int>(bridges[0].size()), b = static_cast<int>(bridges[1].size()), c = static_cast<int>(bridges[2].size());
  if (c == 0 && b + 2 >= 5) return fail("ortho-fused bicyclic system without a supported fusion name");
  R.kind = PKind::VonBaeyer;
  R.vb[0] = a;
  R.vb[1] = b;
  R.vb[2] = c;
  std::vector<int> perm = {0, 1, 2};
  do {
    if (bridges[perm[0]].size() != static_cast<std::size_t>(a) || bridges[perm[1]].size() != static_cast<std::size_t>(b)) continue;
    for (int flip : {0, 1}) {
      auto A = bridges[perm[0]], B = bridges[perm[1]], C = bridges[perm[2]];
      int s1 = h1, s2 = h2;
      if (flip) {
        std::reverse(A.begin(), A.end());
        std::reverse(B.begin(), B.end());
        std::reverse(C.begin(), C.end());
        std::swap(s1, s2);
      }
      std::vector<int> o = {s1};
      o.insert(o.end(), A.begin(), A.end());
      o.push_back(s2);
      o.insert(o.end(), B.rbegin(), B.rend());
      o.insert(o.end(), C.begin(), C.end());
      R.orders.push_back(o);
    }
  } while (std::next_permutation(perm.begin(), perm.end()));
  for (int x : R.atoms)
    for (std::size_t k = 0; k < at_[x].nb.size(); ++k)
      if (ringbond_[x][k] && at_[x].nb[k] > x && at_[x].bo[k] > 1) ++R.unsat;
  return true;
}

// Whether the ring-bond graph on `atoms` has a perfect matching (backtracking; ring systems are small).
bool Namer::perfect_matching(std::vector<int> atoms) const {
  if (atoms.empty()) return true;
  const int a = atoms[0];
  for (std::size_t k = 0; k < at_[a].nb.size(); ++k) {
    if (!ringbond_[a][k]) continue;
    auto it = std::find(atoms.begin() + 1, atoms.end(), at_[a].nb[k]);
    if (it == atoms.end()) continue;
    std::vector<int> rest;
    for (auto j = atoms.begin() + 1; j != atoms.end(); ++j)
      if (j != it) rest.push_back(*j);
    if (perfect_matching(rest)) return true;
  }
  return false;
}

// Indicated hydrogen (lowest feasible locant) and hydro prefixes for a mancude template under numbering L.
void Namer::hydro(const Cand& c, const std::vector<std::string>& L, std::string& ind, std::vector<std::string>& hyd) const {
  ind.clear();
  hyd.clear();
  if (c.kind != PKind::Template) return;
  const RingSys& R = rs_[c.rsys];
  if (!R.mancude) return;
  int pick = -1;
  if (R.ih) {
    std::vector<int> xs = R.X;
    std::sort(xs.begin(), xs.end(), [&](int a, int b) { return locant_key(L[a]) < locant_key(L[b]); });
    for (int x : xs) {
      std::vector<int> rest;
      for (int a : R.matchable)
        if (a != x) rest.push_back(a);
      if (perfect_matching(rest)) {
        pick = x;
        break;
      }
    }
    if (pick < 0 && !xs.empty()) pick = xs[0];
    if (pick >= 0) ind = L[pick];
  }
  for (int x : R.X)
    if (x != pick) hyd.push_back(L[x]);
}

// Prefixes with locants (P-14.5 alphanumerical order, P-16.3 multiplying prefixes, P-16.5 enclosing marks).
std::string Namer::format_prefixes(const std::vector<std::pair<std::string, Named>>& entries, bool omit) const {
  struct Group {
    std::vector<std::string> locs;
    bool compound = false;
    int count = 0;
  };
  std::map<std::string, Group> groups;
  for (const auto& [l, nm] : entries) {
    Group& g = groups[nm.text];
    g.locs.push_back(l);
    g.compound = nm.compound;
    ++g.count;
  }
  std::vector<std::pair<std::string, Group>> gv(groups.begin(), groups.end());
  std::sort(gv.begin(), gv.end(), [](const auto& a, const auto& b) {
    const std::string ka = alpha_key(a.first), kb = alpha_key(b.first);
    return ka != kb ? ka < kb : a.first < b.first;
  });
  std::string prefixes;
  for (std::size_t i = 0; i < gv.size(); ++i) {
    const auto& [t, g] = gv[i];
    std::string body = g.compound ? (g.count > 1 ? multiplier(g.count, true) : "") + enclose(t) : multiplier(g.count, false) + t;
    std::vector<std::string> shown;
    for (const auto& l : g.locs)
      if (!omit || (!l.empty() && l[0] == 'N')) shown.push_back(l);
    if (shown.empty() && i > 0 && !is_enclosed(body) && !(g.compound && g.count > 1)) body = enclose(body);
    append_piece(prefixes, shown.empty() ? body : join_locants(shown) + "-" + body);
  }
  return prefixes;
}

Score Namer::score(const Spec& sp, const Cand& c, const std::vector<int>& o, bool& valid) const {
  Score s;
  valid = true;
  const auto L = locants(c, o);
  if (sp.mode == Mode::Oyl && o[0] != sp.root) valid = false;
  std::vector<std::string> suf, un, ene, pre;
  for (int a : c.suffix_at) suf.push_back(L[a]);
  {
    std::string ind;
    std::vector<std::string> hyd;
    hydro(c, L, ind, hyd);
    if (!ind.empty()) s.ind.push_back(locant_key(ind));
    if (hyd.size() % 2) valid = false;  // hydro prefixes come in pairs; an odd set has no mancude reference
    un.insert(un.end(), hyd.begin(), hyd.end());
  }
  if (sp.mode == Mode::Yl || sp.mode == Mode::Ylidene || sp.mode == Mode::Oyl) suf.push_back(L[sp.root]);
  if (sp.mode == Mode::Carbonyl) suf.push_back(L[sp.link]);
  if (c.kind == PKind::Chain || c.kind == PKind::Carbocycle || c.kind == PKind::VonBaeyer) {
    const std::size_t m = o.size();
    const std::size_t lim = c.kind == PKind::Carbocycle ? m : m - 1;
    int counted = 0;
    for (std::size_t i = 0; i < lim; ++i) {
      const int b = order(o[i], o[(i + 1) % m]);
      if (b < 2) continue;
      if (i + 1 == m) valid = false;  // a ring multiple bond must not span the last and first atoms
      ++counted;
      un.push_back(L[o[i]]);
      if (b == 2) ene.push_back(L[o[i]]);
    }
    if (c.kind == PKind::VonBaeyer || (c.kind == PKind::Carbocycle && rs_[c.rsys].nhetero > 0)) {
      if (c.kind == PKind::VonBaeyer && counted != rs_[c.rsys].unsat) valid = false;  // would need a compound locant
      std::vector<std::string> het, by_kind;
      for (int z : {8, 16, 7})
        for (int x : o)
          if (at_[x].z == z) het.push_back(L[x]);
      s.het = sorted_keys(het);
      for (int z : {8, 16, 7}) {
        std::vector<std::string> l;
        for (int x : o)
          if (at_[x].z == z) l.push_back(L[x]);
        const auto k = sorted_keys(l);
        s.het.insert(s.het.end(), k.begin(), k.end());
      }
    }
  }
  std::vector<std::pair<std::string, std::string>> alpha;
  for (const auto& [p, nm] : c.pre)
    if (p >= 0) {
      pre.push_back(L[p]);
      alpha.push_back({alpha_key(nm.text) + "\x01" + nm.text, L[p]});
    }
  std::stable_sort(alpha.begin(), alpha.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  std::vector<double> al;
  for (std::size_t i = 0; i < alpha.size();) {
    std::size_t j = i;
    std::vector<double> ks;
    while (j < alpha.size() && alpha[j].first == alpha[i].first) ks.push_back(locant_key(alpha[j++].second));
    std::sort(ks.begin(), ks.end());
    al.insert(al.end(), ks.begin(), ks.end());
    i = j;
  }
  s.suf = sorted_keys(suf);
  s.unsat = sorted_keys(un);
  s.ene = sorted_keys(ene);
  s.pre = sorted_keys(pre);
  s.alpha = al;
  return s;
}

std::string Namer::hydride(const Cand& c, const std::vector<int>& o, bool omit_unsat, bool& saturated) const {
  saturated = true;
  if (c.kind == PKind::Benzene) return "benzene";
  if (c.kind == PKind::Template) {
    std::string ind;
    std::vector<std::string> hyd;
    hydro(c, locants(c, o), ind, hyd);
    std::string s;
    if (!hyd.empty()) s = join_locants(hyd) + "-" + multiplier(static_cast<int>(hyd.size()), false) + "hydro";
    if (!ind.empty()) s += (s.empty() ? "" : "-") + ind + "H-";
    append_piece(s, templates()[rs_[c.rsys].templ].name);
    return s;
  }
  const std::size_t m = o.size();
  const std::size_t lim = c.kind == PKind::Carbocycle ? m : m - 1;
  std::vector<std::string> enes, ynes;
  for (std::size_t i = 0; i < lim; ++i) {
    const int b = order(o[i], o[(i + 1) % m]);
    if (b == 2) enes.push_back(std::to_string(i + 1));
    if (b == 3) ynes.push_back(std::to_string(i + 1));
  }
  std::string stem = (c.kind == PKind::Carbocycle ? "cyclo" : "") + alkane_stem(static_cast<int>(m));
  if (c.kind == PKind::VonBaeyer || (c.kind == PKind::Carbocycle && rs_[c.rsys].nhetero > 0)) {
    const auto L = locants(c, o);
    std::string rep_prefixes;
    for (const auto& [z, word] : std::vector<std::pair<int, std::string>>{{8, "oxa"}, {16, "thia"}, {7, "aza"}}) {
      std::vector<std::string> l;
      for (int x : o)
        if (at_[x].z == z) l.push_back(L[x]);
      if (l.empty()) continue;
      rep_prefixes += (rep_prefixes.empty() ? "" : "-") + join_locants(l) + "-" + multiplier(static_cast<int>(l.size()), false) + word;
    }
    const RingSys& R = rs_[c.rsys];
    if (c.kind == PKind::VonBaeyer)
      stem = rep_prefixes + "bicyclo[" + std::to_string(R.vb[0]) + "." + std::to_string(R.vb[1]) + "." + std::to_string(R.vb[2]) + "]" + stem;
    else
      stem = rep_prefixes + stem;
  }
  if (enes.empty() && ynes.empty()) return stem + "ane";
  saturated = false;
  const bool a = enes.size() > 1 || (enes.empty() && ynes.size() > 1);
  std::string s = stem + (a ? "a" : "");
  if (!enes.empty()) s += (omit_unsat ? "" : "-" + join_locants(enes) + "-") + multiplier(static_cast<int>(enes.size()), false) + "ene";
  if (!ynes.empty()) {
    const std::string w = multiplier(static_cast<int>(ynes.size()), false) + "yne";
    if (!enes.empty() && starts_with_vowel(w)) s.pop_back();
    s += (omit_unsat ? "" : "-" + join_locants(ynes) + "-") + w;
  }
  return s;
}

std::optional<Built> Namer::assemble(const Spec& sp, const Cand& c, const std::vector<int>& o) {
  const bool mol = sp.mode == Mode::Molecule;
  const auto L = locants(c, o);
  const int m = static_cast<int>(o.size());
  int nnum = 0;
  for (const auto& pr : c.pre) nnum += pr.first >= 0;
  const int attachments = nnum + static_cast<int>(c.suffix_at.size()) + (mol ? 0 : 1);
  // P-14.3.4: locants omitted for one-carbon parents, monosubstituted ethane/ethene and
  // monosubstituted benzene and cycloalkanes
  bool omit = false;
  if (c.kind == PKind::Chain && m == 1) omit = true;
  if (c.kind == PKind::Chain && m == 2 && attachments <= 1) omit = true;
  const bool carbo = c.kind == PKind::Carbocycle && rs_[c.rsys].nhetero == 0;
  if ((c.kind == PKind::Benzene || (carbo && c.nmult == 0)) && attachments <= 1) omit = true;
  const bool omit_unsat = (c.kind == PKind::Chain && m <= 2) || (carbo && c.nmult == 1 && attachments == 0);

  Built out;
  // stereodescriptors of this level: parent atoms, parent double bonds, exocyclic double bonds
  std::vector<std::pair<std::string, std::string>> st;
  for (int a : o)
    if (auto it = stereo_.find(atom_key(a)); it != stereo_.end()) {
      st.push_back({L[a], to_string(it->second)});
      out.stereo.push_back(it->first);
    }
  for (int i = 0; i < m; ++i) {
    if (c.kind == PKind::Chain && i + 1 == m) break;
    if (m < 3 && c.kind != PKind::Chain) break;
    const int a = o[i], b = o[(i + 1) % m];
    if (auto it = stereo_.find(bond_key(a, b)); it != stereo_.end()) {
      st.push_back({L[a], to_string(it->second)});
      out.stereo.push_back(it->first);
    }
  }
  if (sp.mode == Mode::Ylidene)
    if (auto it = stereo_.find(bond_key(sp.root, sp.from)); it != stereo_.end()) {
      st.push_back({L[sp.root], to_string(it->second)});
      out.stereo.push_back(it->first);
    }
  for (const auto& [p, q] : c.exo_double)
    if (auto it = stereo_.find(bond_key(p, q)); it != stereo_.end()) {
      st.push_back({L[p], to_string(it->second)});
      out.stereo.push_back(it->first);
    }
  for (const auto& pr : c.pre) out.stereo.insert(out.stereo.end(), pr.second.stereo.begin(), pr.second.stereo.end());
  for (const auto& [o2, at] : c.ester_o) out.ester_o.push_back({o2, L[at]});

  std::vector<std::pair<std::string, Named>> entries;
  for (const auto& [p, nm] : c.pre) entries.push_back({p >= 0 ? L[p] : p == -1 ? "N" : "N'", nm});
  const std::string prefixes = format_prefixes(entries, omit);

  // parent hydride + suffix
  bool saturated = true;
  const std::string hyd = hydride(c, o, omit_unsat, saturated);
  std::string core;
  const int k = static_cast<int>(c.suffix_at.size());
  auto suffix_locs = [&]() {
    if (omit) return std::string();
    std::vector<std::string> l;
    for (int a : c.suffix_at) l.push_back(L[a]);
    return join_locants(l);
  };
  if (mol) {
    if (k == 0 || sp.principal == Cls::None) {
      core = hyd;
    } else {
      std::string word;
      switch (sp.principal) {
        case Cls::Acid: word = c.inchain ? "oic acid" : "carboxylic acid"; break;
        case Cls::SulfonicAcid: word = "sulfonic acid"; break;
        case Cls::Sulfonamide: word = "sulfonamide"; break;
        case Cls::Ester: word = c.inchain ? "oate" : "carboxylate"; break;
        case Cls::Amide: word = c.inchain ? "amide" : "carboxamide"; break;
        case Cls::Nitrile: word = c.inchain ? "nitrile" : "carbonitrile"; break;
        case Cls::Aldehyde: word = c.inchain ? "al" : "carbaldehyde"; break;
        case Cls::Ketone: word = "one"; break;
        case Cls::Thione: word = "thione"; break;
        case Cls::Thioamide: word = c.inchain ? "thioamide" : "carbothioamide"; break;
        case Cls::Hydrazide: word = c.inchain ? "ohydrazide" : "carbohydrazide"; break;
        case Cls::Alcohol: word = "ol"; break;
        case Cls::Thiol: word = "thiol"; break;
        case Cls::Amine: word = "amine"; break;
        case Cls::Carboxylate: word = c.inchain ? "oate" : "carboxylate"; break;
        case Cls::Sulfonate: word = "sulfonate"; break;
        case Cls::Olate: word = "olate"; break;
        case Cls::Ammonium: word = "aminium"; break;
        case Cls::RingCation: word = "ium"; break;
        case Cls::Imine: word = "imine"; break;
        default: break;
      }
      core = attach(hyd, c.inchain ? "" : suffix_locs(), multiplier(k, false) + word);
      // retained names (P-65.1.1.1, P-66.1.1.1, P-66.5.1.1, P-66.6.1.1, P-63.1.1, P-62.2.1.1)
      if (c.kind == PKind::Chain && c.inchain && k == 1 && m <= 2) {
        static const std::map<Cls, std::pair<const char*, const char*>> chain_names = {
            {Cls::Acid, {"formic acid", "acetic acid"}},     {Cls::Ester, {"formate", "acetate"}}, {Cls::Carboxylate, {"formate", "acetate"}},
            {Cls::Amide, {"formamide", "acetamide"}},        {Cls::Hydrazide, {"formohydrazide", "acetohydrazide"}},        {Cls::Aldehyde, {"formaldehyde", "acetaldehyde"}},
            {Cls::Nitrile, {"formonitrile", "acetonitrile"}}};
        if (auto it = chain_names.find(sp.principal); it != chain_names.end()) core = m == 1 ? it->second.first : it->second.second;
      }
      if (c.kind == PKind::Benzene && k == 1) {
        static const std::map<Cls, const char*> benz = {
            {Cls::Acid, "benzoic acid"},  {Cls::Ester, "benzoate"},        {Cls::Amide, "benzamide"}, {Cls::Nitrile, "benzonitrile"},
            {Cls::Aldehyde, "benzaldehyde"}, {Cls::Alcohol, "phenol"}, {Cls::Amine, "aniline"}, {Cls::Carboxylate, "benzoate"},
            {Cls::Olate, "phenolate"}, {Cls::Ammonium, "anilinium"}, {Cls::Hydrazide, "benzohydrazide"}};
        if (auto it = benz.find(sp.principal); it != benz.end()) core = it->second;
      }
    }
  } else {
    const std::string fv = sp.mode == Mode::Carbonyl ? L[sp.link] : L[sp.root];
    switch (sp.mode) {
      case Mode::Yl:
      case Mode::Ylidene: {
        const std::string w = sp.mode == Mode::Yl ? "yl" : "ylidene";
        if (c.kind == PKind::Benzene) {
          if (sp.mode != Mode::Yl) {
            fail("benzene ylidene");
            return std::nullopt;
          }
          core = "phenyl";
        } else if ((c.kind == PKind::Chain || carbo) && saturated && fv == "1") {
          core = (c.kind == PKind::Carbocycle ? "cyclo" : "") + alkane_stem(m) + w;
        } else {
          core = attach(hyd, omit ? "" : fv, w);
        }
        break;
      }
      case Mode::Oyl: core = m == 1 ? "formyl" : m == 2 ? "acetyl" : attach(hyd, "", "oyl"); break;
      case Mode::Carbonyl: core = c.kind == PKind::Benzene && sp.word == "carbonyl" ? "benzoyl" : attach(hyd, omit ? "" : fv, sp.word); break;
      default: break;
    }
  }
  std::string name = prefixes;
  append_piece(name, core);
  if (!st.empty()) {
    std::stable_sort(st.begin(), st.end(), [](const auto& a, const auto& b) { return locant_key(a.first) < locant_key(b.first); });
    std::string d;
    for (std::size_t i = 0; i < st.size(); ++i) d += (i ? "," : "") + (omit ? "" : st[i].first) + st[i].second;
    name = "(" + d + ")-" + name;
  }
  out.name = name;
  out.compound = !c.pre.empty() || has_digit(name) || !st.empty();
  return out;
}

std::optional<Built> Namer::build(const Spec& sp) {
  const bool mol = sp.mode == Mode::Molecule;
  std::vector<Cand> cands;
  int ring_of_root = -1;
  if (sp.mode == Mode::Yl || sp.mode == Mode::Ylidene) ring_of_root = ringsys_[sp.root];
  if (sp.mode == Mode::Carbonyl) ring_of_root = ringsys_[sp.link];
  for (int i = 0; i < static_cast<int>(rs_.size()); ++i) {
    const RingSys& R = rs_[i];
    if (!sp.region[R.atoms[0]]) continue;
    if (!mol && i != ring_of_root) continue;
    Cand c;
    c.kind = R.kind;
    c.rsys = i;
    c.atoms = R.atoms;
    c.orders = R.orders;
    cands.push_back(std::move(c));
  }
  // chains (P-44.3): paths in the forest of eligible acyclic carbons
  const bool chains = mol || sp.mode == Mode::Oyl || (sp.mode != Mode::Molecule && ring_of_root < 0);
  if (chains) {
    auto PA = [&](int a) { return mol && is_acyl(sp.principal) && acyl_[a] == sp.principal && sp.pok[a]; };
    std::vector<char> E(n(), 0);
    for (int a = 0; a < n(); ++a) {
      if (!sp.region[a] || at_[a].z != 6 || at_[a].ring) continue;
      const Cls k = acyl_[a];
      if ((k == Cls::Acid || k == Cls::Carboxylate || k == Cls::Nitrile || k == Cls::Carbonate || k == Cls::Carbamate || k == Cls::Urea || k == Cls::Thiourea) && !PA(a))
        continue;
      E[a] = 1;
    }
    if (sp.mode == Mode::Oyl) E[sp.root] = 1;
    std::vector<int> ends, normal_deg(n(), 0);
    for (int a = 0; a < n(); ++a) {
      if (!E[a]) continue;
      for (int b : at_[a].nb)
        if (E[b] && !PA(b)) ++normal_deg[a];
      if (PA(a) || normal_deg[a] <= 1 || (sp.mode == Mode::Oyl && a == sp.root)) ends.push_back(a);
    }
    for (std::size_t iu = 0; iu < ends.size(); ++iu) {
      const int u = ends[iu];
      std::vector<int> parent(n(), -2);
      std::vector<int> queue = {u};
      parent[u] = -1;
      for (std::size_t h = 0; h < queue.size(); ++h) {
        const int a = queue[h];
        if (a != u && PA(a)) continue;  // principal acyl carbons only end chains
        for (int b : at_[a].nb)
          if (E[b] && parent[b] == -2) parent[b] = a, queue.push_back(b);
      }
      for (std::size_t iv = iu; iv < ends.size(); ++iv) {
        const int v = ends[iv];
        if (parent[v] == -2) continue;
        if (u == v && normal_deg[u] > 0) continue;
        std::vector<int> path;
        for (int x = v; x != -1; x = parent[x]) path.push_back(x);
        std::reverse(path.begin(), path.end());
        bool ok = true;
        for (std::size_t i = 1; i + 1 < path.size(); ++i) ok &= !PA(path[i]);
        if (sp.mode == Mode::Yl || sp.mode == Mode::Ylidene) ok &= std::find(path.begin(), path.end(), sp.root) != path.end();
        if (sp.mode == Mode::Carbonyl) ok &= std::find(path.begin(), path.end(), sp.link) != path.end();
        if (sp.mode == Mode::Oyl) ok &= path.front() == sp.root || path.back() == sp.root;
        for (int a : path)
          if (PA(a))
            for (int b : at_[a].nb)
              if (at_[b].z == 6 && std::find(path.begin(), path.end(), b) == path.end()) ok = false;
        if (!ok) continue;
        Cand c;
        c.kind = PKind::Chain;
        c.atoms = path;
        c.orders.push_back(path);
        if (path.size() > 1) c.orders.push_back(std::vector<int>(path.rbegin(), path.rend()));
        cands.push_back(std::move(c));
      }
    }
  }
  if (cands.empty()) {
    fail("no parent hydride (no carbon skeleton)");
    return std::nullopt;
  }
  if (cands.size() > 20000) {  // bounds the work on very large, highly branched skeletons
    fail("too many parent chain candidates");
    return std::nullopt;
  }

  struct Ranked {
    std::vector<double> senior;  // larger is better
    Score sc;
    int npre = 0;
    Built b;
  };
  std::optional<Ranked> best;
  for (Cand& c : cands) {
    if (!analyze(sp, c)) continue;
    std::optional<std::pair<Score, std::size_t>> bo;
    for (std::size_t i = 0; i < c.orders.size(); ++i) {
      bool valid = true;
      Score s = score(sp, c, c.orders[i], valid);
      if (!valid) continue;
      if (!bo || compare_score(s, bo->first) < 0) bo = std::make_pair(s, i);
    }
    if (!bo) continue;
    auto b = assemble(sp, c, c.orders[bo->second]);
    if (!b) continue;
    Ranked r;
    r.sc = bo->first;
    r.b = *b;
    for (const auto& pr : c.pre) r.npre += pr.first >= 0;
    const bool ring = c.kind != PKind::Chain;
    r.senior.push_back(mol ? static_cast<double>(c.suffix_at.size()) : 0.0);  // P-44.1.1
    r.senior.push_back(ring ? 1 : 0);                                            // P-44.1.2.2
    if (ring) {                                                                  // P-44.2
      const RingSys& R = rs_[c.rsys];
      r.senior.insert(r.senior.end(), {double(R.nhetero > 0), double(R.nN > 0), double(R.nrings), double(R.atoms.size()),
                                       double(R.nhetero), double(R.variety), double(R.nO), double(R.nS), double(R.nN), double(R.unsat), 0, 0, 0});
    } else {  // P-44.3
      r.senior.insert(r.senior.end(), {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, double(c.atoms.size()), double(c.nmult), double(c.ndouble)});
    }
    bool better = !best;
    if (!better) {
      int cmp = 0;
      for (std::size_t i = 0; i < r.senior.size() && !cmp; ++i)
        if (r.senior[i] != best->senior[i]) cmp = r.senior[i] > best->senior[i] ? -1 : 1;
      if (!cmp) cmp = compare_keys(r.sc.suf, best->sc.suf);
      if (!cmp) cmp = compare_keys(r.sc.unsat, best->sc.unsat);
      if (!cmp) cmp = compare_keys(r.sc.ene, best->sc.ene);
      if (!cmp && r.npre != best->npre) cmp = r.npre > best->npre ? -1 : 1;  // P-45.2.1
      if (!cmp) cmp = compare_score(r.sc, best->sc);                         // P-45.2.2, P-45.2.3
      if (!cmp) cmp = r.b.name < best->b.name ? -1 : (r.b.name > best->b.name ? 1 : 0);
      better = cmp < 0;
    }
    if (better) best = std::move(r);
  }
  if (!best) return std::nullopt;
  return best->b;
}

std::optional<std::string> Namer::run() {
  if (!load() || !perceive() || !rings()) return std::nullopt;
  Cls principal = Cls::None;
  auto senior = [&](Cls c) {
    if (static_cast<int>(c) < static_cast<int>(principal)) principal = c;
  };
  for (int a = 0; a < n(); ++a) {
    if (acyl_[a] != Cls::None && static_cast<int>(acyl_[a]) < 50) senior(acyl_[a]);
    if (is_OH(a) && at_[at_[a].nb[0]].z == 6 && acyl_[at_[a].nb[0]] != Cls::Acid) senior(Cls::Alcohol);
    if (is_SH(a)) senior(Cls::Thiol);
    if (amineN_[a]) senior(Cls::Amine);
    if (cationN_[a]) senior(Cls::Ammonium);
    if (ringCat_[a]) senior(Cls::RingCation);
    if (anionO_[a] && at_[at_[a].nb[0]].z == 6 && acyl_[at_[a].nb[0]] != Cls::Carboxylate) senior(Cls::Olate);
    if (imineC_[a]) senior(Cls::Imine);
  }
  std::optional<std::string> result;
  std::vector<Key> cited;
  auto finish = [&](const Built& b, const std::string& name, const std::vector<Key>& extra) {
    if (!result || name.size() < result->size() || (name.size() == result->size() && name < *result)) {
      result = name;
      cited = b.stereo;
      cited.insert(cited.end(), extra.begin(), extra.end());
    }
  };
  // the alcohol / organyl component of a functional class name ("ethyl", "[(2R)-butan-2-yl]")
  auto organyl = [](const Named& R) { return R.text[0] == '(' ? enclose(R.text) : R.text; };
  if (principal == Cls::PhosMono || principal == Cls::Phosphonic) {
    // "R dihydrogen phosphate" (P-67.1.3.2) and "R-phosphonic acid" (P-67.1.1.1)
    for (int P = 0; P < n(); ++P) {
      if (acyl_[P] != principal) continue;
      int c = -1, o = -1;
      for (int x : at_[P].nb) {
        if (at_[x].z == 6) c = x, o = P;
        if (at_[x].z == 8 && deg(x) == 2) o = x, c = at_[x].nb[0] == P ? at_[x].nb[1] : at_[x].nb[0];
      }
      auto R = sub(c, o, 1);
      if (!R) continue;
      const std::string name = principal == Cls::PhosMono ? organyl(*R) + " dihydrogen phosphate"
                                                           : (R->compound ? enclose(R->text) : R->text) + "phosphonic acid";
      finish(Built{}, name, R->stereo);
    }
    if (!result) {
      if (err_.empty()) err_ = "phosphorus acid derivative could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::Carbonate || principal == Cls::Carbamate || principal == Cls::Urea || principal == Cls::Thiourea) {
    // functional class names of carbonic acid derivatives (P-65.2.1): "dimethyl carbonate",
    // "ethyl N-methylcarbamate", "N-methyl-N'-phenylurea"
    for (int C = 0; C < n(); ++C) {
      if (acyl_[C] != principal) continue;
      std::vector<int> os, ns;
      for (int x : at_[C].nb) {
        if (at_[x].z == 8 && deg(x) == 2) os.push_back(x);
        if (at_[x].z == 7) ns.push_back(x);
      }
      std::vector<Key> extra;
      std::vector<std::string> organyls;
      bool ok = true;
      for (int o : os) {
        const int r = at_[o].nb[0] == C ? at_[o].nb[1] : at_[o].nb[0];
        auto R = sub(r, o, 1);
        if (!R) { ok = false; break; }
        organyls.push_back(organyl(*R));
        extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
      }
      // N substituents, with N / N' locants (unprimed to the nitrogen giving lower locants, P-31.1.4.3.4)
      std::vector<std::vector<Named>> nsubs;
      for (int N : ns) {
        nsubs.emplace_back();
        for (int x : at_[N].nb) {
          if (x == C) continue;
          auto R = sub(x, N, 1);
          if (!R) { ok = false; break; }
          extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
          nsubs.back().push_back(*R);
        }
      }
      if (!ok) continue;
      std::string npart;
      if (principal == Cls::Urea || principal == Cls::Thiourea) {
        int first = 0;
        if (nsubs[1].size() > nsubs[0].size()) first = 1;
        else if (nsubs[1].size() == nsubs[0].size() && !nsubs[0].empty()) {
          auto least = [](const std::vector<Named>& v) {
            std::string m = "~";
            for (const auto& x : v) m = std::min(m, alpha_key(x.text));
            return m;
          };
          if (least(nsubs[1]) < least(nsubs[0])) first = 1;
        }
        std::vector<std::pair<std::string, Named>> e;
        for (const auto& nm : nsubs[first]) e.push_back({"N", nm});
        for (const auto& nm : nsubs[1 - first]) e.push_back({"N'", nm});
        const std::string best = format_prefixes(e, false);
        npart = best;
        finish(Built{}, npart + (principal == Cls::Urea ? "urea" : "thiourea"), extra);
        continue;
      }
      if (principal == Cls::Carbamate) {
        std::vector<std::pair<std::string, Named>> e;
        for (const auto& nm : nsubs[0]) e.push_back({"N", nm});
        npart = format_prefixes(e, false);
      }
      std::sort(organyls.begin(), organyls.end(), [](const std::string& a, const std::string& b) { return alpha_key(a) < alpha_key(b); });
      std::string alk;
      if (organyls.size() == 2 && organyls[0] == organyls[1]) {
        alk = (organyls[0][0] == '(' || organyls[0][0] == '[' || has_digit(organyls[0])) ? "bis" + enclose(organyls[0]) : "di" + organyls[0];
      } else {
        for (const auto& o : organyls) alk += (alk.empty() ? "" : " ") + o;
      }
      finish(Built{}, alk + " " + npart + (principal == Cls::Carbonate ? "carbonate" : "carbamate"), extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "carbonic acid derivative could not be named";
      return std::nullopt;
    }
  } else if (principal != Cls::Ester) {
    Spec sp;
    sp.region.assign(n(), 1);
    sp.pok.assign(n(), 1);
    sp.principal = principal;
    auto b = build(sp);
    if (!b) {
      if (err_.empty()) err_ = "no parent could be named";
      return std::nullopt;
    }
    finish(*b, b->name, {});
  } else {
    // esters: functional class names "alkyl ...oate" (P-65.6.3.2)
    std::vector<int> E;
    for (int a = 0; a < n(); ++a)
      if (acyl_[a] == Cls::Ester) E.push_back(a);
    std::vector<std::vector<int>> attempts;
    {
      std::vector<char> blocked(n(), 0);
      for (int a : E) blocked[ester_O(a)] = 1;
      const auto reg = region_from(E[0], blocked);
      bool all = true;
      for (int a : E) all &= reg[a] != 0;
      if (all) attempts.push_back(E);
      else
        for (int a : E) attempts.push_back({a});
    }
    for (const auto& att : attempts) {
      std::vector<char> blocked(n(), 0);
      Spec sp;
      sp.pok.assign(n(), 0);
      for (int a : att) {
        blocked[ester_O(a)] = 1;
        sp.pok[a] = 1;
      }
      sp.region = region_from(att[0], blocked);
      sp.principal = Cls::Ester;
      auto b = build(sp);
      if (!b || b->ester_o.empty()) continue;
      std::map<std::string, std::pair<int, bool>> alk;
      std::map<std::string, std::vector<std::string>> alk_locs;
      std::vector<Key> extra;
      bool ok = true;
      for (const auto& [o, loc] : b->ester_o) {
        const int x0 = at_[o].nb[0], x1 = at_[o].nb[1];
        const int r = (acyl_[x0] == Cls::Ester && sp.pok[x0]) ? x1 : x0;
        auto R = sub(r, o, 1);
        if (!R) {
          ok = false;
          break;
        }
        auto& g = alk[organyl(*R)];
        ++g.first;
        g.second = R->compound;
        alk_locs[organyl(*R)].push_back(loc);
        extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
      }
      if (!ok) continue;
      std::vector<std::pair<std::string, std::pair<int, bool>>> v(alk.begin(), alk.end());
      std::sort(v.begin(), v.end(), [](const auto& a, const auto& b2) { return alpha_key(a.first) < alpha_key(b2.first); });
      // different organyl groups on a polyester: each is located on its acid position (P-65.6.3.3.3)
      const bool locate = v.size() > 1 && b->ester_o.size() > 1;
      std::string s;
      for (const auto& [t, g] : v) {
        if (!s.empty()) s += " ";
        const bool comp = g.second || t[0] == '(' || t[0] == '[';
        if (locate) {
          const std::string body = g.first > 1 ? multiplier(g.first, comp) + (comp ? enclose(t) : t) : (comp ? enclose(t) : t);
          s += join_locants(alk_locs[t]) + "-" + body;
        } else {
          s += g.first > 1 ? (comp ? multiplier(g.first, true) + enclose(t) : multiplier(g.first, false) + t) : t;
        }
      }
      finish(*b, s + " " + b->name, extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "ester could not be named";
      return std::nullopt;
    }
  }
  const std::set<Key> c(cited.begin(), cited.end());
  for (const auto& kv : stereo_)
    if (!c.count(kv.first)) {
      err_ = "a stereo element could not be placed in the name";
      return std::nullopt;
    }
  return result;
}

}  // namespace

namespace {

// Salts, hydrates and acid adducts (P-72, P-74): each component is named on its own and the names
// are combined: cations, then anions, each set in alphanumerical order ("disodium butanedioate",
// "1-methylpyridin-1-ium iodide"), or "name hydrochloride" for a neutral base with hydrogen halide.
std::optional<std::string> name_components(const ChemicalSystemState& st, std::uint32_t sp, std::string& err) {
  const Species& s = st.species[sp];
  std::vector<int> comp(s.atom_count, -1);
  int ncomp = 0;
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    if (comp[i] >= 0) continue;
    std::vector<std::uint32_t> stack = {s.atom_begin + i};
    comp[i] = ncomp;
    while (!stack.empty()) {
      const std::uint32_t g = stack.back();
      stack.pop_back();
      for (std::uint32_t b : st.adjacency.neighbors(g))
        if (comp[b - s.atom_begin] < 0) comp[b - s.atom_begin] = ncomp, stack.push_back(b);
    }
    ++ncomp;
  }
  if (ncomp == 1) {
    Namer nm(st, sp);
    auto name = nm.run();
    if (!name) err = nm.error();
    return name;
  }
  struct Ion {
    std::string name;
    int charge = 0;
    bool compound = false;
  };
  std::vector<Ion> cations, anions, neutrals;
  std::map<std::string, int> adducts;
  int waters = 0;
  for (int c = 0; c < ncomp; ++c) {
    std::vector<std::uint32_t> atoms;
    int heavy = -1, nheavy = 0, charge = 0, hydrogens = 0;
    for (std::uint32_t i = 0; i < s.atom_count; ++i) {
      if (comp[i] != c) continue;
      const std::uint32_t g = s.atom_begin + i;
      atoms.push_back(g);
      charge += st.atoms.formal_charge[g];
      if (st.atoms.element[g] == 1) ++hydrogens;
      else ++nheavy, heavy = static_cast<int>(g), hydrogens += st.atoms.implicit_h[g];
    }
    if (nheavy == 1) {
      const int z = st.atoms.element[static_cast<std::uint32_t>(heavy)];
      static const std::map<std::pair<int, int>, const char*> cat = {
          {{3, 1}, "lithium"}, {{11, 1}, "sodium"},   {{19, 1}, "potassium"}, {{37, 1}, "rubidium"}, {{55, 1}, "caesium"},
          {{12, 2}, "magnesium"}, {{20, 2}, "calcium"}, {{38, 2}, "strontium"}, {{56, 2}, "barium"}, {{30, 2}, "zinc"}};
      static const std::map<int, std::pair<const char*, const char*>> hal = {
          {9, {"fluoride", "hydrofluoride"}}, {17, {"chloride", "hydrochloride"}}, {35, {"bromide", "hydrobromide"}}, {53, {"iodide", "hydroiodide"}}};
      if (auto it = cat.find({z, charge}); it != cat.end() && hydrogens == 0) {
        cations.push_back({it->second, charge, false});
        continue;
      }
      if (auto it = hal.find(z); it != hal.end()) {
        if (charge == -1 && hydrogens == 0) {
          anions.push_back({it->second.first, -1, false});
          continue;
        }
        if (charge == 0 && hydrogens == 1) {
          ++adducts[it->second.second];
          continue;
        }
      }
      if (z == 8 && charge == 0 && hydrogens == 2) {
        ++waters;
        continue;
      }
      if (z == 7 && charge == 1 && hydrogens == 4) {
        cations.push_back({"azanium", 1, false});
        continue;
      }
      if (z != 6) {
        err = std::string("inorganic component ") + std::string(element_symbol(static_cast<std::uint8_t>(z))) + " is not covered";
        return std::nullopt;
      }
    }
    Namer nm(st, sp, atoms);
    auto name = nm.run();
    if (!name) {
      err = nm.error();
      return std::nullopt;
    }
    Ion ion{*name, charge, has_digit(*name) || name->find(' ') != std::string::npos || name->find('(') != std::string::npos};
    (charge > 0 ? cations : charge < 0 ? anions : neutrals).push_back(ion);
  }
  int total = 0;
  for (const auto& x : cations) total += x.charge;
  for (const auto& x : anions) total += x.charge;
  for (const auto& [n, k] : adducts) total += 0 * k;
  if (total != 0) {
    err = "components do not balance in charge";
    return std::nullopt;
  }
  std::string out;
  if (!neutrals.empty()) {
    if (neutrals.size() != 1 || !cations.empty() || !anions.empty() || (adducts.empty() && waters == 0)) {
      err = "mixture of several neutral components";
      return std::nullopt;
    }
    out = neutrals[0].name;
    for (const auto& [n, k] : adducts) out += " " + multiplier(k, false) + n;
  } else {
    if (!adducts.empty()) {
      err = "hydrogen halide alongside a salt";
      return std::nullopt;
    }
    auto group = [](std::vector<Ion> v) {
      std::map<std::string, std::pair<int, bool>> g;
      for (const auto& x : v) ++g[x.name].first, g[x.name].second = x.compound;
      std::vector<std::pair<std::string, std::pair<int, bool>>> gv(g.begin(), g.end());
      std::sort(gv.begin(), gv.end(), [](const auto& a, const auto& b) { return alpha_key(a.first) < alpha_key(b.first); });
      std::string r;
      for (const auto& [n, kc] : gv) {
        if (!r.empty()) r += " ";
        if (kc.first == 1) r += n;
        else r += kc.second ? multiplier(kc.first, true) + enclose(n) : multiplier(kc.first, false) + n;
      }
      return r;
    };
    if (cations.empty() || anions.empty()) {
      err = "unbalanced ionic components";
      return std::nullopt;
    }
    out = group(cations) + " " + group(anions);
  }
  if (waters) out += " " + multiplier(waters, false) + "hydrate";
  return out;
}

}  // namespace

NameResult iupac_name(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag) {
  NameResult r;
  if (species >= state.species.size()) {
    diag.error(DiagCode::NamingUnsupported, "no IUPAC name: no such species");
    return r;
  }
  std::string err;
  auto name = name_components(state, species, err);
  if (!name) {
    diag.error(DiagCode::NamingUnsupported, "no IUPAC name: " + (err.empty() ? std::string("unsupported structure") : err)).on_species(species);
    return r;
  }
  r.ok = true;
  r.name = *name;
  return r;
}

NameResult iupac_name_from_smiles(const std::string& smiles, Diagnostics& diag) {
  ChemicalSystemState st;
  const std::int32_t sp = parse_smiles(st, smiles, diag);
  if (sp < 0) return {};
  return iupac_name(st, static_cast<std::uint32_t>(sp), diag);
}

}  // namespace chemlarp
