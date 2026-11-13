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
#include <cmath>
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
#include "chemlarp/rings.hpp"
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
  v += 0.0001 * static_cast<double>(std::count(l.begin(), l.end(), '\''));  // 2 < 2' < 2''
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
  if (n <= 20) return simple[n];
  if (n >= 100) return "";
  // P-14.2.1.2: units before tens (henicosa, docosa, tricosa, ... hentriaconta, dotriaconta, ...)
  static const char* units[] = {"", "hen", "do", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona"};
  static const char* tens[] = {"", "", "cosa", "triaconta", "tetraconta", "pentaconta", "hexaconta", "heptaconta", "octaconta", "nonaconta"};
  if (n % 10 == 0) return n == 20 ? "icosa" : tens[n / 10];
  return std::string(units[n % 10]) + (n / 10 == 2 ? "cosa" : tens[n / 10]);
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
      if (++run > 4) return false;
    } else if ((c >= '0' && c <= '9') || c == ',' || c == '.' || c == '-' || c == '\'' || c == ':' || c == '^' || c == '{' || c == '}') {
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
      {"quinolizine", "c1cccn2ccccc12", {"1", "2", "3", "4", "5", "6", "7", "8", "9", "9a"}, true},
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
      std::string p = t.mancude ? generic_smarts(t.pattern) : std::string(t.pattern);
      if (!t.mancude)  // saturated parents also carry ring cations: piperidin-1-ium
        for (std::size_t k = p.find("[NX3]"); k != std::string::npos; k = p.find("[NX3]", k + 1)) p.replace(k, 5, "[N;X3+0,X4+,X2-]");
      if (!t.mancude)  // and ring sulfur oxides: 1,1-dioxo-1lambda6-thiolane
        for (std::size_t k = p.find("[SX2]"); k != std::string::npos; k = p.find("[SX2]", k + 1)) p.replace(k, 5, "[S;X2,X3,X4]");
      v.push_back(parse_smarts(p, d));
    }
    return v;
  }();
  return q;
}

// ---- data ------------------------------------------------------------------------------------

// Seniority of classes (P-41), lower = more senior.
// Sulfonyl and Sulfinyl mark sulfones/sulfoxides (never suffixes).
enum class Cls : int {
  Carboxylate = -60, CarbonicAnion = -55, Azanide = -44, Borate = -43, Nitrate = -42, SulfateAnion = -48, SulfamateAnion = -47, SulfuricInorg = -46, PhosAnion = -45, Sulfonate = -50, RingAnion = -41, Olate = -40, Ammonium = -30, Sulfonium = -25, Phosphonium = -24, Diazonium = -22, RingCation = -20,  // anions, then cations (P-41: ions first)
  Acid = 10, CarbonicAcid = 15, IsothioAcid = 16, SulfonicAcid = 20, SulfamicAcid = 21, SulfinicAcid = 22, Phosphonic = 30, PhosphinicAcid = 31, Anhydride = 35, Ester = 40, ThioEster = 42,
  SulfEster = 44, Isothiourea = 45, SulfamateEster = 46, Carbonate = 50, Carbamate = 60, PhosAcidInorg = 69, PhosMono = 70, PolyP = 71, PhosEster = 72, SulfateMono = 73, AcylHalide = 74,
  PhosThio = 75, SulfonylHalide = 76, Amide = 80, Semicarbazide = 85, Thioamide = 90, ThioSemicarbazide = 95, Hydrazide = 100, Amidine = 110, Sulfonamide = 120, SulfonoHydrazide = 122, SulfuricDiamide = 125, Urea = 130,
  Thiourea = 140, Guanidine = 150, Nitrile = 160, Aldehyde = 170, Ketone = 180, PhosOxide = 181, Thione = 190, Alcohol = 200,
  Thiol = 210, Peroxol = 215, Amine = 220, Imine = 230, Hydrazine = 240, Phosphane = 400, Sulfonyl = 500, Sulfinyl = 510, AcylSub = 520, RingP = 530, None = 999
};

bool is_acyl(Cls c) {
  return c == Cls::Acid || c == Cls::Carboxylate || c == Cls::Ester || c == Cls::ThioEster || c == Cls::Amide || c == Cls::Thioamide || c == Cls::Hydrazide ||
         c == Cls::Amidine ||
         c == Cls::Nitrile ||
         c == Cls::Aldehyde;
}
bool is_sulfur_suffix(Cls c) {
  return c == Cls::SulfonicAcid || c == Cls::SulfinicAcid || c == Cls::Sulfonate || c == Cls::SulfEster || c == Cls::Sulfonamide ||
         c == Cls::SulfonoHydrazide;
}

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
  bool spiro = false;           // VonBaeyer kind used for a monospiro system: spiro[vb0.vb1]
  std::vector<std::string> names;               // Template kind: parent name per numbering in `orders`
  std::vector<std::vector<std::string>> locs;   // Template kind: locants per numbering
  bool het_numbering = false;                   // numbering chosen by heteroatom locants first (Hantzsch-Widman)
  bool spiro_mancude = false;                   // saturation already worked out per component (spiro_fused)
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

struct Component;

// A fused ring system: rings, shared bonds, periphery.
struct FusedSys {
  struct Edge {
    int r, s, u, v;
  };
  std::vector<int> atoms;
  std::vector<std::vector<int>> rings;
  std::map<int, std::vector<int>> rings_of;
  std::vector<Edge> edges;
  std::vector<std::vector<int>> radj;
  std::map<std::pair<int, int>, int> bond_rings;  // rings containing each bond
  std::set<int> interior;                         // atoms with no peripheral bond
};

struct SpiroPart {
  std::vector<int> atoms;
  std::vector<std::vector<int>> orders;
  std::vector<std::vector<std::string>> locs;
  std::vector<std::string> names;
  bool mancude = false;
};

struct FusedNumbering {
  std::vector<int> o;          // peripheral atoms clockwise from locant 1, then interior atoms
  std::vector<std::string> L;  // their locants
  std::vector<double> key;
};

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
  std::vector<char> hydrN_, hydrN2_, hydroxamO_;
  std::vector<char> amidN_;  // imino nitrogen of an amidine or guanidine
  std::vector<char> lambdaS_;   // ring sulfur with 1 or 2 oxo groups (lambda4 / lambda6); ring P(=X) (5: lambda5)
  std::string lambda_tag(int a) const { return lambdaS_[a] == 5 ? "λ5" : lambdaS_[a] == 1 ? "λ4" : "λ6"; }
  std::vector<char> thioP_;     // sulfur of a thiophosphate (P=S or P-S-C)
  std::vector<int> hzPartner_;  // hydrazine / hydrazone nitrogens: the other nitrogen, else -1
  std::vector<int> azoPartner_;  // R-N=N-R': the other nitrogen, else -1
  std::vector<char> nO_;        // oxygen of a hydroxylamine (N-O-H or N-O-C, N not acyl)
  // prefix-only groups (P-61.5, P-63.4.2): the anchor atom bonded to the parent holds the kind,
  // the other atoms of the group are flagged in pseudoIn_
  enum Pseudo : char { PNone, Nitroso, Azido, Isocyanato, Isothiocyanato, Cyanato, Thiocyanato, Hydroperoxy, Peroxy };
  std::vector<char> pseudo_, pseudoIn_;
  // semicarbazides R2N-C(=O/S)-NR-NR2: the carbonyl carbon; acyl groups on their hydrazine nitrogens
  std::vector<char> semC_, acylSub_;
  // centres of simple ions named by ligands or suffixes: diazonium N, sulfonium S, phosphonium P,
  // borate B, sulfonyl azanide N, nitrate N (Cls value stored); their terminal atoms go in pseudoIn_
  std::vector<Cls> ionC_;
  std::vector<char> ringAn_;  // ring N- or C- with two single ring bonds: "-ide" (P-72.2.2.2)
  std::vector<char> carbX_;  // O and S atoms on a carbonic or carbamic acid carbon (acid, ester or anion)
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
  bool plain_acyl(int y) const {  // acyclic C(=O) whose only other heteroatom is the one it is reached from
    if (at_[y].z != 6 || at_[y].ring || dO(y) < 0) return false;
    int het = 0;
    for (int w : at_[y].nb) het += at_[w].z != 6 && w != dO(y);
    return het == 1;
  }
  bool is_Ominus(int o) const { return at_[o].z == 8 && at_[o].q == -1 && deg(o) == 1 && at_[o].h == 0; }
  bool is_OH(int o) const { return at_[o].z == 8 && at_[o].h == 1 && deg(o) == 1 && at_[o].q == 0; }
  bool is_SH(int s) const { return at_[s].z == 16 && at_[s].h == 1 && deg(s) == 1; }
  int ester_O(int c) const {  // the ester heteroatom of a carboxylic, thio or sulfonic ester
    for (int x : at_[c].nb)
      if ((at_[x].z == 8 || (acyl_[c] == Cls::ThioEster && at_[x].z == 16)) && deg(x) == 2) return x;
    return -1;
  }
  std::optional<Named> phosphoryl_name(int P, int from);
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
  bool von_baeyer_poly(RingSys& R);
  std::vector<std::pair<std::string, int>> ring_multiple_bonds(const Cand& c, const std::vector<int>& o) const;
  bool spiro_system(RingSys& R);
  bool hantzsch_widman(RingSys& R, bool& mancude);
  bool fusion_bicycle(RingSys& R);
  bool fusion_polycycle(RingSys& R);
  bool spiro_fused(RingSys& R);
  bool spiro_component(const std::vector<int>& atoms, int nrings, int centre, SpiroPart& out);
  bool build_fused(const std::vector<std::vector<int>>& rings, FusedSys& F);
  bool number_fused(const FusedSys& F, std::vector<FusedNumbering>& out);
  std::vector<std::vector<int>> sssr_;
  Component ring_component(const std::vector<int>& cyc) const;
  int order_index(const Cand& c, const std::vector<int>& o) const;
  void hydro(const Cand& c, const std::vector<std::string>& L, std::string& ind, std::vector<std::string>& hyd) const;
  std::optional<Built> assemble(const Spec& sp, const Cand& c, const std::vector<int>& o);
  std::string hydride(const Cand& c, const std::vector<int>& o, bool omit_unsat, bool& saturated) const;
  std::string lambda_name(std::string name, const Cand& c, const std::vector<int>& o) const;
  std::optional<Named> sub(int q, int from, int bo);
  std::optional<Named> sub_impl(int q, int from, int bo);
  std::optional<Named> build_sub(Mode mode, int root, int from, int link = -1, const std::string& word = "carbonyl");
  std::optional<Named> sulfonyl_name(int S, int from);
  std::optional<Named> acyl_name(int c, int from);
  char carbonic_kind(int c, int h) const;
  std::optional<std::pair<std::string, std::vector<Key>>> carbonic_name(int c);
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
  // prefix-only groups on a carbon without a double bond to a heteroatom: nitroso, azido,
  // iso(thio)cyanato, (thio)cyanato, hydroperoxy and alkylperoxy
  pseudo_.assign(n(), PNone);
  pseudoIn_.assign(n(), 0);
  {
    auto plainC = [&](int c) {
      if (at_[c].z != 6) return false;
      const NAtom& C = at_[c];
      for (std::size_t k = 0; k < C.nb.size(); ++k)
        if (at_[C.nb[k]].z != 6 && C.bo[k] != 1 && !C.ab[k]) return false;
      return true;
    };
    auto other = [&](int x, int not_) { return at_[x].nb[0] == not_ ? at_[x].nb[1] : at_[x].nb[0]; };
    for (int a = 0; a < n(); ++a) {
      const NAtom& A = at_[a];
      if (A.ring || A.arom || A.q != 0 || A.h != 0 || deg(a) != 2 || nitro_[a] || nitroO[a]) continue;
      int s1 = -1, d1 = -1, t1 = -1;
      for (std::size_t k = 0; k < A.nb.size(); ++k) {
        if (A.bo[k] == 1) (s1 < 0 ? s1 : t1) = A.nb[k];
        else if (A.bo[k] == 2) d1 = A.nb[k];
      }
      if (A.z == 7 && s1 >= 0 && d1 >= 0) {
        const NAtom& D = at_[d1];
        if (D.ring) continue;
        // C-N=O on a plain carbon, or N-N=O on an acyclic or ring nitrogen
        const bool host_ok = plainC(s1) || (at_[s1].z == 7 && at_[s1].q == 0 && !nitro_[s1]) ||
                             (at_[s1].z == 16 && deg(s1) == 2 && at_[s1].q == 0 && !at_[s1].ring && plainC(other(s1, a)));  // thionitrites
        if (D.z == 8 && deg(d1) == 1 && D.q == 0 && host_ok) {
          pseudo_[a] = Nitroso, pseudoIn_[d1] = 1;
        } else if (D.z == 7 && D.q == 1 && deg(d1) == 2 && plainC(s1)) {
          const int e = other(d1, a);
          if (at_[e].z == 7 && at_[e].q == -1 && deg(e) == 1 && order(d1, e) == 2) pseudo_[a] = Azido, pseudoIn_[d1] = pseudoIn_[e] = 1;
        } else if (D.z == 6 && D.q == 0 && deg(d1) == 2 && D.h == 0 && plainC(s1)) {
          const int e = other(d1, a);
          if ((at_[e].z == 8 || at_[e].z == 16) && deg(e) == 1 && at_[e].q == 0 && order(d1, e) == 2)
            pseudo_[a] = at_[e].z == 8 ? Isocyanato : Isothiocyanato, pseudoIn_[d1] = pseudoIn_[e] = 1;
        }
      } else if ((A.z == 8 || A.z == 16) && s1 >= 0 && t1 >= 0) {
        for (int pass = 0; pass < 2; ++pass) {
          const int c = pass ? t1 : s1, g = pass ? s1 : t1;
          if (!plainC(c)) continue;
          const NAtom& G = at_[g];
          if (G.z == 6 && !G.ring && deg(g) == 2 && G.q == 0) {  // -O-C#N, -S-C#N
            const int e = other(g, a);
            if (at_[e].z == 7 && deg(e) == 1 && order(g, e) == 3 && !pseudoIn_[a]) {
              pseudo_[a] = A.z == 8 ? Cyanato : Thiocyanato, pseudoIn_[g] = pseudoIn_[e] = 1;
              break;
            }
          }
          if (A.z == 8 && G.z == 8 && !G.ring && G.q == 0) {
            if (deg(g) == 1 && G.h == 1) {
              pseudo_[a] = Hydroperoxy, pseudoIn_[g] = 1;
              break;
            }
            if (deg(g) == 2 && plainC(other(g, a))) {
              pseudo_[a] = Peroxy;  // each oxygen of R-O-O-R' is an anchor
              break;
            }
          }
        }
      }
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
      int oxido = 0;
      for (int b : A.nb) {
        if (is_Ominus(b)) {  // N-oxide: an oxido prefix on the cationic nitrogen
          ++oxido;
          continue;
        }
        ok &= at_[b].z == 6 || (A.arom && at_[b].z == 7 && at_[b].arom);
      }
      ok &= oxido <= 1;
      if (ok && !A.ring && std::all_of(A.bo.begin(), A.bo.end(), [](int o) { return o == 1; })) cationN_[a] = 1;
      else if (ok && (A.arom || A.ring)) ringCat_[a] = 1;
    }
    if (is_Ominus(a) && A.bo[0] == 1 && (at_[A.nb[0]].z == 6 || at_[A.nb[0]].z == 16)) anionO_[a] = 1;
  }
  for (int a = 0; a < n(); ++a)
    if (is_Ominus(a) && at_[a].bo[0] == 1 && (cationN_[at_[a].nb[0]] || ringCat_[at_[a].nb[0]])) anionO_[a] = 2;  // N-oxide oxygen
  ionC_.assign(n(), Cls::None);
  ringAn_.assign(n(), 0);
  for (int a = 0; a < n(); ++a) {  // ring anions: indol-1-ide, 9H-fluoren-9-ide, piperidin-1-ide
    const NAtom& A = at_[a];
    if (!A.ring || A.q != -1 || (A.z != 7 && A.z != 6)) continue;
    int rb = 0, other = 0;
    for (std::size_t k = 0; k < A.nb.size(); ++k) {
      if (ringbond_[a][k] && A.bo[k] == 1) ++rb;
      else if (ringbond_[a][k] && A.ab[k]) ++rb;
      else ++other;
    }
    if (rb == 2 && (A.z == 7 ? other == 0 && A.h == 0 : other + A.h == 1)) ringAn_[a] = 1;
  }
  for (int a = 0; a < n(); ++a) {
    const NAtom& A = at_[a];
    if (A.ring || A.arom) continue;
    const bool singles = std::all_of(A.bo.begin(), A.bo.end(), [](int o) { return o == 1; });
    auto all_C = [&] { return std::all_of(A.nb.begin(), A.nb.end(), [&](int b) { return at_[b].z == 6; }); };
    if (A.z == 7 && A.q == 1 && deg(a) == 2 && A.h == 0) {  // C-N+#N
      int c = -1, t = -1;
      for (std::size_t k = 0; k < A.nb.size(); ++k) {
        const int b = A.nb[k];
        if (A.bo[k] == 1 && at_[b].z == 6) c = b;
        if (A.bo[k] == 3 && at_[b].z == 7 && deg(b) == 1 && at_[b].q == 0) t = b;
      }
      if (c >= 0 && t >= 0) ionC_[a] = Cls::Diazonium, pseudoIn_[t] = 1;
    } else if (A.z == 16 && A.q == 1 && deg(a) == 3 && A.h == 0 && singles && all_C()) {
      ionC_[a] = Cls::Sulfonium;
    } else if (A.z == 15 && A.q == 1 && deg(a) == 4 && A.h == 0 && singles && all_C()) {
      ionC_[a] = Cls::Phosphonium;
    } else if (A.z == 5 && A.q == -1 && deg(a) == 4 && A.h == 0 && singles) {  // borates: F, Cl, OH, C ligands
      bool ok = true;
      for (int b : A.nb) ok &= at_[b].z == 6 || ((at_[b].z == 9 || at_[b].z == 17) && deg(b) == 1) || is_OH(b);
      if (ok) {
        ionC_[a] = Cls::Borate;
        for (int b : A.nb)
          if (at_[b].z != 6) pseudoIn_[b] = 1;
      }
    } else if (A.z == 7 && A.q == -1 && deg(a) == 2 && A.h == 0 && singles) {  // sulfonyl azanides R-SO2-N(-)-R'
      int S = -1, C = -1;
      for (int b : A.nb) {
        if (at_[b].z == 16) S = b;
        else if (at_[b].z == 6) C = b;
      }
      int dOs = 0, cS = 0;
      if (S >= 0)
        for (std::size_t k = 0; k < at_[S].nb.size(); ++k) {
          const int y = at_[S].nb[k];
          dOs += at_[S].bo[k] == 2 && at_[y].z == 8 && deg(y) == 1;
          cS += at_[S].bo[k] == 1 && at_[y].z == 6;
        }
      if (S >= 0 && C >= 0 && dOs == 2 && cS == 1 && deg(S) == 4) ionC_[a] = Cls::Azanide;
    } else if (A.z == 7 && A.q == 1 && deg(a) == 3) {  // nitric acid and nitrate
      int dOn = 0, oh = 0, om = 0;
      for (std::size_t k = 0; k < A.nb.size(); ++k) {
        const int b = A.nb[k];
        if (A.bo[k] == 2 && at_[b].z == 8 && deg(b) == 1) ++dOn;
        else if (is_OH(b)) ++oh;
        else if (is_Ominus(b)) ++om;
      }
      if (dOn == 1 && oh + om == 2 && om >= 1) {
        ionC_[a] = Cls::Nitrate;
        for (int b : A.nb) pseudoIn_[b] = 1;
      }
    }
  }
  // hydrazides C(=O)-NR-NR2 / C(=O)-NR-N=CR2 and N-hydroxy amides C(=O)-NR-OH
  hydrN_.assign(n(), 0);
  amidN_.assign(n(), 0);
  hzPartner_.assign(n(), -1);
  azoPartner_.assign(n(), -1);
  lambdaS_.assign(n(), 0);
  thioP_.assign(n(), 0);
  nO_.assign(n(), 0);
  hydrN2_.assign(n(), 0);
  hydroxamO_.assign(n(), 0);
  semC_.assign(n(), 0);
  acylSub_.assign(n(), 0);
  // H2N-C(=O)- or H2N-C(=S)- on nitrogen X, its other nitrogen acyclic with only plain carbons
  auto carbamoylC = [&](int y, int X) {
    if (at_[y].z != 6 || at_[y].ring || (dO(y) < 0 && dS(y) < 0)) return -1;
    const int dx = dO(y) >= 0 ? dO(y) : dS(y);
    int na = -1, het = 0;
    for (int w : at_[y].nb) {
      if (w == dx || at_[w].z == 6) continue;
      ++het;
      if (w != X) na = w;
    }
    if (het != 2 || na < 0 || at_[na].z != 7 || at_[na].ring || at_[na].arom || at_[na].q != 0) return -1;
    for (std::size_t k = 0; k < at_[na].nb.size(); ++k) {
      const int w = at_[na].nb[k];
      if (at_[na].bo[k] != 1) return -1;
      if (w != y && (at_[w].z != 6 || dO(w) >= 0 || dS(w) >= 0)) return -1;
    }
    return na;
  };
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
    const bool alkoxy = at_[other].z == 8 && deg(other) == 2 && at_[other].q == 0 && [&] {
      const int y = at_[other].nb[0] == N ? at_[other].nb[1] : at_[other].nb[0];
      return at_[y].z == 6 && dO(y) < 0 && dS(y) < 0;
    }();
    if (is_OH(other) || alkoxy) {  // N-hydroxy and N-alkoxy amides
      hydroxamO_[other] = 1;
      amideN_[N] = 1;
      continue;
    }
    if (at_[other].z != 7 || at_[other].ring || at_[other].arom || at_[other].q != 0) continue;
    bool semi = false;  // 1-acyl semicarbazide: named from the carboxamide (P-41: amides before hydrazides)
    for (int y : at_[other].nb) semi |= y != N && carbamoylC(y, other) >= 0;
    if (semi) continue;
    bool ok = true;
    for (std::size_t k = 0; k < at_[other].nb.size(); ++k) {
      const int y = at_[other].nb[k];
      if (y == N) continue;
      if (at_[y].z != 6 || (dO(y) >= 0 && !plain_acyl(y))) ok = false;
      if (at_[other].bo[k] == 2)
        for (int w : at_[y].nb) ok &= w == other || at_[w].z == 6;
    }
    if (!ok) continue;
    hydrN_[N] = 1;
    hydrN2_[other] = 1;
  }
  // hydrazines R2N-NR2 and hydrazones R2C=N-NR2 (no acyl groups); hydroxylamines R2N-OR
  auto plain_sp3_N = [&](int x) {
    return at_[x].z == 7 && !at_[x].ring && !at_[x].arom && at_[x].q == 0 && !hydrN_[x] && !hydrN2_[x] && !nitro_[x];
  };
  for (int N = 0; N < n(); ++N) {
    if (!plain_sp3_N(N)) continue;
    for (std::size_t k = 0; k < at_[N].nb.size(); ++k) {
      const int M = at_[N].nb[k];
      if (at_[N].bo[k] != 1 || M < N || !plain_sp3_N(M)) continue;
      bool ok = true;
      int dbl = 0, nsemi = 0, nacyl = 0;
      for (int X : {N, M})
        for (std::size_t j = 0; j < at_[X].nb.size(); ++j) {
          const int y = at_[X].nb[j];
          if (y == N || y == M) continue;
          if (at_[X].bo[j] == 1 && carbamoylC(y, X) >= 0) {
            ++nsemi;
            continue;
          }
          if (at_[X].bo[j] == 1 && plain_acyl(y)) {
            ++nacyl;
            continue;
          }
          if (at_[y].z != 6 || dO(y) >= 0 || dS(y) >= 0) ok = false;
          if (at_[X].bo[j] == 2) {
            ++dbl;
            for (std::size_t i = 0; i < at_[y].nb.size(); ++i)
              if (at_[y].nb[i] != X && (at_[at_[y].nb[i]].z != 6 || at_[y].bo[i] != 1)) ok = false;
          } else if (at_[X].bo[j] != 1) {
            ok = false;
          }
        }
      if (!ok || dbl > 1 || nsemi > 1 || (nacyl > 0 && nsemi == 0)) continue;
      hzPartner_[N] = M;
      hzPartner_[M] = N;
      for (int X : {N, M})
        for (std::size_t j = 0; j < at_[X].nb.size(); ++j) {
          const int y = at_[X].nb[j];
          if (y == N || y == M || at_[X].bo[j] != 1) continue;
          if (carbamoylC(y, X) >= 0) semC_[y] = 1;
          else if (plain_acyl(y)) acylSub_[y] = 1;
        }
    }
  }
  for (int N = 0; N < n(); ++N) {  // azo compounds R-N=N-R'
    const NAtom& A = at_[N];
    if (A.z != 7 || A.ring || A.arom || A.q != 0 || deg(N) > 2) continue;
    for (std::size_t k = 0; k < A.nb.size(); ++k) {
      const int M = A.nb[k];
      if (A.bo[k] != 2 || at_[M].z != 7 || at_[M].ring || at_[M].arom || at_[M].q != 0 || deg(M) > 2) continue;
      bool ok = true;
      for (int X : {N, M})
        for (std::size_t j = 0; j < at_[X].nb.size(); ++j) {
          const int y = at_[X].nb[j];
          if (y == N || y == M) continue;
          if (at_[X].bo[j] != 1 || at_[y].z != 6 || dO(y) >= 0 || dS(y) >= 0) ok = false;
        }
      if (ok) azoPartner_[N] = M;
    }
  }
  for (int N = 0; N < n(); ++N) {
    if (!plain_sp3_N(N) || hzPartner_[N] >= 0 || amideN_[N]) continue;
    for (std::size_t k = 0; k < at_[N].nb.size(); ++k) {
      const int O = at_[N].nb[k];
      if (at_[O].z != 8 || at_[O].ring || at_[O].q != 0 || at_[N].bo[k] != 1) continue;
      bool ok = is_OH(O);
      if (deg(O) == 2) {
        const int y = at_[O].nb[0] == N ? at_[O].nb[1] : at_[O].nb[0];
        ok = at_[y].z == 6 && dO(y) < 0 && dS(y) < 0;
      }
      if (ok) nO_[O] = 1;
    }
  }
  // imines C=N-X marked before the carbons are read, so that the atom X on the far side knows it
  // borders an imine: ring carbons with an exocyclic =N-X (2-(phenylimino)-1,3-thiazolidin-4-one),
  // and acyclic carbons with only carbon or hydrogen neighbours (amidines are read with the carbon)
  for (int a = 0; a < n(); ++a) {
    const NAtom& A = at_[a];
    if (A.z != 6 || (A.arom && !A.ring)) continue;
    int imN = -1, other_double = 0, single_het = 0;
    for (std::size_t k = 0; k < A.nb.size(); ++k) {
      const int b = A.nb[k];
      if (A.bo[k] == 2 && at_[b].z == 7 && !at_[b].ring && !at_[b].arom && at_[b].q == 0 && deg(b) <= 2 && hzPartner_[b] < 0 && !hydrN2_[b])
        imN = b;
      else if (A.bo[k] >= 2) ++other_double;
      else if (at_[b].z != 6 && !ringbond_[a][k]) ++single_het;
    }
    if (imN < 0 || other_double || (!A.ring && single_het)) continue;
    bool ok = true;
    for (int y : at_[imN].nb) {
      if (y == a) continue;
      if (order(imN, y) != 1) ok = false;
      else if (at_[y].z == 6) ok &= dO(y) < 0 || plain_acyl(y);
      else if (at_[y].z == 8) ok &= is_OH(y) || (deg(y) == 2 && [&] {
                                   const int w = at_[y].nb[0] == imN ? at_[y].nb[1] : at_[y].nb[0];
                                   return at_[w].z == 6 && dO(w) < 0;
                                 }());
      else if (at_[y].z == 7) ok &= at_[y].q == 0 && std::all_of(at_[y].bo.begin(), at_[y].bo.end(), [](int o) { return o == 1; }) && !nitro_[y];
      else ok = false;
    }
    if (!ok) continue;
    imineC_[a] = 1;
    imineN_[imN] = 1;
  }
  // sulfur and phosphorus oxo groups
  std::vector<char> okO(n(), 0);
  for (int a = 0; a < n(); ++a) {  // silanes and silyl ethers: four single bonds to C, O or H
    if (at_[a].z != 14) continue;
    int val = at_[a].h;
    for (int o : at_[a].bo) val += o;
    if (val != 4 || at_[a].q != 0 || at_[a].ring) return fail("unsupported silicon at " + sym(a));
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k) {
      const int b = at_[a].nb[k];
      if (at_[a].bo[k] != 1) return fail("unsupported silicon at " + sym(a));
      if (at_[b].z == 8) {
        if (!(is_OH(b) || (deg(b) == 2 && [&] {
              const int y = at_[b].nb[0] == a ? at_[b].nb[1] : at_[b].nb[0];
              return at_[y].z == 6 && dO(y) < 0;
            }())))
          return fail("siloxane or silyl ester at " + sym(a));
        okO[b] = 1;
      } else if (at_[b].z != 6) {
        return fail("unsupported silicon bond at " + sym(a));
      }
    }
  }
  for (int a = 0; a < n(); ++a) {
    const NAtom& A = at_[a];
    if ((A.z != 16 && A.z != 15) || A.arom || A.q != 0) continue;
    if (A.z == 16 && deg(a) == 1 && A.bo[0] == 2 && (at_[A.nb[0]].z == 6 || at_[A.nb[0]].z == 15)) continue;  // C=S / P=S
    std::vector<int> dO_, single, dS_;
    for (std::size_t k = 0; k < A.nb.size(); ++k) {
      const int b = A.nb[k];
      if (A.bo[k] == 2 && at_[b].z == 8 && deg(b) == 1) dO_.push_back(b);
      else if (A.bo[k] == 1) single.push_back(b);
      else if (A.z == 15 && A.bo[k] == 2 && at_[b].z == 16 && deg(b) == 1 && dS_.empty()) dS_.push_back(b);
      else return fail("unsupported multiple bond at " + sym(a));
    }
    if (A.z == 15) {  // phosphanes, phosphine oxides, phosphinic acids, thiophosphates
      if (A.ring) {
        // ring phosphorus with two single ring bonds and one exocyclic substituent (lambda3), or with
        // an exocyclic =O / =S as well (lambda5, "2-oxo-1,3,2lambda5-dioxaphosphinan-2-ol"); in an
        // unsaturated ring system only between two ring O/S atoms, where it can carry no double bond
        int rb = 0, exo = 0;
        bool between_chalcogens = true;
        for (std::size_t k = 0; k < A.nb.size(); ++k) {
          const int b = A.nb[k];
          if (ringbond_[a][k]) {
            ++rb;
            if (A.bo[k] != 1 || A.ab[k]) return fail("cyclic phosphorus compound at " + sym(a));
            between_chalcogens &= at_[b].z == 8 || at_[b].z == 16;
          } else if (A.bo[k] == 1) {
            ++exo;
            if (at_[b].z == 8 && deg(b) == 2 && at_[b].ring) return fail("cyclic phosphorus compound at " + sym(a));
          }
        }
        bool unsat = false;
        {
          std::vector<char> seen(n(), 0);
          std::vector<int> st = {a};
          seen[a] = 1;
          while (!st.empty()) {
            const int x = st.back();
            st.pop_back();
            for (std::size_t k = 0; k < at_[x].nb.size(); ++k) {
              if (!ringbond_[x][k]) continue;
              unsat |= at_[x].bo[k] != 1 || at_[x].ab[k];
              const int y = at_[x].nb[k];
              if (!seen[y]) seen[y] = 1, st.push_back(y);
            }
          }
        }
        if (rb != 2 || exo + A.h != 1 || dO_.size() + dS_.size() > 1 || (unsat && !between_chalcogens) || A.q != 0)
          return fail("cyclic phosphorus compound at " + sym(a));
        if (!dO_.empty() || !dS_.empty()) lambdaS_[a] = 5;
        for (int b : dO_) okO[b] = 1;
        for (int b : dS_) thioP_[b] = 1;
        for (int b : single) {
          if (ringbond_[a][static_cast<std::size_t>(std::find(A.nb.begin(), A.nb.end(), b) - A.nb.begin())]) continue;
          if (is_Ominus(b)) anionO_[b] = 3;
          else if (at_[b].z == 8) okO[b] = 1;
          else if (at_[b].z == 16) thioP_[b] = 1;
        }
        acyl_[a] = Cls::RingP;
        continue;
      }
      int nC = 0, nOH = 0, nOR = 0, nSR = 0, other = 0;
      for (int b : single) {
        if (at_[b].z == 6) ++nC;
        else if (is_OH(b)) ++nOH;
        else if (at_[b].z == 8 && deg(b) == 2 && [&] {
                   const int y = at_[b].nb[0] == a ? at_[b].nb[1] : at_[b].nb[0];
                   return at_[y].z == 6 && dO(y) < 0 && !at_[b].ring;
                 }())
          ++nOR;
        else if (at_[b].z == 16 && deg(b) == 2 && [&] {
                   const int y = at_[b].nb[0] == a ? at_[b].nb[1] : at_[b].nb[0];
                   return at_[y].z == 6 && dO(y) < 0 && !at_[b].ring;
                 }())
          ++nSR;
        else ++other;
      }
      Cls k = Cls::None;
      if (dO_.empty() && dS_.empty() && !other && nOH == 0 && nOR == 0 && nSR == 0 && nC + A.h == 3 && nC >= 1) k = Cls::Phosphane;
      else if (dO_.size() == 1 && dS_.empty() && single.size() == 3 && nC == 3) k = Cls::PhosOxide;
      else if (dO_.size() == 1 && dS_.empty() && single.size() == 3 && nC == 2 && nOH == 1) k = Cls::PhosphinicAcid;
      else if (single.size() == 3 && !other && nC == 0 && nOH == 0 && nOR + nSR == 3 &&
               ((dS_.size() == 1 && dO_.empty() && nSR <= 1) || (dO_.size() == 1 && dS_.empty() && nSR == 1)))
        k = Cls::PhosThio;  // phosphorothioate and phosphorodithioate triesters
      if (k != Cls::None) {
        acyl_[a] = k;
        for (int b : dO_) okO[b] = 1;
        for (int b : single)
          if (at_[b].z == 8) okO[b] = 1;
        for (int b : single)
          if (at_[b].z == 16) thioP_[b] = 1;
        for (int b : dS_) thioP_[b] = 1;
        continue;
      }
      if (!dS_.empty()) return fail("unsupported thiophosphorus group at " + sym(a));
    }
    if (A.z == 16) {
      if (dO_.empty()) continue;
      if (A.h) return fail("hydrido sulfur oxide at " + sym(a));
      if (A.ring) {  // ring sulfoxide / sulfone: oxo prefixes on a lambda4 / lambda6 ring sulfur (P-14.1.3)
        if (dO_.size() > 2 || single.size() != 2) return fail("unsupported cyclic sulfur oxide at " + sym(a));
        for (int b : single)
          if (!ringbond_[a][static_cast<std::size_t>(std::find(A.nb.begin(), A.nb.end(), b) - A.nb.begin())])
            return fail("unsupported cyclic sulfur oxide at " + sym(a));
        lambdaS_[a] = static_cast<char>(dO_.size());
        for (int b : dO_) okO[b] = 1;
        continue;
      }
      int nC = 0, nOH = 0, nOm = 0, nOR = 0, nX = 0;
      int N = -1;
      for (int b : single) {
        if (at_[b].z == 6) ++nC;
        else if (is_OH(b)) ++nOH;
        else if (is_Ominus(b)) ++nOm;
        else if (at_[b].z == 7 && !nitro_[b]) N = b;
        else if (at_[b].z == 8 && deg(b) == 2 && !at_[b].ring && [&] {
                   const int y = at_[b].nb[0] == a ? at_[b].nb[1] : at_[b].nb[0];
                   return at_[y].z == 6 && dO(y) < 0 && dS(y) < 0;
                 }())
          ++nOR;
        else if ((at_[b].z == 9 || at_[b].z == 17 || at_[b].z == 35 || at_[b].z == 53) && deg(b) == 1) ++nX;
        else return fail("sulfur ester or halide at " + sym(a));
      }
      if (dO_.size() == 2 && single.size() == 2 && nC == 0 && nOR == 1 && nOH == 1) {  // sulfuric acid monoester
        acyl_[a] = Cls::SulfateMono;
        for (int b : dO_) okO[b] = 1;
        for (int b : single) okO[b] = 1;
        continue;
      }
      if (dO_.size() == 1 && single.size() == 2 && nC == 1 && nOH == 1) {  // sulfinic acid
        acyl_[a] = Cls::SulfinicAcid;
        for (int b : dO_) okO[b] = 1;
        for (int b : single)
          if (at_[b].z == 8) okO[b] = 1;
        continue;
      }
      if (dO_.size() == 2 && single.size() == 2 && nC == 1 && (nOR == 1 || nX == 1)) {
        acyl_[a] = nOR ? Cls::SulfEster : Cls::SulfonylHalide;
        for (int b : dO_) okO[b] = 1;
        for (int b : single)
          if (at_[b].z == 8) okO[b] = 1;
        continue;
      }
      if (dO_.size() == 2 && single.size() == 2 && nC == 0) {
        // sulfuric acid, hydrogen sulfate and sulfate; sulfate monoester anions; sulfamic acids and sulfamates
        auto plainN = [&](int x) {
          if (at_[x].ring || at_[x].arom || at_[x].q != 0 || nitro_[x] || pseudo_[x]) return false;
          for (std::size_t k = 0; k < at_[x].nb.size(); ++k) {
            const int y = at_[x].nb[k];
            if (y == a) continue;
            if (at_[x].bo[k] != 1 || at_[y].z != 6 || dO(y) >= 0 || dS(y) >= 0) return false;
          }
          return true;
        };
        Cls k = Cls::None;
        if (nOR == 1 && nOm == 1) k = Cls::SulfateAnion;
        else if (N >= 0 && plainN(N) && nOH == 1) k = Cls::SulfamicAcid;
        else if (N >= 0 && plainN(N) && nOm == 1) k = Cls::SulfamateAnion;
        else if (nOH + nOm == 2) k = Cls::SulfuricInorg;
        if (k != Cls::None) {
          acyl_[a] = k;
          if (N >= 0) amideN_[N] = 1;
          for (int b : dO_) okO[b] = 1;
          for (int b : single)
            if (at_[b].z == 8) okO[b] = 1;
          continue;
        }
      }
      if (dO_.size() == 2 && single.size() == 2 && nC == 0) {
        // a ring nitrogen as the attachment atom ("indole-1-sulfonic acid", "piperidine-1-sulfonamide",
        // "1-(piperidine-1-sulfonyl)piperazine"); sulfuric diamides; sulfamate esters
        auto plainN = [&](int x) {
          if (at_[x].z != 7 || at_[x].ring || at_[x].arom || at_[x].q != 0 || nitro_[x] || pseudo_[x]) return false;
          for (std::size_t k = 0; k < at_[x].nb.size(); ++k) {
            const int y = at_[x].nb[k];
            if (y == a) continue;
            if (at_[x].bo[k] != 1 || at_[y].z != 6 || dO(y) >= 0 || dS(y) >= 0) return false;
          }
          return true;
        };
        auto ringN = [&](int x) { return at_[x].z == 7 && at_[x].ring && at_[x].q == 0 && deg(x) == 3; };
        const int x = single[0], y = single[1];
        Cls k = Cls::None;
        int amN = -1, amN2 = -1;
        for (int sw = 0; sw < 2 && k == Cls::None; ++sw) {
          const int r = sw ? y : x, o = sw ? x : y;
          if (!ringN(r)) continue;
          if (is_OH(o)) k = Cls::SulfonicAcid;
          else if (is_Ominus(o)) k = Cls::Sulfonate;
          else if (plainN(o)) k = Cls::Sulfonamide, amN = o;
          else if (ringN(o)) k = Cls::Sulfonyl;
        }
        if (k == Cls::None && plainN(x) && plainN(y)) k = Cls::SulfuricDiamide, amN = x, amN2 = y;
        if (k == Cls::None && nOR == 1 && (plainN(x) || plainN(y))) k = Cls::SulfamateEster, amN = plainN(x) ? x : y;
        if (k != Cls::None) {
          acyl_[a] = k;
          if (amN >= 0) amideN_[amN] = 1;
          if (amN2 >= 0) amideN_[amN2] = 1;
          for (int b : dO_) okO[b] = 1;
          for (int b : single)
            if (at_[b].z == 8) okO[b] = 1;
          continue;
        }
      }
      if (dO_.size() == 1 && nC == 2 && single.size() == 2) acyl_[a] = Cls::Sulfinyl;
      else if (dO_.size() == 2 && single.size() == 2 && nC == 2) acyl_[a] = Cls::Sulfonyl;
      else if (dO_.size() == 2 && single.size() == 2 && nC == 1 && nOH == 1) acyl_[a] = Cls::SulfonicAcid;
      else if (dO_.size() == 2 && single.size() == 2 && nC == 1 && nOm == 1) acyl_[a] = Cls::Sulfonate;
      else if (dO_.size() == 2 && single.size() == 2 && nC == 1 && N >= 0) {
        if (at_[N].ring) acyl_[a] = Cls::Sulfonyl;
        else {
          int hz = -1;  // sulfonohydrazides R-SO2-NH-NR2
          for (int y : at_[N].nb) {
            if (y == a) continue;
            if ((at_[y].z == 17 || at_[y].z == 35) && deg(y) == 1) continue;  // N-chloro, N-bromo
            if (at_[y].z == 7 && !at_[y].ring && !at_[y].arom && at_[y].q == 0 && hz < 0 && [&] {
                  for (std::size_t k = 0; k < at_[y].nb.size(); ++k) {
                    const int w = at_[y].nb[k];
                    if (w == N) continue;
                    if (at_[y].bo[k] != 1 || at_[w].z != 6 || dO(w) >= 0 || dS(w) >= 0) return false;
                  }
                  return true;
                }()) {
              hz = y;
              continue;
            }
            if (at_[y].z != 6) return fail("N-heteroatom sulfonamide at " + sym(N));
            if (dO(y) >= 0 && !plain_acyl(y)) {  // sulfonylureas and sulfonyl carbamates: C(=O) between N and N/O
              int het = 0;
              bool ok2 = !at_[y].ring;
              for (int w : at_[y].nb)
                if (at_[w].z != 6 && w != dO(y)) ++het, ok2 &= at_[w].z == 7 || at_[w].z == 8;
              if (het != 2 || !ok2) return fail("N-acyl sulfonamide at " + sym(N));
            }
          }
          if (hz >= 0) {
            acyl_[a] = Cls::SulfonoHydrazide;
            hydrN_[N] = 1;
            hydrN2_[hz] = 1;
          } else {
            acyl_[a] = Cls::Sulfonamide;
            amideN_[N] = 1;
          }
        }
      } else return fail("unsupported sulfur oxide at " + sym(a));
      for (int b : dO_) okO[b] = 1;
      for (int b : single)
        if (at_[b].z == 8) okO[b] = 1;
    } else {
      int nOH = 0, nC = 0, bridge = -1, nPOP = 0, nOm = 0;
      for (int b : single) {
        if (is_OH(b)) ++nOH;
        else if (is_Ominus(b)) ++nOm;
        else if (at_[b].z == 6) ++nC;
        else if (at_[b].z == 8 && deg(b) == 2) {
          const int y = at_[b].nb[0] == a ? at_[b].nb[1] : at_[b].nb[0];
          if (at_[b].ring) return fail("cyclic phosphate at " + sym(a));
          if (at_[y].z == 15) ++nPOP;  // P-O-P: di- and triphosphates
          else if (at_[y].z != 6 || (dO(y) >= 0 && !plain_acyl(y))) return fail("mixed phosphoric anhydride at " + sym(a));
          else bridge = b;
        } else return fail("unsupported phosphorus group at " + sym(a));
      }
      if (nPOP) {
        if (dO_.size() != 1 || single.size() != 3 || nC || nOm) return fail("unsupported polyphosphate at " + sym(a));
        acyl_[a] = Cls::PolyP;
        okO[dO_[0]] = 1;
        for (int b : single)
          if (at_[b].z == 8) okO[b] = 1;
        continue;
      }
      int nOR = 0;
      for (int b : single)
        if (at_[b].z == 8 && deg(b) == 2) ++nOR;
      if (dO_.size() != 1 || single.size() != 3 || nC > 1) return fail("phosphorus group other than a phosphoric or phosphonic acid derivative at " + sym(a));
      if (nOm) {  // phosphate and phosphonate anions: "ethyl hydrogen phosphate", "methylphosphonate"
        acyl_[a] = Cls::PhosAnion;
        for (int b : single)
          if (is_Ominus(b)) anionO_[b] = 3;
      } else if (nOH == 3) acyl_[a] = Cls::PhosAcidInorg;
      else if (nOH == 2 && nC == 1) acyl_[a] = Cls::Phosphonic;
      else if (nOH == 2 && bridge >= 0) acyl_[a] = Cls::PhosMono;
      else if (nOR >= 1) acyl_[a] = Cls::PhosEster;  // di-/triesters, phosphonate esters
      else return fail("unsupported phosphorus group at " + sym(a));
      okO[dO_[0]] = 1;
      for (int b : single)
        if (at_[b].z == 8) okO[b] = 1;
    }
  }
  // carbons first: their groups decide how the heteroatoms around them are read
  std::vector<int> visit;
  for (int a = 0; a < n(); ++a)
    if (at_[a].z == 6) visit.push_back(a);
  for (int a = 0; a < n(); ++a)
    if (at_[a].z != 6) visit.push_back(a);
  carbX_.assign(n(), 0);
  // carbonic and carbamic acids, their esters, thioesters and anions: C(=O/S)(X)(Y), X and Y each
  // a plain amine nitrogen, OH/SH, O-/S- or OR/SR, with at most one nitrogen
  auto classify_carbonic = [&](int c, const std::vector<int>& hs) {
    if (at_[c].ring || hs.size() != 2) return false;
    int nN = 0, nH = 0, nM = 0;
    for (int h : hs) {
      const char k = carbonic_kind(c, h);
      if (!k) return false;
      nN += k == 'N', nH += k == 'H', nM += k == '-';
    }
    if (nN > 1) return false;
    acyl_[c] = nM ? Cls::CarbonicAnion : nH ? Cls::CarbonicAcid : nN ? Cls::Carbamate : Cls::Carbonate;
    for (int h : hs) {
      if (at_[h].z == 7) amideN_[h] = 1;
      else carbX_[h] = 1;
    }
    carbX_[dO(c) >= 0 ? dO(c) : dS(c)] = 1;
    return true;
  };
  for (int a : visit) {
    const NAtom& A = at_[a];
    int val = A.h;
    for (int o : A.bo) val += o;
    if (pseudoIn_[a]) continue;
    if (pseudo_[a]) continue;
    if (ionC_[a] != Cls::None) continue;
    if (ringAn_[a]) continue;
    if (A.q != 0 && !nitro_[a] && !nitroO[a] && !cationN_[a] && !ringCat_[a] && !anionO_[a] && !(carbX_[a] && A.q == -1))
      return fail("charged atom " + sym(a) + " (ion type not covered)");
    switch (A.z) {
      case 6: {
        if (val != 4) return fail("carbon " + sym(a) + " with unusual valence (radical or carbene)");
        if (semC_[a]) {
          acyl_[a] = dS(a) >= 0 ? Cls::ThioSemicarbazide : Cls::Semicarbazide;
          for (int b : A.nb)
            if (at_[b].z == 7 && hzPartner_[b] < 0) amideN_[b] = 1;
          break;
        }
        if (acylSub_[a]) {
          acyl_[a] = Cls::AcylSub;
          break;
        }
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
              if (y != a && !imineN_[y] && pseudo_[y] != Nitroso && !(at_[y].z == 16 && acyl_[y] == Cls::Sulfonamide) &&
                  (at_[y].z != 6 || (dO(y) >= 0 && !plain_acyl(y)) || dS(y) >= 0))
                return false;
            return true;
          };
          if (A.ring || hs.empty()) acyl_[a] = Cls::Thione;
          else if (hs.size() == 1 && at_[hs[0]].z == 7 && at_[hs[0]].ring) acyl_[a] = Cls::Thione;  // N-thioacyl heterocycle
          else if (hs.size() == 1 && plainN(hs[0])) acyl_[a] = Cls::Thioamide, amideN_[hs[0]] = 1;
          else if (hs.size() == 2 && plainN(hs[0]) && plainN(hs[1])) acyl_[a] = Cls::Thiourea, amideN_[hs[0]] = amideN_[hs[1]] = 1;
          else if (hs.size() == 2 && ((plainN(hs[0]) && at_[hs[1]].z == 7 && at_[hs[1]].ring) || (plainN(hs[1]) && at_[hs[0]].z == 7 && at_[hs[0]].ring)))
            acyl_[a] = Cls::Thioamide, amideN_[at_[hs[0]].ring ? hs[1] : hs[0]] = 1;
          else if (!classify_carbonic(a, hs)) return fail("thiocarbonyl ester, acid or halide at " + sym(a));
          if (acyl_[a] == Cls::Thioamide && A.h > 1) return fail("unsupported thioformamide at " + sym(a));
          break;
        }
        if (imN >= 0 && (hydrN2_[imN] || hzPartner_[imN] >= 0)) {  // hydrazone carbon: an ylidene on its nitrogen
          if (nd || nt) return fail("unsupported carbon " + sym(a));
          break;
        }
        if (imineC_[a]) break;  // read above
        if (imN >= 0) {  // imines and oximes (C=N-H, C=N-R, C=N-OH, C=N-OR); amidines and guanidines
          if (nd || nt) return fail("unsupported carbon " + sym(a));
          // a plain acyclic nitrogen (H or non-acyl carbon neighbours only)
          auto plainNx = [&](int b) {
            if (at_[b].z != 7 || at_[b].ring || at_[b].arom || at_[b].q != 0) return false;
            for (std::size_t k = 0; k < at_[b].nb.size(); ++k) {
              const int y = at_[b].nb[k];
              if (y == a) continue;
              if (at_[b].bo[k] != 1 || at_[y].z != 6 || dO(y) >= 0 || dS(y) >= 0) return false;
            }
            return true;
          };
          // isothioureas R-S-C(=NR)-NR2 (carbamimidothioates) and their free acids (P-66.4.1.2)
          if (!A.ring && plainNx(imN) && deg(imN) <= 2) {
            int nN = -1, sS = -1, other = 0;
            for (int b : A.nb) {
              if (b == imN) continue;
              if (at_[b].z == 7 && order(a, b) == 1 && plainNx(b)) nN = b;
              else if (at_[b].z == 16 && order(a, b) == 1 && !at_[b].ring && at_[b].q == 0 &&
                       (is_SH(b) || (deg(b) == 2 && [&] {
                          const int y = at_[b].nb[0] == a ? at_[b].nb[1] : at_[b].nb[0];
                          return at_[y].z == 6 && dO(y) < 0 && dS(y) < 0;
                        }())))
                sS = b;
              else ++other;
            }
            if (nN >= 0 && sS >= 0 && !other) {
              acyl_[a] = is_SH(sS) ? Cls::IsothioAcid : Cls::Isothiourea;
              amidN_[imN] = 1;
              amideN_[nN] = 1;
              carbX_[sS] = 1;
              break;
            }
          }
          // acylguanidines R-CO-N=C(NR2)2: the amide is senior, the guanidine carbon a
          // "(diaminomethylidene)" substituent on its nitrogen
          if (!A.ring && !at_[imN].ring && deg(imN) == 2) {
            int acylC = -1, nplain = 0, other = 0;
            for (int y : at_[imN].nb)
              if (y != a && plain_acyl(y)) acylC = y;
            for (int b : A.nb) {
              if (b == imN) continue;
              if (order(a, b) == 1 && plainNx(b)) ++nplain;
              else ++other;
            }
            if (acylC >= 0 && nplain == 2 && !other) {
              amidN_[imN] = 1;
              break;
            }
          }
          // ring-N carboximidamides: ring N-C(=NR)-NR2, "piperidine-1-carboximidamide"
          if (!A.ring && plainNx(imN) && deg(imN) <= 2) {
            int ringN = -1, nN = -1, other = 0;
            for (int b : A.nb) {
              if (b == imN) continue;
              if (at_[b].z == 7 && at_[b].ring && !at_[b].arom && at_[b].q == 0 && order(a, b) == 1 && ringN < 0) ringN = b;
              else if (order(a, b) == 1 && plainNx(b)) nN = b;
              else ++other;
            }
            if (ringN >= 0 && nN >= 0 && !other) {
              acyl_[a] = Cls::Amidine;
              amidN_[imN] = 1;
              amideN_[nN] = 1;
              break;
            }
          }
          std::vector<int> ns;
          for (int b : A.nb) {
            if (b == imN || at_[b].z == 6) continue;
            if (at_[b].z != 7 || at_[b].ring || at_[b].arom || at_[b].q != 0 || order(a, b) != 1)
              return fail("imidate, amidoxime or related C(=N)X group at " + sym(a));
            for (int y : at_[b].nb)
              if (y != a && !imineN_[y] && (at_[y].z != 6 || dO(y) >= 0)) return fail("substituted amidine nitrogen at " + sym(b));
            ns.push_back(b);
          }
          if (!ns.empty()) {
            if (A.ring || at_[imN].ring || deg(imN) > 2) return fail("cyclic or substituted amidine at " + sym(a));
            for (int y : at_[imN].nb)
              if (y != a && !((at_[y].z == 6 && dO(y) < 0) || is_OH(y))) return fail("amidine imino nitrogen with a heteroatom or acyl at " + sym(imN));
            acyl_[a] = ns.size() == 1 ? Cls::Amidine : Cls::Guanidine;
            amidN_[imN] = 1;
            for (int b : ns) amideN_[b] = 1;
            for (int y : at_[imN].nb)
              if (is_OH(y)) hydroxamO_[y] = 1;  // amidoxime: N'-hydroxy
            break;
          }
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
                if (at_[y].z == 15 && acyl_[y] != Cls::None && acyl_[y] != Cls::PolyP) break;  // acyl phosphate: oxo + (phosphonooxy)
                if (at_[y].z != 6) return fail("anhydride at " + sym(a));
                if (dO(y) >= 0) {  // carboxylic anhydride R-CO-O-CO-R'
                  int het = 0;
                  for (int w : at_[y].nb) het += at_[w].z != 6 && w != dO(y);
                  if (het != 1 || at_[y].ring) return fail("anhydride at " + sym(a));
                  acyl_[a] = Cls::Anhydride;
                } else {
                  acyl_[a] = Cls::Ester;
                }
              } else return fail("unsupported acyl group at " + sym(a));
            } else if (at_[x].z == 7) {
              if (at_[x].ring) {
                acyl_[a] = Cls::Ketone;  // N-acyl heterocycle: named as a ketone
              } else if (nitro_[x]) {
                return fail("unsupported acyl group at " + sym(a));
              } else {
                for (int y : at_[x].nb) {
                  if (y == a || hydroxamO_[y] || hydrN2_[y] || imineN_[y] || pseudo_[y] == Nitroso) continue;
                  if (at_[y].z == 16 && acyl_[y] == Cls::Sulfonamide) continue;  // N-sulfonyl amide
                  if (at_[y].z != 6) return fail("N-heteroatom amide at " + sym(x));
                  if (dO(y) >= 0) {  // diacylamine: one acyl is the amide, the other an N-acyl prefix
                    int het = 0;
                    for (int w : at_[y].nb) het += at_[w].z != 6 && w != dO(y);
                    bool urea = het == 2 && !at_[y].ring;
                    for (int w : at_[y].nb) urea &= w == dO(y) || at_[w].z == 7;
                    if ((het != 1 && !urea) || at_[y].ring) return fail("imide or diacylamine at " + sym(x));
                  }
                }
                acyl_[a] = hydrN_[x] ? Cls::Hydrazide : Cls::Amide;
                amideN_[x] = 1;
              }
            } else if (at_[x].z == 9 || at_[x].z == 17 || at_[x].z == 35 || at_[x].z == 53) {
              acyl_[a] = Cls::AcylHalide;
            } else if (at_[x].z == 16 && deg(x) == 2 && !at_[x].ring && acyl_[x] == Cls::None) {
              const int y = at_[x].nb[0] == a ? at_[x].nb[1] : at_[x].nb[0];
              if (at_[y].z != 6 || dO(y) >= 0 || dS(y) >= 0) return fail("unsupported thioester at " + sym(a));
              acyl_[a] = Cls::ThioEster;
            } else return fail("unsupported acyl group at " + sym(a));
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
                if (y != a && !imineN_[y] && pseudo_[y] != Nitroso && !(at_[y].z == 16 && acyl_[y] == Cls::Sulfonamide) &&
                  (at_[y].z != 6 || (dO(y) >= 0 && !plain_acyl(y)) || dS(y) >= 0))
                return false;
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
            } else if (!classify_carbonic(a, hs)) return fail("carbonic acid derivative (carbonic or carbamic acid, halide, ...) at " + sym(a));
          } else return fail("unsupported carbonyl group at " + sym(a));
        } else if (nd > 1) return fail("unsupported carbon " + sym(a));
        break;
      }
      case 7: {
        if (nitro_[a] || A.arom || imineN_[a] || cationN_[a] || ringCat_[a] || hydrN_[a] || hydrN2_[a] || amidN_[a] || hzPartner_[a] >= 0 ||
            azoPartner_[a] >= 0)
          break;
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
          const bool sulfonyl_N = at_[b].z == 16 && (acyl_[b] == Cls::Sulfonamide || acyl_[b] == Cls::Sulfonyl || acyl_[b] == Cls::SulfamicAcid ||
                                                    acyl_[b] == Cls::SulfamateAnion || acyl_[b] == Cls::SulfuricDiamide ||
                                                    acyl_[b] == Cls::SulfamateEster);
          const bool n_amino = (at_[b].z == 7 || acyl_[b] == Cls::RingP) && at_[b].ring && A.bo[k] == 1 && !A.ring;  // N-amino heterocycle
          const bool halo_N = amideN_[a] && (at_[b].z == 17 || at_[b].z == 35) && deg(b) == 1;
          if (at_[b].z != 6 && !ringbond_[a][k] && !sulfonyl_N && !nO_[b] && !n_amino && !imineN_[b] && pseudo_[b] != Nitroso && !halo_N) return fail("N-heteroatom bond at " + sym(a) + " (hydrazines, hydroxylamines, ...)");
        }
        if (!A.ring) amineN_[a] = 1;
        break;
      }
      case 8: {
        if (nitroO[a] || okO[a] || anionO_[a] || hydroxamO_[a] || nO_[a] || carbX_[a]) break;
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
      case 14:
        break;  // checked above
      case 16: {
        if (A.arom || acyl_[a] != Cls::None || lambdaS_[a] || thioP_[a] || carbX_[a]) break;
        if (deg(a) == 1 && A.bo[0] == 2 && at_[A.nb[0]].z == 6) break;  // thiocarbonyl
        if (val != 2) return fail("sulfur " + sym(a) + " in a higher oxidation state (sulfoxides, sulfones, sulfonic acids)");
        if (deg(a) == 0) return fail("hydrogen sulfide");
        for (std::size_t k = 0; k < A.nb.size(); ++k)
          if (A.bo[k] != 1 || (at_[A.nb[k]].z != 6 && pseudo_[A.nb[k]] != Nitroso && !(at_[A.nb[k]].z == 16 && deg(A.nb[k]) == 2 && !at_[A.nb[k]].arom &&
                                                         acyl_[A.nb[k]] == Cls::None && !lambdaS_[A.nb[k]] && !thioP_[A.nb[k]])))
            return fail("unsupported sulfur bond at " + sym(a));
        break;
      }
      case 9: case 17: case 35: case 53:
        if (deg(a) != 1 || A.bo[0] != 1 || (at_[A.nb[0]].z != 6 && acyl_[A.nb[0]] != Cls::SulfonylHalide && !(at_[A.nb[0]].z == 7 && amideN_[A.nb[0]])))
          return fail("halogen " + sym(a) + " not bonded to carbon");
        break;
      default:
        return fail("element " + std::string(element_symbol(static_cast<std::uint8_t>(A.z))) + " is outside the naming scope (C, H, N, O, S, halogens)");
    }
  }
  for (int a = 0; a < n(); ++a)
    if (amideN_[a] || nitro_[a] || imineN_[a] || cationN_[a] || hydrN_[a] || hydrN2_[a] || amidN_[a] || hzPartner_[a] >= 0 || azoPartner_[a] >= 0 ||
        at_[a].ring ||
        at_[a].arom)
      amineN_[a] = 0;
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
            R.names.push_back(templates()[t].name);
            R.locs.push_back(templates()[t].locants);
          }
        }
      }
      if (R.orders.empty() && R.nrings == 2 && (spiro_system(R) || von_baeyer(R))) {
        const int id = static_cast<int>(rs_.size());
        for (int a : atoms) ringsys_[a] = id;
        rs_.push_back(std::move(R));
        continue;
      }
      if (R.orders.empty() && R.nrings >= 3 && spiro_fused(R)) {
        const int id = static_cast<int>(rs_.size());
        for (int a : atoms) ringsys_[a] = id;
        rs_.push_back(std::move(R));
        continue;
      }
      bool generic_mancude = false;
      if (R.orders.empty() && R.nrings == 1 && !hantzsch_widman(R, generic_mancude)) return false;
      if (R.orders.empty() && R.nrings == 2 && fusion_bicycle(R)) generic_mancude = true;
      if (R.orders.empty() && R.nrings >= 3 && fusion_polycycle(R)) generic_mancude = true;
      if (R.orders.empty() && R.nrings >= 2 && von_baeyer_poly(R)) {
        const int id = static_cast<int>(rs_.size());
        for (int a : atoms) ringsys_[a] = id;
        rs_.push_back(std::move(R));
        continue;
      }
      if (R.orders.empty())
        return fail(std::string(R.nrings == 1 ? "monocycle" : "ring system") + " of " + std::to_string(atoms.size()) +
                    " atoms without a supported parent name (fused, bridged, spiro, partially hydrogenated or uncommon heterocycle)");
      R.kind = PKind::Template;
      R.mancude = R.templ >= 0 ? templates()[R.templ].mancude : generic_mancude;
      for (int a : atoms)
        for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
          if (ringbond[a][k] && at_[a].nb[k] > a && at_[a].bo[k] == 2) ++R.unsat;
      if (R.mancude) {
        for (int a : atoms) {
          int rdeg = 0;
          bool dbl = false;
          for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
            if (ringbond[a][k]) ++rdeg, dbl |= at_[a].bo[k] == 2;
          const bool inherent = ((at_[a].z == 8 || at_[a].z == 16 || at_[a].z == 15) && rdeg == 2) || (at_[a].z == 7 && rdeg == 3);
          if (inherent) continue;
          R.matchable.push_back(a);
          if (!dbl) R.X.push_back(a);
        }
        R.ih = static_cast<int>(R.matchable.size() % 2);
        if (!R.ih && !perfect_matching(R.matchable)) R.ih = 2;  // two indicated hydrogens: 1H,5H-benzo[ij]quinolizine
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
// one single-bonded neighbour h of a carbonic acid carbon c: 'N' plain amine nitrogen, 'H' OH or SH,
// '-' O- or S-, 'R' OR or SR with a plain carbon R; 0 for anything else
char Namer::carbonic_kind(int c, int h) const {
  const NAtom& H = at_[h];
  if (H.ring || H.arom) return 0;
  if (H.z == 7) {
    if (H.q != 0 || nitro_[h] || pseudo_[h] || hzPartner_[h] >= 0) return 0;
    for (std::size_t k = 0; k < H.nb.size(); ++k) {
      const int y = H.nb[k];
      if (y == c) continue;
      if (H.bo[k] != 1 || at_[y].z != 6 || dO(y) >= 0 || dS(y) >= 0) return 0;
    }
    return 'N';
  }
  if (H.z != 8 && H.z != 16) return 0;
  if (deg(h) == 1 && H.q == 0 && H.h == 1) return 'H';
  if (deg(h) == 1 && H.q == -1 && H.h == 0) return '-';
  if (deg(h) == 2 && H.q == 0) {
    const int y = H.nb[0] == c ? H.nb[1] : H.nb[0];
    if (at_[y].z == 6 && dO(y) < 0 && dS(y) < 0 && pseudoIn_[y] == 0) return 'R';
  }
  return 0;
}

// functional class names of carbonic and carbamic acids and their thio analogues (P-65.2.1.2,
// P-65.2.1.5): "N-phenylcarbamic acid", "benzyl hydrogen carbonate", "methyl dimethylcarbamodithioate",
// "S-methyl N,N-dimethylcarbamothioate", "O-ethyl S-hydrogen carbonodithioate", "N,N-dimethylcarbamodithioate"
std::optional<std::pair<std::string, std::vector<Key>>> Namer::carbonic_name(int c) {
  const int X = dO(c) >= 0 ? dO(c) : dS(c);
  int N = -1, nS = at_[X].z == 16, nch = 1;
  std::vector<int> chal;
  for (int h : at_[c].nb) {
    if (h == X || at_[h].z == 6) continue;
    if (at_[h].z == 7) N = h;
    else chal.push_back(h), ++nch, nS += at_[h].z == 16;
  }
  const bool letters = nS > 0 && nS < nch;
  static const char* thio[] = {"", "othio", "odithio", "otrithio"};
  const std::string base = std::string(N >= 0 ? "carbam" : "carbon") + thio[nS];
  std::vector<Key> extra;
  std::string npart;
  if (N >= 0) {
    std::vector<std::pair<std::string, Named>> e;
    for (int x : at_[N].nb) {
      if (x == c) continue;
      auto R = sub(x, N, 1);
      if (!R) return std::nullopt;
      extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
      e.push_back({"N", *R});
    }
    npart = format_prefixes(e, false);
  }
  std::map<std::string, std::vector<std::string>> by;  // organyl -> locant letters
  std::vector<std::string> hyd;
  bool anion = false;
  for (int h : chal) {
    const std::string L = at_[h].z == 8 ? "O" : "S";
    if (at_[h].q == -1) {
      anion = true;
    } else if (deg(h) == 1) {
      hyd.push_back(L);
    } else {
      const int r = at_[h].nb[0] == c ? at_[h].nb[1] : at_[h].nb[0];
      auto R = sub(r, h, 1);
      if (!R) return std::nullopt;
      extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
      by[R->text[0] == '(' ? enclose(R->text) : R->text].push_back(L);
    }
  }
  if (by.empty() && !anion) {  // the free acid
    if (letters) {
      fail("thio acid with O/S isomers at " + sym(c));
      return std::nullopt;
    }
    return std::pair{npart + base + "ic acid", extra};
  }
  std::vector<std::string> parts;
  std::vector<std::pair<std::string, std::vector<std::string>>> orgs(by.begin(), by.end());
  std::sort(orgs.begin(), orgs.end(), [](const auto& a, const auto& b) { return alpha_key(a.first) < alpha_key(b.first); });
  for (auto& [org, ls] : orgs) {
    std::sort(ls.begin(), ls.end());
    std::string loc;
    for (const auto& l : ls) loc += (loc.empty() ? "" : ",") + l;
    std::string t = org;
    if (ls.size() == 2) t = (org[0] == '(' || org[0] == '[' || has_digit(org)) ? "bis" + enclose(org) : "di" + org;
    parts.push_back(letters ? loc + "-" + t : t);
  }
  if (!hyd.empty()) {
    std::string loc;
    for (const auto& l : hyd) loc += (loc.empty() ? "" : ",") + l;
    parts.push_back((letters ? loc + "-" : "") + (hyd.size() == 2 ? "dihydrogen" : "hydrogen"));
  }
  std::string out;
  for (const auto& p : parts) out += p + " ";
  return std::pair{out + npart + base + "ate", extra};
}

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
  if (is_OH(x) && acyl_[c] == Cls::CarbonicAcid) return Named{"carboxy", false, {}};  // -NH-COOH: (carboxyamino)
  if (at_[x].z == 16 && at_[x].q == 0 && (acyl_[c] == Cls::CarbonicAcid || acyl_[c] == Cls::Carbonate || acyl_[c] == Cls::Carbamate)) {
    auto R = sub(x, c, 1);  // (methylsulfanyl)carbonyl, sulfanylcarbonyl
    if (!R) return std::nullopt;
    return Named{(R->compound ? enclose(R->text) : R->text) + "carbonyl", true, R->stereo};
  }
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
  if (at_[x].z == 7 && hzPartner_[x] >= 0 && acyl_[c] == Cls::Semicarbazide) {  // hydrazinecarbonyl
    auto R = sub(x, c, 1);
    if (!R) return std::nullopt;
    if (R->text == "hydrazinyl") return Named{"hydrazinecarbonyl", true, R->stereo};
    if (ends_with(R->text, "hydrazin-1-yl")) return Named{R->text.substr(0, R->text.size() - 13) + "hydrazine-1-carbonyl", true, R->stereo};
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

// "dimethoxyphosphoryl", "[hydroxy(methyl)phosphoryl]": phosphorus P reached from `from`.
std::optional<Named> Namer::phosphoryl_name(int P, int from) {
  std::vector<int> xs;
  bool thio = false;
  for (std::size_t k = 0; k < at_[P].nb.size(); ++k) {
    if (at_[P].nb[k] != from && at_[P].bo[k] == 1) xs.push_back(at_[P].nb[k]);
    if (at_[P].bo[k] == 2 && at_[at_[P].nb[k]].z == 16) thio = true;
  }
  return ligands(xs, P, thio ? "phosphorothioyl" : "phosphoryl");
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
    if (A.z == 7 && hzPartner_[q] >= 0) {  // =N-NR2: (2,2-dimethylhydrazinylidene)
      const int M = hzPartner_[q];
      std::vector<std::pair<std::string, Named>> e;
      Named out;
      for (std::size_t k = 0; k < at_[M].nb.size(); ++k) {
        if (at_[M].nb[k] == q) continue;
        auto R = sub(at_[M].nb[k], M, at_[M].bo[k]);
        if (!R) return std::nullopt;
        out.stereo.insert(out.stereo.end(), R->stereo.begin(), R->stereo.end());
        e.push_back({"2", *R});
      }
      out.text = format_prefixes(e, false) + "hydrazinylidene";
      out.compound = !e.empty();
      return out;
    }
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
  switch (pseudo_[q]) {
    case Nitroso: return Named{"nitroso", false, {}};
    case Azido: return Named{"azido", false, {}};
    case Isocyanato: return Named{"isocyanato", false, {}};
    case Isothiocyanato: return Named{"isothiocyanato", false, {}};
    case Cyanato: return Named{"cyanato", false, {}};
    case Thiocyanato: return Named{"thiocyanato", false, {}};
    case Hydroperoxy: return Named{"hydroperoxy", false, {}};
    case Peroxy: {  // -O-O-R: (methylperoxy)
      int g = -1;
      for (int y : A.nb)
        if (y != from) g = y;
      if (g < 0 || at_[g].z != 8 || from < 0 || at_[from].z != 6) return fail_n("unsupported peroxide at " + sym(q));
      const int r = at_[g].nb[0] == q ? at_[g].nb[1] : at_[g].nb[0];
      auto R = sub(r, g, 1);
      if (!R) return std::nullopt;
      return Named{(R->compound ? enclose(R->text) : R->text) + "peroxy", true, R->stereo};
    }
    default: break;
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
      if (others.size() == 1 && at_[others[0]].z == 15 && acyl_[others[0]] == Cls::PhosMono) return Named{"phosphonooxy", true, {}};
      if (others.size() == 1 && at_[others[0]].z == 16 && acyl_[others[0]] == Cls::SulfateMono) return Named{"sulfooxy", true, {}};
      if (others.size() == 1 && at_[others[0]].z == 16 && acyl_[others[0]] == Cls::SulfateAnion) return Named{"sulfonatooxy", true, {}};
      if (others.size() == 1 && at_[others[0]].z == 14) {  // [(trimethylsilyl)oxy]
        auto R = sub(others[0], q, 1);
        if (!R) return std::nullopt;
        return Named{(R->compound ? enclose(R->text) : R->text) + "oxy", true, R->stereo};
      }
      if (others.size() == 1 && at_[others[0]].z == 16 && acyl_[others[0]] == Cls::SulfamateEster) {  // (sulfamoyloxy)
        auto R = sub(others[0], q, 1);
        if (!R) return std::nullopt;
        return Named{(R->compound ? enclose(R->text) : R->text) + "oxy", true, R->stereo};
      }
      if (others.size() == 1 && at_[others[0]].z == 16 && acyl_[others[0]] == Cls::SulfEster) {  // (methanesulfonyl)oxy
        auto R = sulfonyl_name(others[0], q);
        if (!R) return std::nullopt;
        return Named{enclose(R->text) + "oxy", true, R->stereo};
      }
      if (others.size() == 1 && at_[others[0]].z == 15 && acyl_[others[0]] == Cls::PolyP) {
        const int Pp = others[0];
        int oh = 0, singles = 0;
        for (std::size_t k = 0; k < at_[Pp].nb.size(); ++k)
          if (at_[Pp].bo[k] == 1 && at_[Pp].nb[k] != q) ++singles, oh += is_OH(at_[Pp].nb[k]);
        if (oh == singles) return Named{"phosphonooxy", true, {}};
        auto R = phosphoryl_name(Pp, q);
        if (!R) return std::nullopt;
        return Named{enclose(R->text) + "oxy", true, R->stereo};
      }
      if (others.size() == 1 && at_[others[0]].z == 15 &&
          (acyl_[others[0]] == Cls::PhosEster || acyl_[others[0]] == Cls::PhosThio || acyl_[others[0]] == Cls::PhosphinicAcid ||
           acyl_[others[0]] == Cls::PhosAnion ||
           acyl_[others[0]] == Cls::PhosOxide)) {  // (dimethoxyphosphoryl)oxy
        auto R = phosphoryl_name(others[0], q);
        if (!R) return std::nullopt;
        return Named{enclose(R->text) + "oxy", true, R->stereo};
      }
      if (others.size() == 1 && nO_[q] && at_[others[0]].z == 7) {  // -O-NR2: (aminooxy)
        auto R = sub(others[0], q, 1);
        if (!R) return std::nullopt;
        return Named{(R->compound ? enclose(R->text) : R->text) + "oxy", true, R->stereo};
      }
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
    case 14:
      if (others.empty()) return Named{"silyl", false, {}};
      return ligands(others, q, "silyl");
    case 15:
      if (ionC_[q] == Cls::Phosphonium) return ligands(others, q, "phosphaniumyl");
      if (acyl_[q] == Cls::Phosphonic) return Named{"phosphono", false, {}};
      if (acyl_[q] == Cls::PhosEster || acyl_[q] == Cls::PhosAnion || acyl_[q] == Cls::PhosMono || acyl_[q] == Cls::PolyP || acyl_[q] == Cls::PhosOxide ||
          acyl_[q] == Cls::PhosphinicAcid || acyl_[q] == Cls::PhosThio)
        return phosphoryl_name(q, from);
      if (acyl_[q] == Cls::Phosphane) {
        if (others.empty()) return Named{"phosphanyl", false, {}};
        return ligands(others, q, "phosphanyl");
      }
      return fail_n("unsupported phosphorus substituent at " + sym(q));
    case 16: {
      if (ionC_[q] == Cls::Sulfonium) return ligands(others, q, "sulfaniumyl");
      if (thioP_[q] && others.size() == 1 && at_[others[0]].z == 15) {  // [(dimethoxyphosphorothioyl)sulfanyl]
        auto R = phosphoryl_name(others[0], q);
        if (!R) return std::nullopt;
        return Named{enclose(R->text) + "sulfanyl", true, R->stereo};
      }
      switch (acyl_[q]) {
        case Cls::SulfonicAcid: return Named{"sulfo", false, {}};
        case Cls::SulfamicAcid: return Named{"sulfo", false, {}};
        case Cls::SulfinicAcid: return Named{"sulfino", false, {}};
        case Cls::SulfEster: {  // (methoxysulfonyl)
          auto alk = sub(ester_O(q), q, 1);
          if (!alk) return std::nullopt;
          return Named{(alk->compound ? enclose(alk->text) : alk->text) + "sulfonyl", true, alk->stereo};
        }
        case Cls::SulfonylHalide:
          for (int x : A.nb)
            if (at_[x].z == 9 || at_[x].z == 17 || at_[x].z == 35 || at_[x].z == 53)
              return Named{std::string(at_[x].z == 9 ? "fluoro" : at_[x].z == 17 ? "chloro" : at_[x].z == 35 ? "bromo" : "iodo") + "sulfonyl", true, {}};
          return fail_n("unsupported sulfonyl halide at " + sym(q));
        case Cls::Sulfonate: return Named{"sulfonato", false, {}};
        case Cls::SulfamateEster: {
          int N = -1, O = -1;
          for (int x : A.nb) {
            if (at_[x].z == 7) N = x;
            if (at_[x].z == 8 && deg(x) == 2) O = x;
          }
          if (from == O) {  // sulfamoyl, then "oxy" on the oxygen
            std::vector<int> xs;
            for (int x : at_[N].nb)
              if (x != q) xs.push_back(x);
            if (xs.empty()) return Named{"sulfamoyl", false, {}};
            return ligands(xs, N, "sulfamoyl");
          }
          auto alk = sub(O, q, 1);  // (methoxysulfonyl), reached from the nitrogen
          if (!alk) return std::nullopt;
          return Named{(alk->compound ? enclose(alk->text) : alk->text) + "sulfonyl", true, alk->stereo};
        }
        case Cls::SulfuricDiamide:
        case Cls::Sulfonamide: {
          int N = -1;
          for (int x : A.nb)
            if (at_[x].z == 7 && amideN_[x] && (N < 0 || x != from)) N = x;
          if (N < 0) return fail_n("unsupported sulfonamide substituent at " + sym(q));
          if (N == from && acyl_[q] == Cls::SulfuricDiamide) {
            for (int x : A.nb)
              if (at_[x].z == 7 && x != from) N = x;
          } else if (N == from) return sulfonyl_name(q, from);  // N-sulfonyl amides and ureas
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
      if (others.size() == 1 && at_[others[0]].z == 16) {  // -S-S-R: (R-disulfanyl)
        const int S2 = others[0];
        int r = -1;
        for (int y : at_[S2].nb)
          if (y != q) r = y;
        if (r < 0) return Named{"disulfanyl", false, {}};
        if (at_[r].z != 6) return fail_n("unsupported polysulfide at " + sym(q));
        auto R = sub(r, S2, 1);
        if (!R) return std::nullopt;
        return Named{(R->compound ? enclose(R->text) : R->text) + "disulfanyl", true, R->stereo};
      }
      if (others.size() == 1 && pseudo_[others[0]] == Nitroso) return Named{"nitrososulfanyl", true, {}};
      if (others.size() != 1 || at_[others[0]].z != 6) return fail_n("unsupported sulfur substituent at " + sym(q));
      auto R = sub(others[0], q, 1);
      if (!R) return std::nullopt;
      return Named{(R->compound ? enclose(R->text) : R->text) + "sulfanyl", true, R->stereo};
    }
    case 7: {
      if (nitro_[q]) return Named{"nitro", false, {}};
      if (ionC_[q] == Cls::Diazonium) return Named{"diazonio", false, {}};
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
      if (azoPartner_[q] >= 0) {  // -N=N-R: (phenyldiazenyl), with its E/Z
        const int M = azoPartner_[q];
        Named out;
        out.compound = true;
        std::string t = "diazenyl";
        for (int y : at_[M].nb) {
          if (y == q) continue;
          auto R = sub(y, M, 1);
          if (!R) return std::nullopt;
          out.stereo = R->stereo;
          t = (R->compound ? enclose(R->text) : R->text) + "diazenyl";
        }
        if (auto it = stereo_.find(bond_key(q, M)); it != stereo_.end()) {
          t = "(" + std::string(to_string(it->second)) + ")-" + t;
          out.stereo.push_back(it->first);
        }
        out.text = t;
        return out;
      }
      if (hzPartner_[q] >= 0) {  // -NR-NR2: (2-phenylhydrazin-1-yl)
        const int M = hzPartner_[q];
        std::vector<std::pair<std::string, Named>> e;
        Named out;
        for (int X : {q, M})
          for (std::size_t k = 0; k < at_[X].nb.size(); ++k) {
            const int y = at_[X].nb[k];
            if (y == from || y == q || y == M) continue;
            auto R = sub(y, X, at_[X].bo[k]);
            if (!R) return std::nullopt;
            out.stereo.insert(out.stereo.end(), R->stereo.begin(), R->stereo.end());
            e.push_back({X == q ? "1" : "2", *R});
          }
        if (e.empty()) return Named{"hydrazinyl", false, {}};
        out.text = format_prefixes(e, false) + "hydrazin-1-yl";
        out.compound = true;
        return out;
      }
      if (others.size() == 1 && acyl_[others[0]] == Cls::Guanidine) {
        bool plain = true;
        for (int x : at_[others[0]].nb)
          if (x != q && deg(x) != 1) plain = false;
        if (plain) return Named{"carbamimidamido", false, {}};
      }
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
        else if (ends_with(t, "carbonyl") && acyl_[others[0]] != Cls::CarbonicAcid && acyl_[others[0]] != Cls::Carbonate &&
                 acyl_[others[0]] != Cls::Carbamate)
          t = t.substr(0, t.size() - 8) + "carboxamido";
        else if (t == "carboxy") return Named{"carboxyamino", true, ac->stereo};
        else return Named{enclose(t) + "amino", true, ac->stereo};
        return Named{t, ac->compound, ac->stereo};
      }
      return ligands(others, q, "amino");
    }
    case 6: {
      if (acyl_[q] == Cls::Nitrile) return Named{"cyano", false, {}};
      if (acyl_[q] == Cls::Amidine || acyl_[q] == Cls::Guanidine) {  // carbamimidoyl, -C(=NH)NH2 only
        int plainN = 0;
        for (int x : A.nb) {
          if (x == from) continue;
          if (at_[x].z != 7 || deg(x) != 1) return fail_n("substituted carbamimidoyl group at " + sym(q));
          ++plainN;
        }
        if (plainN != 2) return fail_n("unsupported amidine substituent at " + sym(q));
        return Named{"carbamimidoyl", false, {}};
      }
      if ((acyl_[q] == Cls::AcylSub || acyl_[q] == Cls::Semicarbazide || acyl_[q] == Cls::ThioSemicarbazide) && at_[from].z == 6)
        return fail_n("hydrazine acyl group reached from carbon at " + sym(q));
      if (dS(q) >= 0 && (acyl_[q] == Cls::Carbonate || acyl_[q] == Cls::CarbonicAcid || acyl_[q] == Cls::Carbamate)) {
        // (methylsulfanyl)carbonothioyl, ethoxycarbonothioyl, sulfanylcarbonothioyl
        int h = -1;
        for (int x : A.nb)
          if (x != from && x != dS(q) && at_[x].z != 7) h = x;
        if (h >= 0) {
          auto R = sub(h, q, 1);
          if (!R) return std::nullopt;
          return Named{(R->compound ? enclose(R->text) : R->text) + "carbonothioyl", true, R->stereo};
        }
      }
      if (acyl_[q] == Cls::Thioamide || acyl_[q] == Cls::Thiourea || acyl_[q] == Cls::ThioSemicarbazide ||
          (dS(q) >= 0 && acyl_[q] == Cls::Carbamate)) {  // carbamothioyl, from the carbon side only
        int N = -1;
        for (int x : A.nb)
          if (x != from && at_[x].z == 7 && amideN_[x]) N = x;
        if (N < 0) return fail_n("thioacyl substituent at " + sym(q));
        std::vector<int> xs;
        for (int x : at_[N].nb)
          if (x != q) xs.push_back(x);
        if (xs.empty()) return Named{"carbamothioyl", false, {}};
        return ligands(xs, N, "carbamothioyl");
      }
      if (dO(q) >= 0 && at_[from].z != 6 && acyl_[q] != Cls::Carbamate && acyl_[q] != Cls::Urea) return acyl_name(q, from);
      if (dO(q) >= 0 && acyl_[q] == Cls::AcylHalide) {
        for (int x : A.nb)
          if (x != from && (at_[x].z == 9 || at_[x].z == 17 || at_[x].z == 35 || at_[x].z == 53))
            return Named{std::string("carbono") + (at_[x].z == 9 ? "fluor" : at_[x].z == 17 ? "chlor" : at_[x].z == 35 ? "brom" : "iod") + "idoyl", false, {}};
      }
      if (dO(q) >= 0 && acyl_[q] == Cls::ThioEster) {
        int S = -1;
        for (int x : A.nb)
          if (at_[x].z == 16 && x != from) S = x;
        auto R = sub(S, q, 1);
        if (!R) return std::nullopt;
        return Named{(R->compound ? enclose(R->text) : R->text) + "carbonyl", true, R->stereo};
      }
      if (dO(q) >= 0) {
        switch (acyl_[q]) {
          case Cls::Acid: return Named{"carboxy", false, {}};
          case Cls::Anhydride: {  // [(acetyloxy)carbonyl]
            auto alk = sub(ester_O(q), q, 1);
            if (!alk) return std::nullopt;
            return Named{enclose(alk->text) + "carbonyl", true, alk->stereo};
          }
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
  bool nprefix = false;
  auto consume_acyl = [&](int A, int at) -> bool {
    for (int b : at_[A].nb) {
      if (c.in[b] || b == sp.from) continue;
      if (at_[b].z == 6) return fail("acyl carbon " + sym(A) + " with a second carbon");
      if (at_[b].ring) return fail("acyl group on a ring atom outside the parent at " + sym(A));
      used[b] = 1;
      if ((acyl_[A] == Cls::Ester || acyl_[A] == Cls::SulfEster || acyl_[A] == Cls::ThioEster) && b == ester_O(A)) c.ester_o.push_back({b, at});
      if ((acyl_[A] == Cls::Hydrazide || acyl_[A] == Cls::SulfonoHydrazide) && at_[b].z == 7)
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
      if (acyl_[A] == Cls::Amidine && at_[b].z == 7) {
        if (amidN_[b]) c.exo_double.push_back({at, b});
        for (int x : at_[b].nb) {
          if (x == A) continue;
          auto s = sub(x, b, order(x, b));
          if (!s) return false;
          c.pre.push_back({amidN_[b] ? -2 : -1, *s});
          nprefix = true;
        }
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
    if (mol && sp.principal == Cls::RingAnion && ringAn_[p]) c.suffix_at.push_back(p);
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
            (sp.principal == Cls::Peroxol && pseudo_[q] == Hydroperoxy) ||
            (sp.principal == Cls::Olate && is_Ominus(q) && at_[p].z == 6 && acyl_[p] != Cls::Carboxylate)) {
          used[q] = 1;
          c.suffix_at.push_back(p);
          continue;
        }
        if (sp.principal == Cls::Diazonium && ionC_[q] == Cls::Diazonium) {
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
  if (nprefix && (sp.principal == Cls::Amine || sp.principal == Cls::Ammonium || sp.principal == Cls::Amide || sp.principal == Cls::Hydrazide ||
                  sp.principal == Cls::Amidine || sp.principal == Cls::Thioamide || sp.principal == Cls::Sulfonamide ||
                  sp.principal == Cls::SulfonoHydrazide ||
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

int Namer::order_index(const Cand& c, const std::vector<int>& o) const {
  const auto& os = rs_[c.rsys].orders;
  for (std::size_t i = 0; i < os.size(); ++i)
    if (os[i] == o) return static_cast<int>(i);
  return 0;
}

std::vector<std::string> Namer::locants(const Cand& c, const std::vector<int>& o) const {
  std::vector<std::string> L(n());
  if (c.kind == PKind::Template) {
    const auto& lc = rs_[c.rsys].locs[static_cast<std::size_t>(order_index(c, o))];
    for (std::size_t i = 0; i < o.size(); ++i) L[o[i]] = lc[i];
  } else {
    for (std::size_t i = 0; i < o.size(); ++i) L[o[i]] = std::to_string(i + 1);
  }
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

// Monospiro systems of two monocycles (P-24.2.1): spiro[4.5]decane, 1-oxa-8-azaspiro[4.5]decane.
bool Namer::spiro_system(RingSys& R) {
  int centre = -1;
  auto rdeg = [&](int a) {
    int d = 0;
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k) d += ringbond_[a][k];
    return d;
  };
  for (int a : R.atoms) {
    const int d = rdeg(a);
    if (d == 4 && centre < 0) centre = a;
    else if (d != 2) return false;
  }
  if (centre < 0) return false;
  for (int a : R.atoms)
    if (at_[a].arom) return fail("spiro system with an aromatic component at " + sym(a));
  // the two rings, each walked from the centre
  std::vector<std::vector<int>> rings;
  std::set<int> seen;
  for (std::size_t k = 0; k < at_[centre].nb.size(); ++k) {
    const int first = at_[centre].nb[k];
    if (!ringbond_[centre][k] || seen.count(first)) continue;
    std::vector<int> path;
    int prev = centre, cur = first;
    while (cur != centre) {
      path.push_back(cur);
      seen.insert(cur);
      int next = -1;
      for (std::size_t j = 0; j < at_[cur].nb.size(); ++j)
        if (ringbond_[cur][j] && at_[cur].nb[j] != prev) next = at_[cur].nb[j];
      if (next < 0 || path.size() > R.atoms.size()) return false;
      prev = cur;
      cur = next;
    }
    rings.push_back(path);
  }
  if (rings.size() != 2) return false;
  if (rings[0].size() > rings[1].size()) std::swap(rings[0], rings[1]);
  R.kind = PKind::VonBaeyer;
  R.spiro = true;
  R.vb[0] = static_cast<int>(rings[0].size());
  R.vb[1] = static_cast<int>(rings[1].size());
  for (int swap : {0, 1}) {
    if (swap && rings[0].size() != rings[1].size()) break;
    const auto& small = rings[swap];
    const auto& large = rings[1 - swap];
    for (int d1 : {0, 1})
      for (int d2 : {0, 1}) {
        std::vector<int> o;
        if (d1) o.insert(o.end(), small.rbegin(), small.rend());
        else o.insert(o.end(), small.begin(), small.end());
        o.push_back(centre);
        if (d2) o.insert(o.end(), large.rbegin(), large.rend());
        else o.insert(o.end(), large.begin(), large.end());
        R.orders.push_back(o);
      }
  }
  for (int x : R.atoms)
    for (std::size_t k = 0; k < at_[x].nb.size(); ++k)
      if (ringbond_[x][k] && at_[x].nb[k] > x && at_[x].bo[k] > 1) ++R.unsat;
  return true;
}

// Hantzsch-Widman names of heteromonocycles with 3-10 members (P-22.2.2) that have no retained name.
bool Namer::hantzsch_widman(RingSys& R, bool& mancude) {
  const int m = static_cast<int>(R.atoms.size());
  if (m < 3 || m > 10) return true;  // nothing to add; the caller reports the ring
  std::vector<int> cyc = {R.atoms[0]};
  {
    std::set<int> in(R.atoms.begin(), R.atoms.end());
    int prev = -1, cur = R.atoms[0];
    while (static_cast<int>(cyc.size()) <= m) {
      int next = -1;
      for (std::size_t k = 0; k < at_[cur].nb.size(); ++k)
        if (ringbond_[cur][k] && at_[cur].nb[k] != prev) {
          next = at_[cur].nb[k];
          break;
        }
      if (next < 0 || next == R.atoms[0]) break;
      cyc.push_back(next);
      prev = cur;
      cur = next;
    }
    if (static_cast<int>(cyc.size()) != m) return true;
  }
  bool saturated = true, hasN = false, hasP = false;
  for (int i = 0; i < m; ++i) {
    const int a = cyc[i];
    if (at_[a].z != 6 && at_[a].z != 7 && at_[a].z != 8 && at_[a].z != 16 && at_[a].z != 15) return true;
    hasN |= at_[a].z == 7;
    hasP |= at_[a].z == 15;
    if (order(a, cyc[(i + 1) % m]) > 1) saturated = false;
  }
  int senior_z = 15;  // O > S > N > P (P-22.2.2.1.3)
  for (int a : cyc) {
    if (at_[a].z == 8) senior_z = 8;
    if (at_[a].z == 16 && senior_z != 8) senior_z = 16;
    if (at_[a].z == 7 && senior_z == 15) senior_z = 7;
  }
  static const char* sat_n[] = {"", "", "", "iridine", "etidine", "olidine", "inane", "epane", "ocane", "onane", "ecane"};
  static const char* sat_o[] = {"", "", "", "irane", "etane", "olane", "ane", "epane", "ocane", "onane", "ecane"};
  static const char* uns_n[] = {"", "", "", "irine", "ete", "ole", "ine", "epine", "ocine", "onine", "ecine"};
  static const char* uns_o[] = {"", "", "", "irene", "ete", "ole", "ine", "epine", "ocine", "onine", "ecine"};
  std::string stem = saturated ? (hasN ? sat_n[m] : sat_o[m]) : (hasN ? uns_n[m] : uns_o[m]);
  if (m == 6 && hasP) stem = saturated ? "inane" : "inine";  // six-membered ring, last-cited P (6C)
  for (int st = 0; st < m; ++st)
    for (int dir : {1, -1}) {
      std::vector<int> o;
      for (int i = 0; i < m; ++i) o.push_back(cyc[((st + dir * i) % m + m) % m]);
      if (at_[o[0]].z != senior_z) continue;  // locant 1 to the senior heteroatom (P-22.2.2.1.3)
      std::vector<std::string> L;
      for (int i = 0; i < m; ++i) L.push_back(std::to_string(i + 1));
      std::string locs, pre;
      int nhet = 0;
      for (const auto& [z, word] : std::vector<std::pair<int, std::string>>{{8, "oxa"}, {16, "thia"}, {7, "aza"}, {15, "phospha"}}) {
        std::vector<std::string> l;
        for (int i = 0; i < m; ++i)
          if (at_[o[i]].z == z) l.push_back(L[i]);
        if (l.empty()) continue;
        nhet += static_cast<int>(l.size());
        std::sort(l.begin(), l.end(), [](const std::string& x, const std::string& y) { return locant_key(x) < locant_key(y); });
        for (const auto& x : l) locs += (locs.empty() ? "" : ",") + x;
        std::string w = multiplier(static_cast<int>(l.size()), false) + word;
        if (!pre.empty() && pre.back() == 'a' && starts_with_vowel(w)) pre.pop_back();
        pre += w;
      }
      if (pre.back() == 'a' && starts_with_vowel(stem)) pre.pop_back();
      R.orders.push_back(o);
      R.locs.push_back(L);
      R.names.push_back((nhet > 1 ? locs + "-" : "") + pre + stem);
    }
  R.het_numbering = true;
  mancude = !saturated;
  return true;
}

// ---- fusion nomenclature for ortho-fused bicycles (P-25.3, FR-3 to FR-5) -----------------------

// One ring of a fused system named as a monocyclic component.
struct Component {
  std::vector<int> cyc;                       // atoms in ring order
  std::vector<std::vector<int>> numberings;   // standard numberings (cyclic orders, locant 1 first)
  bool carbocycle = true, benzene = false;
  int size = 0, nN = 0, nO = 0, nS = 0, nP = 0;
  std::string name, prefix;                   // "pyridine"/"pyrido", "1,3-thiazole"/"[1,3]thiazolo"
  std::string fusion_base;                    // the name as a base component: "pyridine", "[1,3]thiazole"
  std::string hw_bare;                        // name without locants, for benzo names: "oxazole"
  std::vector<int> het_locs;                  // heteroatom locants of the standard numbering, by element O, S, N
};

// A ring as a fusion component: its standard numberings, name, prefix and base forms.
Component Namer::ring_component(const std::vector<int>& cyc) const {
  Component c;
  c.cyc = cyc;
  c.size = static_cast<int>(c.cyc.size());
  for (int a : c.cyc) {
    c.nN += at_[a].z == 7;
    c.nO += at_[a].z == 8;
    c.nS += at_[a].z == 16;
    c.nP += at_[a].z == 15;
  }
  c.carbocycle = c.nN + c.nO + c.nS + c.nP == 0;
  c.benzene = c.carbocycle && c.size == 6;
  // standard numberings: heteroatoms as low as possible, then O, S, N in that order
  std::vector<std::vector<double>> keys;
  std::vector<std::vector<int>> all;
  for (int st = 0; st < c.size; ++st)
    for (int dir : {1, -1}) {
      std::vector<int> o;
      for (int i = 0; i < c.size; ++i) o.push_back(c.cyc[((st + dir * i) % c.size + c.size) % c.size]);
      std::vector<double> key;
      {
        const int senior_z = c.nO ? 8 : c.nS ? 16 : c.nN ? 7 : 15;  // locant 1 to the senior heteroatom (P-22.2.2.1.3)
        key.push_back(c.carbocycle || at_[o[0]].z == senior_z ? 0 : 1);
      }
      for (int i = 0; i < c.size; ++i)
        if (at_[o[i]].z != 6) key.push_back(i + 1);
      for (int z : {8, 16, 7, 15})
        for (int i = 0; i < c.size; ++i)
          if (at_[o[i]].z == z) key.push_back(i + 1);
      keys.push_back(key);
      all.push_back(o);
    }
  const auto best = *std::min_element(keys.begin(), keys.end());
  for (std::size_t i = 0; i < all.size(); ++i)
    if (keys[i] == best) c.numberings.push_back(all[i]);
  const auto& o = c.numberings[0];
  std::map<int, std::vector<int>> hl;
  for (int i = 0; i < c.size; ++i)
    if (at_[o[i]].z != 6) hl[at_[o[i]].z].push_back(i + 1);
  for (int z : {8, 16, 7, 15})
    for (int l : hl[z]) c.het_locs.push_back(l);
  if (c.carbocycle) {
    c.name = c.benzene ? "benzene" : "cyclo" + alkane_stem(c.size) + "ene";
    c.prefix = c.benzene ? "benzo" : "cyclo" + alkane_stem(c.size) + "a";
    c.fusion_base = c.name;
    return c;
  }
  std::string key = std::to_string(c.size) + ":";
  for (int z : {8, 16, 7, 15})
    for (int l : hl[z]) key += std::string(1, z == 8 ? 'O' : z == 16 ? 'S' : z == 7 ? 'N' : 'P') + std::to_string(l);
  static const std::map<std::string, std::pair<const char*, const char*>> retained = {
      {"6:N1", {"pyridine", "pyrido"}},     {"6:N1N3", {"pyrimidine", "pyrimido"}}, {"6:N1N4", {"pyrazine", "pyrazino"}},
      {"6:N1N2", {"pyridazine", "pyridazino"}}, {"6:O1", {"pyran", "pyrano"}},   {"6:S1", {"thiopyran", "thiopyrano"}},
      {"5:N1", {"pyrrole", "pyrrolo"}},     {"5:O1", {"furan", "furo"}},           {"5:S1", {"thiophene", "thieno"}},
      {"5:N1N3", {"imidazole", "imidazo"}}, {"5:N1N2", {"pyrazole", "pyrazolo"}}, {"5:N1N2N3N4", {"tetrazole", "tetrazolo"}}};
  if (auto it = retained.find(key); it != retained.end()) {
    c.name = c.fusion_base = it->second.first;
    c.prefix = it->second.second;
    c.hw_bare = c.name;
    return c;
  }
  // Hantzsch-Widman name of the mancude ring
  static const char* uns_n[] = {"", "", "", "irine", "ete", "ole", "ine", "epine", "ocine"};
  static const char* uns_o[] = {"", "", "", "irene", "ete", "ole", "ine", "epine", "ocine"};
  const std::string stem = c.size == 6 && c.nP ? "inine" : c.nN ? uns_n[c.size] : uns_o[c.size];
  std::string pre, locs;
  int nhet = 0;
  for (const auto& [z, word] : std::vector<std::pair<int, std::string>>{{8, "oxa"}, {16, "thia"}, {7, "aza"}, {15, "phospha"}}) {
    if (hl[z].empty()) continue;
    for (int l : hl[z]) locs += (locs.empty() ? "" : ",") + std::to_string(l), ++nhet;
    std::string w = multiplier(static_cast<int>(hl[z].size()), false) + word;
    if (!pre.empty() && pre.back() == 'a' && starts_with_vowel(w)) pre.pop_back();
    pre += w;
  }
  if (pre.back() == 'a' && starts_with_vowel(stem)) pre.pop_back();
  c.hw_bare = pre + stem;
  c.name = (nhet > 1 ? locs + "-" : "") + c.hw_bare;
  c.fusion_base = (nhet > 1 ? "[" + locs + "]" : "") + c.hw_bare;
  std::string pf = c.hw_bare;
  if (pf.back() == 'e') pf.pop_back();
  c.prefix = (nhet > 1 ? "[" + locs + "]" : "") + pf + "o";
  return c;
}

bool Namer::fusion_bicycle(RingSys& R) {
  const std::set<int> members(R.atoms.begin(), R.atoms.end());
  // two rings sharing one bond, members C, N, O, S, ring sizes 5-8
  std::vector<int> heads;
  for (int a : R.atoms) {
    int d = 0;
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k) d += ringbond_[a][k] && members.count(at_[a].nb[k]);
    if (d == 3) heads.push_back(a);
    else if (d != 2) return false;
    if (at_[a].z != 6 && at_[a].z != 7 && at_[a].z != 8 && at_[a].z != 16 && at_[a].z != 15) return false;
  }
  if (heads.size() != 2 || order(heads[0], heads[1]) == 0) return false;
  const int f1 = heads[0], f2 = heads[1];
  std::vector<std::vector<int>> sides;  // atoms from f1 to f2, excluding both
  for (std::size_t k = 0; k < at_[f1].nb.size(); ++k) {
    if (!ringbond_[f1][k] || at_[f1].nb[k] == f2 || !members.count(at_[f1].nb[k])) continue;
    std::vector<int> path;
    int prev = f1, cur = at_[f1].nb[k];
    while (cur != f2) {
      path.push_back(cur);
      int next = -1;
      for (std::size_t j = 0; j < at_[cur].nb.size(); ++j)
        if (ringbond_[cur][j] && at_[cur].nb[j] != prev && members.count(at_[cur].nb[j])) next = at_[cur].nb[j];
      if (next < 0 || path.size() > R.atoms.size()) return false;
      prev = cur;
      cur = next;
    }
    sides.push_back(path);
  }
  if (sides.size() != 2) return false;
  for (const auto& sd : sides)
    if (sd.size() + 2 < 5 || sd.size() + 2 > 8) return fail("fused ring of fewer than five or more than eight members");

  auto make = [&](const std::vector<int>& side) {
    std::vector<int> cyc = {f1};
    cyc.insert(cyc.end(), side.begin(), side.end());
    cyc.push_back(f2);
    return ring_component(cyc);
  };
  Component A = make(sides[0]), B = make(sides[1]);

  // base component (FR-3.1): nitrogenous; senior heteroatom; larger ring; more heteroatoms; variety;
  // more heteroatoms early in O, S, N; lower heteroatom locants
  auto seniority = [](const Component& c) {
    std::vector<double> v = {double(c.nN > 0), double(c.nN == 0 && c.nO > 0), double(c.nN == 0 && c.nO == 0 && c.nS > 0),
                             double(c.nN == 0 && c.nO == 0 && c.nS == 0 && c.nP > 0), double(c.size),
                             double(c.nN + c.nO + c.nS + c.nP), double((c.nN > 0) + (c.nO > 0) + (c.nS > 0) + (c.nP > 0)),
                             double(c.nO), double(c.nS), double(c.nN), double(c.nP)};
    for (int l : c.het_locs) v.push_back(-l);
    return v;
  };
  std::string name;  // built below once the numbering is known (benzo names need it)
  bool benzo_name = false;
  const Component* base = &A;
  const Component* att = &B;
  if (A.carbocycle && B.carbocycle) {
    const Component& big = A.size >= B.size ? A : B;
    const Component& small = A.size >= B.size ? B : A;
    if (small.benzene && big.size >= 7) name = "benzo[" + std::to_string(big.size) + "]annulene";
    else if (small.size == 5 && big.size == 5) name = "pentalene";
    else if (small.size == 7 && big.size == 7) name = "heptalene";
    else return fail("fused carbocyclic system without a supported name");
  } else if ((A.benzene && !B.carbocycle) || (B.benzene && !A.carbocycle)) {
    benzo_name = true;
    base = A.benzene ? &B : &A;
    att = A.benzene ? &A : &B;
  } else {
    if (seniority(B) > seniority(A)) std::swap(base, att);
    if (seniority(*base) == seniority(*att) && base->carbocycle == att->carbocycle && !base->carbocycle && base->size == att->size &&
        base->het_locs == att->het_locs)
      ;  // identical components: either may be the base
    // fusion descriptor: earliest base letter, then lowest attached locants (FR-4.5, FR-4.6)
    std::string bestdesc;
    int bestletter = 1 << 20;
    std::vector<int> bestatt;
    for (const auto& bn : base->numberings) {
      const int m = static_cast<int>(bn.size());
      const int i1 = static_cast<int>(std::find(bn.begin(), bn.end(), f1) - bn.begin());
      const int i2 = static_cast<int>(std::find(bn.begin(), bn.end(), f2) - bn.begin());
      int letter, firstAtom, secondAtom;
      if ((i1 + 1) % m == i2) letter = i1, firstAtom = f1, secondAtom = f2;
      else letter = i2, firstAtom = f2, secondAtom = f1;
      std::vector<int> attl;
      if (!att->carbocycle) {
        std::vector<int> bestpair;
        for (const auto& an : att->numberings) {
          const int l1 = static_cast<int>(std::find(an.begin(), an.end(), firstAtom) - an.begin()) + 1;
          const int l2 = static_cast<int>(std::find(an.begin(), an.end(), secondAtom) - an.begin()) + 1;
          std::vector<int> pr = {l1, l2};
          auto sorted_pr = pr;
          std::sort(sorted_pr.begin(), sorted_pr.end());
          if (bestpair.empty()) bestpair = pr;
          else {
            auto sb = bestpair;
            std::sort(sb.begin(), sb.end());
            if (sorted_pr < sb || (sorted_pr == sb && pr < bestpair)) bestpair = pr;
          }
        }
        attl = bestpair;
      }
      bool better = letter < bestletter;
      if (letter == bestletter && !attl.empty()) {
        auto s1 = attl, s2 = bestatt;
        std::sort(s1.begin(), s1.end());
        std::sort(s2.begin(), s2.end());
        better = s1 < s2 || (s1 == s2 && attl < bestatt);
      }
      if (better) {
        bestletter = letter;
        bestatt = attl;
      }
    }
    std::string desc = "[";
    if (!bestatt.empty()) desc += std::to_string(bestatt[0]) + "," + std::to_string(bestatt[1]) + "-";
    desc += static_cast<char>('a' + bestletter);
    desc += "]";
    name = att->prefix + desc + base->fusion_base;
  }

  // peripheral numbering (FR-5.3, FR-5.4): start next to a fusion atom in either ring, clockwise;
  // fusion carbons take the previous number with a letter, fusion heteroatoms a number
  struct Num {
    std::vector<int> o;
    std::vector<std::string> L;
    std::vector<double> key;
  };
  std::vector<Num> cands;
  for (int right = 0; right < 2; ++right)
    for (int top = 0; top < 2; ++top) {
      const int ft = top ? f2 : f1, fb = top ? f1 : f2;
      auto r = sides[right], l = sides[1 - right];  // each runs f1 -> f2
      if (ft == f2) std::reverse(r.begin(), r.end()), std::reverse(l.begin(), l.end());
      // right ring from ft to fb, then left ring from fb back to ft
      std::vector<int> seq(r.begin(), r.end());
      seq.push_back(fb);
      seq.insert(seq.end(), l.rbegin(), l.rend());
      seq.push_back(ft);
      Num nm;
      nm.o = seq;
      int counter = 0;
      for (int a : seq) {
        const bool fusion = a == f1 || a == f2;
        if (fusion && at_[a].z == 6) nm.L.push_back(std::to_string(counter) + "a");
        else nm.L.push_back(std::to_string(++counter));
      }
      std::vector<std::string> het, fc, fh;
      for (std::size_t i = 0; i < seq.size(); ++i) {
        if (at_[seq[i]].z != 6) het.push_back(nm.L[i]);
        if ((seq[i] == f1 || seq[i] == f2) && at_[seq[i]].z == 6) fc.push_back(nm.L[i]);
      }
      auto add = [&](const std::vector<std::string>& v) {
        auto k = sorted_keys(v);
        nm.key.insert(nm.key.end(), k.begin(), k.end());
        nm.key.push_back(-1);
      };
      add(het);
      for (int z : {8, 16, 7, 15}) {
        std::vector<std::string> v;
        for (std::size_t i = 0; i < seq.size(); ++i)
          if (at_[seq[i]].z == z) v.push_back(nm.L[i]);
        add(v);
      }
      add(fc);
      cands.push_back(nm);
    }
  std::vector<double> best = cands[0].key;
  for (const auto& c : cands) best = std::min(best, c.key);
  for (const auto& c : cands) {
    if (c.key != best) continue;
    std::string nm = name;
    if (benzo_name) {
      const Component& h = *base;
      std::string locs;
      for (int z : {8, 16, 7, 15}) {
        std::vector<std::string> v;
        for (std::size_t i = 0; i < c.o.size(); ++i)
          if (at_[c.o[i]].z == z) v.push_back(c.L[i]);
        std::sort(v.begin(), v.end(), [](const std::string& x, const std::string& y) { return locant_key(x) < locant_key(y); });
        for (const auto& x : v) locs += (locs.empty() ? "" : ",") + x;
      }
      nm = locs + "-" + (starts_with_vowel(h.hw_bare) ? "benz" : "benzo") + h.hw_bare;
    }
    R.orders.push_back(c.o);
    R.locs.push_back(c.L);
    R.names.push_back(nm);
  }
  return true;
}

// ---- fusion nomenclature for fused systems of three or more rings (P-25.3, FR-3 to FR-5) ----------
//
// Covered: ortho- and ortho-peri-fused systems (any two rings share at most one bond; an atom may
// belong to three rings). The name is built from a base component (FR-3.1) and attached components up to second
// order; the numbering comes from the orientation rules of FR-5.2 applied to every placement of
// the rings on a hexagonal grid, then the lowest-locant rules of FR-5.4.

namespace {

struct FusionUnit {
  std::vector<int> rings;                        // ring indices, sorted
  std::vector<std::vector<int>> orders;          // atoms along the periphery from locant 1 towards 2
  std::vector<std::vector<std::string>> locs;    // locants per order
  std::set<int> atomset;
  std::string base, prefix;                      // "quinoline" / "quino", "[1]benzofuran" / "[1]benzofuro"
  bool prefix_locants = true;                    // false for benzo and monocarbocycles
  bool usable_base = true;
  std::vector<double> senior;
};

// Shapes a ring may take on the drawing grid (FR-5.1): the directions, in 45-degree units
// counter-clockwise from east, in which its sides face. A hexagon faces E, NE, NW, W, SW, SE; smaller
// rings leave directions out; a heptagon adds north or south and an octagon both.
const std::vector<std::vector<int>>& ring_shapes(int n) {
  static std::map<int, std::vector<std::vector<int>>> cache;
  auto it = cache.find(n);
  if (it != cache.end()) return it->second;
  std::vector<std::vector<int>> out;
  const std::vector<int> hex = {0, 1, 3, 4, 5, 7};
  if (n == 8) out.push_back({0, 1, 2, 3, 4, 5, 6, 7});
  else if (n == 7) {
    out.push_back({0, 1, 3, 4, 5, 6, 7});
    out.push_back({0, 1, 2, 3, 4, 5, 7});
  } else if (n == 6) {
    out.push_back(hex);
    out.push_back({1, 2, 3, 5, 6, 7});  // turned hexagon: needed above or below a seven- or eight-membered ring
  } else if (n == 5) {
    for (const std::vector<int>& six : {hex, std::vector<int>{1, 2, 3, 5, 6, 7}})
      for (int drop : six) {
        std::vector<int> v;
        for (int d : six)
          if (d != drop) v.push_back(d);
        out.push_back(v);
      }
  } else if (n == 4) {
    out.push_back({0, 2, 4, 6});
    for (std::size_t a = 0; a < 6; ++a)
      for (std::size_t b = a + 1; b < 6; ++b) {
        std::vector<int> v;
        for (std::size_t k = 0; k < 6; ++k)
          if (k != a && k != b) v.push_back(hex[k]);
        out.push_back(v);
      }
  } else if (n == 3) {
    out.push_back({1, 4, 7});
    out.push_back({0, 3, 5});
  }
  return cache[n] = out;
}

std::string bracket_locants(const std::string& name) {  // "1,3-benzothiazole" -> "[1,3]benzothiazole"
  std::size_t i = 0;
  while (i < name.size() && ((name[i] >= '0' && name[i] <= '9') || name[i] == ',')) ++i;
  if (i == 0 || i >= name.size() || name[i] != '-') return name;
  return "[" + name.substr(0, i) + "]" + name.substr(i + 1);
}

std::string fusion_prefix(const std::string& base) {  // FR-2.1 contractions, otherwise -o
  static const std::vector<std::pair<std::string, std::string>> special = {
      {"naphthalene", "naphtho"}, {"anthracene", "anthra"}, {"phenanthrene", "phenanthro"}, {"isoquinoline", "isoquino"},
      {"quinoline", "quino"},     {"benzofuran", "benzofuro"}, {"benzothiophene", "benzothieno"}, {"benzimidazole", "benzimidazo"}};
  for (const auto& [full, pre] : special)
    if (ends_with(base, full)) return base.substr(0, base.size() - full.size()) + pre;
  std::string p = base;
  if (p.back() == 'e') p.pop_back();
  return p + "o";
}

}  // namespace

// Ring analysis of a fused system: its rings, the bonds they share and what lies on the periphery.
bool Namer::build_fused(const std::vector<std::vector<int>>& rings, FusedSys& F) {
  F = FusedSys{};
  F.rings = rings;
  const int nr = static_cast<int>(rings.size());
  for (int r = 0; r < nr; ++r) {
    const auto& c = rings[r];
    if (c.size() < 3 || c.size() > 8) return fail("fused ring of more than eight members");
    for (std::size_t i = 0; i < c.size(); ++i) {
      F.rings_of[c[i]].push_back(r);
      const int a = c[i], b = c[(i + 1) % c.size()];
      ++F.bond_rings[{std::min(a, b), std::max(a, b)}];
    }
  }
  for (const auto& [a, rs] : F.rings_of) {
    F.atoms.push_back(a);
    if (rs.size() > 3) return fail("atom common to four rings at " + sym(a));
    if (at_[a].z != 6 && at_[a].z != 7 && at_[a].z != 8 && at_[a].z != 16 && at_[a].z != 15) return false;
  }
  F.radj.assign(nr, {});
  for (int r = 0; r < nr; ++r)
    for (int s2 = r + 1; s2 < nr; ++s2) {
      std::vector<int> common;
      for (int a : rings[r])
        if (std::find(rings[s2].begin(), rings[s2].end(), a) != rings[s2].end()) common.push_back(a);
      if (common.empty()) continue;
      if (common.size() != 2 || order(common[0], common[1]) == 0) return fail("rings sharing more than one bond");
      F.radj[r].push_back(static_cast<int>(F.edges.size()));
      F.radj[s2].push_back(static_cast<int>(F.edges.size()));
      F.edges.push_back({r, s2, common[0], common[1]});
    }
  for (const auto& kv : F.bond_rings)
    if (kv.second > 2) return false;
  for (int a : F.atoms) {
    bool peripheral = false;
    for (int b : at_[a].nb) {
      auto it = F.bond_rings.find({std::min(a, b), std::max(a, b)});
      if (it != F.bond_rings.end() && it->second == 1) peripheral = true;
    }
    if (!peripheral) {
      if (at_[a].z != 6) return fail("interior heteroatom in a fused ring system at " + sym(a));
      F.interior.insert(a);
    }
  }
  // the rings must hang together through shared bonds
  std::vector<char> seen(nr, 0);
  std::vector<int> q = {0};
  seen[0] = 1;
  for (std::size_t h = 0; h < q.size(); ++h)
    for (int ei : F.radj[q[h]]) {
      const int o2 = F.edges[ei].r == q[h] ? F.edges[ei].s : F.edges[ei].r;
      if (!seen[o2]) seen[o2] = 1, q.push_back(o2);
    }
  return static_cast<int>(q.size()) == nr;
}

// Peripheral numbering of a fused system (FR-5.2 to FR-5.5): every placement of the rings on the
// eight-direction grid is scored by the orientation rules; the best ones are numbered clockwise
// from the upper right ring and the lowest-locant rules pick among them.
bool Namer::number_fused(const FusedSys& F, std::vector<FusedNumbering>& out) {
  out.clear();
  const int nr = static_cast<int>(F.rings.size());
  auto bond_count = [&](int a, int b) {
    auto it = F.bond_rings.find({std::min(a, b), std::max(a, b)});
    return it == F.bond_rings.end() ? 0 : it->second;
  };
  auto fusion_atom = [&](int a) { return F.rings_of.at(a).size() >= 2; };
  static const int vx[8] = {2, 1, 0, -1, -2, -1, 0, 1}, vy[8] = {0, 1, 2, 1, 0, -1, -2, -1};
  std::vector<int> bfs = {0}, parent_edge(nr, -1);
  {
    std::vector<char> seen(nr, 0);
    seen[0] = 1;
    for (std::size_t h = 0; h < bfs.size(); ++h)
      for (int ei : F.radj[bfs[h]]) {
        const int o2 = F.edges[ei].r == bfs[h] ? F.edges[ei].s : F.edges[ei].r;
        if (!seen[o2]) seen[o2] = 1, parent_edge[o2] = ei, bfs.push_back(o2);
      }
  }
  std::size_t n_peripheral = 0;
  for (int a : F.atoms) n_peripheral += !F.interior.count(a);
  long evaluations = 0;
  const long kMaxEvaluations = 2000000;
  std::vector<FusedNumbering> found;
  std::vector<double> best_orient;
  std::vector<std::vector<int>> cw(nr), dir(nr);
  std::vector<int> px(nr), py(nr);
  auto bond_index = [&](int r, int u, int v) {
    const int n = static_cast<int>(cw[r].size());
    for (int i = 0; i < n; ++i)
      if ((cw[r][i] == u && cw[r][(i + 1) % n] == v) || (cw[r][i] == v && cw[r][(i + 1) % n] == u)) return i;
    return -1;
  };

  auto evaluate = [&]() {
    std::map<std::pair<int, int>, int> at;
    for (int r = 0; r < nr; ++r)
      if (!at.emplace(std::make_pair(px[r], py[r]), r).second) return;  // overlapping rings
    // every shared bond must face its neighbour from both sides (closes the cycles of peri-fusion)
    for (const auto& e : F.edges) {
      const int i = bond_index(e.r, e.u, e.v), j = bond_index(e.s, e.u, e.v);
      if (i < 0 || j < 0) return;
      const int d = dir[e.r][i];
      if (px[e.s] != px[e.r] + vx[d] || py[e.s] != py[e.r] + vy[d] || dir[e.s][j] != (d + 4) % 8) return;
      if (cw[e.r][i] == cw[e.s][j]) return;  // the two rings must run through the bond in opposite senses
    }
    auto fused = [&](int p, int q2) {
      for (int ei : F.radj[p])
        if (F.edges[ei].r == q2 || F.edges[ei].s == q2) return true;
      return false;
    };
    auto east_of = [&](int r) {
      auto it = at.find({px[r] + 2, py[r]});
      return it != at.end() && fused(r, it->second) ? it->second : -1;
    };
    // FR-5.2: (a) rings in the horizontal row, (b) upper right quadrant, (c) lower left, (d) above the row
    std::vector<double> orient;
    for (int r = 0; r < nr; ++r) {
      {
        auto it = at.find({px[r] - 2, py[r]});
        if (it != at.end() && fused(r, it->second)) continue;
      }
      int len = 1;
      for (int q2 = east_of(r); q2 >= 0; q2 = east_of(q2)) ++len;
      const double cx = px[r] + (len - 1), cy = py[r];
      double ur = 0, ll = 0, above = 0;
      for (int q2 = 0; q2 < nr; ++q2) {
        const double dx = px[q2] - cx, dy = py[q2] - cy;
        ur += (dx > 0 ? 1 : dx == 0 ? 0.5 : 0) * (dy > 0 ? 1 : dy == 0 ? 0.5 : 0);
        ll += (dx < 0 ? 1 : dx == 0 ? 0.5 : 0) * (dy < 0 ? 1 : dy == 0 ? 0.5 : 0);
        above += dy > 0;
      }
      std::vector<double> t = {double(len), ur, -ll, above};
      if (orient.empty() || t > orient) orient = t;
    }
    if (!best_orient.empty() && orient < best_orient) return;
    // FR-5.3: start in the ring farthest up and to the right, at its most counter-clockwise free atom
    int top = 0;
    for (int r = 1; r < nr; ++r)
      if (py[r] > py[top] || (py[r] == py[top] && px[r] > px[top])) top = r;
    const int n = static_cast<int>(cw[top].size());
    int pick = -1;
    double pick_d = 1e9, pick_a = 0;
    for (int i = 0; i < n; ++i) {
      if (fusion_atom(cw[top][i])) continue;
      const int pb = (i + n - 1) % n;  // atom i lies between bond i-1 and bond i, clockwise
      const int gap = ((dir[top][pb] - dir[top][i]) % 8 + 8) % 8;
      double ang = 45.0 * dir[top][pb] - 22.5 * gap;
      while (ang > 180) ang -= 360;
      while (ang <= -180) ang += 360;
      const double d = std::abs(ang - 30.0);
      if (d < pick_d - 1e-9 || (std::abs(d - pick_d) < 1e-9 && ang > pick_a)) pick = i, pick_d = d, pick_a = ang;
    }
    if (pick < 0) return;
    for (int guard = 0; guard < n && !fusion_atom(cw[top][(pick + n - 1) % n]); ++guard) pick = (pick + n - 1) % n;
    FusedNumbering nm;
    int a = cw[top][pick];
    int counter = 0, letter = 0;
    for (std::size_t steps = 0; steps < n_peripheral; ++steps) {
      nm.o.push_back(a);
      if (fusion_atom(a) && at_[a].z == 6) {
        nm.L.push_back(std::to_string(counter) + std::string(1, static_cast<char>('a' + letter++)));
      } else {
        nm.L.push_back(std::to_string(++counter));
        letter = 0;
      }
      int next = -1;
      for (int r : F.rings_of.at(a)) {
        const auto& c = cw[r];
        const std::size_t i = static_cast<std::size_t>(std::find(c.begin(), c.end(), a) - c.begin());
        const int b = c[(i + 1) % c.size()];
        if (bond_count(a, b) == 1) next = b;
      }
      if (next < 0) return;
      a = next;
    }
    if (a != nm.o[0] || std::set<int>(nm.o.begin(), nm.o.end()).size() != n_peripheral) return;
    // FR-5.5: interior atoms continue the letters of the last peripheral locant, nearest first
    {
      std::set<int> todo = F.interior;
      while (!todo.empty()) {
        int pick_atom = -1, pick_rank = -1, ties = 0;
        for (int x : todo) {
          int rank = -1;
          for (int b : at_[x].nb) {
            const auto it = std::find(nm.o.begin(), nm.o.end(), b);
            if (it != nm.o.end() && bond_count(x, b) > 0) rank = std::max(rank, static_cast<int>(it - nm.o.begin()));
          }
          if (rank > pick_rank) pick_atom = x, pick_rank = rank, ties = 1;
          else if (rank == pick_rank) ++ties;
        }
        if (pick_atom < 0 || pick_rank < 0 || ties != 1) return;
        nm.o.push_back(pick_atom);
        nm.L.push_back(std::to_string(counter) + std::string(1, static_cast<char>('a' + letter++)));
        todo.erase(pick_atom);
      }
    }
    if (letter > 25) return;
    auto add = [&](const std::vector<std::string>& v) {
      const auto k = sorted_keys(v);
      nm.key.insert(nm.key.end(), k.begin(), k.end());
      nm.key.push_back(-1);
    };
    std::vector<std::string> het, fc;
    for (std::size_t i = 0; i < nm.o.size(); ++i) {
      if (at_[nm.o[i]].z != 6) het.push_back(nm.L[i]);
      else if (fusion_atom(nm.o[i])) fc.push_back(nm.L[i]);
    }
    add(het);
    for (int z : {8, 16, 7, 15}) {
      std::vector<std::string> v;
      for (std::size_t i = 0; i < nm.o.size(); ++i)
        if (at_[nm.o[i]].z == z) v.push_back(nm.L[i]);
      add(v);
    }
    add(fc);
    if (best_orient.empty() || orient > best_orient) {
      best_orient = orient;
      found.clear();
    }
    for (const auto& f : found)
      if (f.o == nm.o) return;
    found.push_back(std::move(nm));
  };

  auto place = [&](auto&& self, std::size_t idx) -> void {
    if (evaluations > kMaxEvaluations) return;
    if (idx == bfs.size()) {
      ++evaluations;
      evaluate();
      return;
    }
    const int r = bfs[idx];
    const int n = static_cast<int>(cw[r].size());
    int e = 0;
    const int de = dir[r][0];  // entry direction, left in slot 0 by the parent (or the root loop)
    if (parent_edge[r] >= 0) e = bond_index(r, F.edges[parent_edge[r]].u, F.edges[parent_edge[r]].v);
    std::set<std::vector<int>> tried;  // shapes that send the fused neighbours the same way are one case
    const std::vector<int> saved = dir[r];
    for (const auto& shape : ring_shapes(n)) {
      const auto pos_it = std::find(shape.begin(), shape.end(), de);
      if (pos_it == shape.end()) continue;
      const int p = static_cast<int>(pos_it - shape.begin());
      std::vector<int> d(n);
      for (int k = 0; k < n; ++k) d[(e + k) % n] = shape[((p - k) % n + n) % n];  // clockwise: decreasing angle
      {
        std::vector<int> sig;
        for (int i = 0; i < n; ++i)
          if (bond_count(cw[r][i], cw[r][(i + 1) % n]) == 2) sig.push_back(d[i]);
        if (!tried.insert(sig).second) continue;
      }
      dir[r] = d;
      bool clash = false;
      for (int ei : F.radj[r]) {
        const int c = F.edges[ei].r == r ? F.edges[ei].s : F.edges[ei].r;
        if (parent_edge[c] != ei) continue;
        const int bi = bond_index(r, F.edges[ei].u, F.edges[ei].v);
        px[c] = px[r] + vx[d[bi]];
        py[c] = py[r] + vy[d[bi]];
        const int from = cw[r][(bi + 1) % n], to = cw[r][bi];
        std::vector<int> cc = F.rings[c];
        const std::size_t m = cc.size();
        std::size_t fi = static_cast<std::size_t>(std::find(cc.begin(), cc.end(), from) - cc.begin());
        if (cc[(fi + 1) % m] != to) std::reverse(cc.begin(), cc.end());
        cw[c] = cc;
        dir[c] = {(d[bi] + 4) % 8};
      }
      if (!clash) self(self, idx + 1);
    }
    dir[r] = saved;
  };
  for (int mirror = 0; mirror < 2; ++mirror)
    for (int d0 = 0; d0 < 8; ++d0) {
      cw[0] = F.rings[0];
      if (mirror) std::reverse(cw[0].begin(), cw[0].end());
      px[0] = py[0] = 0;
      dir[0] = {d0};
      place(place, 0);
    }
  if (evaluations > kMaxEvaluations) return fail("fused ring system with too many possible orientations");
  if (found.empty()) return fail("fused ring system that cannot be oriented on the grid");
  std::vector<double> best_key = found[0].key;
  for (const auto& f : found) best_key = std::min(best_key, f.key);
  for (auto& f : found)
    if (f.key == best_key) out.push_back(std::move(f));
  return true;
}

bool Namer::fusion_polycycle(RingSys& R) {
  if (sssr_.empty())
    for (const Ring& rg : find_sssr(st_, sp_)) {
      std::vector<int> cyc;
      for (std::uint32_t g : rg.atoms) {
        auto it = local_.find(g);
        if (it == local_.end()) break;
        cyc.push_back(it->second);
      }
      if (cyc.size() == rg.atoms.size()) sssr_.push_back(cyc);
    }
  const std::set<int> in(R.atoms.begin(), R.atoms.end());
  std::vector<std::vector<int>> rings;
  for (const auto& cyc : sssr_)
    if (std::all_of(cyc.begin(), cyc.end(), [&](int a) { return in.count(a) > 0; })) rings.push_back(cyc);
  const int nr = static_cast<int>(rings.size());
  if (nr != R.nrings || nr < 3) return false;
  FusedSys F;
  if (!build_fused(rings, F)) return false;
  const auto& edges = F.edges;
  const auto& radj = F.radj;
  std::vector<FusedNumbering> numberings;
  if (!number_fused(F, numberings)) return false;

  // ---- components -------------------------------------------------------------------------------
  // FR-3.1: nitrogen; else the senior heteroatom; more rings; larger rings at the first difference;
  // more heteroatoms; more variety; more O, S, N in that order; lower heteroatom locants
  auto unit_seniority = [&](const std::map<int, std::string>& loc, const std::vector<int>& rs, double rank) {
    int nN = 0, nO = 0, nS = 0, nP = 0;
    std::vector<double> hl, sizes;
    for (const auto& [a, l] : loc) {
      nN += at_[a].z == 7;
      nO += at_[a].z == 8;
      nS += at_[a].z == 16;
      nP += at_[a].z == 15;
      if (at_[a].z != 6) hl.push_back(-locant_key(l));
    }
    for (int r : rs) sizes.push_back(static_cast<double>(rings[r].size()));
    std::sort(sizes.rbegin(), sizes.rend());
    std::sort(hl.rbegin(), hl.rend());
    std::vector<double> v = {double(nN > 0), double(nN == 0 && nO > 0), double(nN == 0 && nO == 0 && nS > 0),
                             double(nN == 0 && nO == 0 && nS == 0 && nP > 0), double(rs.size())};
    v.insert(v.end(), sizes.begin(), sizes.end());
    v.insert(v.end(), {double(nN + nO + nS + nP), double((nN > 0) + (nO > 0) + (nS > 0) + (nP > 0)), double(nO), double(nS), double(nN),
                       double(nP)});
    v.insert(v.end(), hl.begin(), hl.end());
    v.push_back(rank);
    return v;
  };
  std::vector<FusionUnit> units;
  // adds one numbering of a unit: its periphery from locant 1 towards locant 2
  auto add_order = [&](FusionUnit& u, const std::map<int, std::string>& loc) {
    std::vector<std::vector<int>> rs;
    for (int r : u.rings) rs.push_back(rings[r]);
    std::map<std::pair<int, int>, int> cnt;
    for (const auto& c : rs)
      for (std::size_t i = 0; i < c.size(); ++i) ++cnt[{std::min(c[i], c[(i + 1) % c.size()]), std::max(c[i], c[(i + 1) % c.size()])}];
    std::map<int, std::vector<int>> pn;  // peripheral neighbours
    for (const auto& [b, k] : cnt)
      if (k == 1) pn[b.first].push_back(b.second), pn[b.second].push_back(b.first);
    int start = -1;
    for (const auto& [a, l] : loc)
      if (l == "1") start = a;
    if (start < 0 || pn[start].size() != 2) return false;
    auto lk = [&](int a) {
      auto it = loc.find(a);
      return it == loc.end() ? 1e9 : locant_key(it->second);
    };
    int prev = start, cur = lk(pn[start][0]) <= lk(pn[start][1]) ? pn[start][0] : pn[start][1];
    std::vector<int> o = {start};
    while (cur != start) {
      if (pn[cur].size() != 2 || o.size() > loc.size()) return false;
      o.push_back(cur);
      const int next = pn[cur][0] == prev ? pn[cur][1] : pn[cur][0];
      prev = cur;
      cur = next;
    }
    std::vector<std::string> l;
    for (int a : o) {
      auto it = loc.find(a);
      if (it == loc.end()) return false;
      l.push_back(it->second);
    }
    u.orders.push_back(o);
    u.locs.push_back(l);
    for (const auto& kv : loc) u.atomset.insert(kv.first);
    return true;
  };
  auto rings_inside = [&](const std::set<int>& os, std::vector<int>& rs) {
    std::set<int> covered;
    rs.clear();
    for (int r = 0; r < nr; ++r)
      if (std::all_of(rings[r].begin(), rings[r].end(), [&](int a) { return os.count(a) > 0; })) {
        rs.push_back(r);
        covered.insert(rings[r].begin(), rings[r].end());
      }
    return covered == os;
  };
  static const std::map<std::string, double> carbo_rank = {
      {"indene", 2}, {"naphthalene", 3}, {"azulene", 4}, {"acenaphthylene", 9}, {"fluorene", 10}, {"phenalene", 11}, {"phenanthrene", 12},
      {"anthracene", 13}, {"fluoranthene", 14}, {"pyrene", 18}, {"perylene", 23}};
  for (int r = 0; r < nr; ++r) {
    const Component c = ring_component(rings[r]);
    FusionUnit u;
    u.rings = {r};
    u.base = c.fusion_base;
    u.prefix = c.prefix;
    u.prefix_locants = !c.carbocycle;
    u.usable_base = !c.carbocycle;
    std::map<int, std::string> first;
    for (const auto& o : c.numberings) {
      std::map<int, std::string> loc;
      for (int i = 0; i < c.size; ++i) loc[o[i]] = std::to_string(i + 1);
      if (first.empty()) first = loc;
      add_order(u, loc);
    }
    u.senior = unit_seniority(first, u.rings, 0);
    units.push_back(std::move(u));
  }
  {
    const auto& qs = template_queries();
    SmartsMatchOptions mo;
    mo.unique = false;
    mo.max_matches = 0;
    for (std::size_t t = 0; t < qs.size(); ++t) {
      const Template& T = templates()[t];
      const int k = pattern_rings(T.pattern);
      const std::string tn = T.name;
      if (!T.mancude || !qs[t].valid() || k < 2 || k >= nr || tn.find('[') != std::string::npos || tn == "purine") continue;
      std::map<std::vector<int>, std::size_t> by_rings;
      for (const auto& mt : match_smarts(st_, sp_, qs[t], mo)) {
        std::map<int, std::string> loc;
        bool ok = true;
        for (std::size_t i = 0; i < mt.size() && ok; ++i) {
          auto it = local_.find(mt[i]);
          if (it == local_.end() || !in.count(it->second)) ok = false;
          else loc[it->second] = T.locants[i];
        }
        if (!ok || loc.size() != mt.size()) continue;
        std::set<int> os;
        for (const auto& kv : loc) os.insert(kv.first);
        std::vector<int> rs;
        if (!rings_inside(os, rs) || static_cast<int>(rs.size()) != k) continue;
        auto it = by_rings.find(rs);
        if (it == by_rings.end()) {
          FusionUnit u;
          u.rings = rs;
          u.base = bracket_locants(tn);
          u.prefix = fusion_prefix(u.base);
          auto rk = carbo_rank.find(tn);
          u.senior = unit_seniority(loc, rs, rk == carbo_rank.end() ? 0 : rk->second);
          it = by_rings.emplace(rs, units.size()).first;
          units.push_back(std::move(u));
        }
        add_order(units[it->second], loc);
      }
    }
  }
  // benzo-heteromonocycles are components in their own right (P-25.2.2.4): "[1,5]benzothiazepine"
  for (const auto& e : edges) {
    auto benzene = [&](int r) { return rings[r].size() == 6 && std::all_of(rings[r].begin(), rings[r].end(), [&](int a) { return at_[a].z == 6; }); };
    if (benzene(e.r) == benzene(e.s)) continue;
    const int het = benzene(e.r) ? e.s : e.r;
    if (std::all_of(rings[het].begin(), rings[het].end(), [&](int a) { return at_[a].z == 6; })) continue;
    std::vector<int> rs = {std::min(e.r, e.s), std::max(e.r, e.s)};
    bool known = false;
    for (const auto& u : units) known |= u.rings == rs;
    if (known) continue;
    // the pair must be ortho-fused on its own: no atom of the pair in both rings beyond the shared bond
    RingSys T;
    T.nrings = 2;
    std::set<int> ta(rings[e.r].begin(), rings[e.r].end());
    ta.insert(rings[e.s].begin(), rings[e.s].end());
    T.atoms.assign(ta.begin(), ta.end());
    const std::string saved_err = err_;
    const bool ok = fusion_bicycle(T);
    err_ = saved_err;
    if (!ok || T.names.empty() || T.names[0].find('[') != std::string::npos) continue;
    FusionUnit u;
    u.rings = rs;
    u.base = bracket_locants(T.names[0]);
    u.prefix = fusion_prefix(u.base);
    std::map<int, std::string> first;
    for (std::size_t i = 0; i < T.orders.size(); ++i) {
      std::map<int, std::string> loc;
      for (std::size_t j = 0; j < T.orders[i].size(); ++j) loc[T.orders[i][j]] = T.locs[i][j];
      if (first.empty()) first = loc;
      add_order(u, loc);
    }
    if (u.orders.empty()) continue;
    u.senior = unit_seniority(first, rs, 0);
    units.push_back(std::move(u));
  }
  // retained peri-fused hydrocarbons (P-25.1.1), numbered by the general rules
  {
    static const std::vector<std::pair<const char*, const char*>> peri = {
        {"acenaphthylene", "c1cc2cccc3cccc1c23"}, {"phenalene", "c1cc2cccc3cccc(c1)c23"}, {"pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"},
        {"fluoranthene", "c1ccc2c(c1)c1cccc3cccc2c13"}, {"perylene", "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54"}};
    static const std::vector<SmartsQuery> pq = [] {
      std::vector<SmartsQuery> v;
      for (const auto& p : peri) {
        Diagnostics d;
        v.push_back(parse_smarts(generic_smarts(p.second), d));
      }
      return v;
    }();
    SmartsMatchOptions mo;
    mo.max_matches = 0;
    for (std::size_t t = 0; t < peri.size(); ++t) {
      if (!pq[t].valid() || !F.interior.size()) continue;
      std::set<std::vector<int>> done;
      for (const auto& mt : match_smarts(st_, sp_, pq[t], mo)) {
        std::set<int> os;
        bool ok = true;
        for (std::uint32_t g : mt) {
          auto it = local_.find(g);
          if (it == local_.end() || !in.count(it->second)) ok = false;
          else os.insert(it->second);
        }
        std::vector<int> rs;
        if (!ok || !rings_inside(os, rs) || static_cast<int>(rs.size()) != pattern_rings(peri[t].second) || !done.insert(rs).second) continue;
        std::vector<std::vector<int>> sub;
        for (int r : rs) sub.push_back(rings[r]);
        FusedSys G;
        std::vector<FusedNumbering> nums;
        const std::string saved_err = err_;
        const bool good = build_fused(sub, G) && number_fused(G, nums);
        err_ = saved_err;
        if (!good) continue;
        FusionUnit u;
        u.rings = rs;
        u.base = peri[t].first;
        u.prefix = fusion_prefix(u.base);
        std::map<int, std::string> first;
        for (const auto& nmb : nums) {
          std::map<int, std::string> loc;
          for (std::size_t j = 0; j < nmb.o.size(); ++j) loc[nmb.o[j]] = nmb.L[j];
          if (first.empty()) first = loc;
          add_order(u, loc);
        }
        if (u.orders.empty()) continue;
        u.senior = unit_seniority(first, rs, carbo_rank.at(peri[t].first));
        units.push_back(std::move(u));
      }
    }
  }

  // base component (FR-3.1)
  int base = -1;
  for (std::size_t i = 0; i < units.size(); ++i)
    if (units[i].usable_base && !units[i].orders.empty() && (base < 0 || units[i].senior > units[static_cast<std::size_t>(base)].senior))
      base = static_cast<int>(i);
  if (base < 0) return fail("fused carbocyclic system without a named base component");
  const FusionUnit& B = units[static_cast<std::size_t>(base)];
  std::vector<char> in_base(nr, 0);
  for (int r : B.rings) in_base[r] = 1;

  // attached components: first-order units on the base, second-order units on those
  struct Attached {
    int unit = -1;
    std::vector<Attached> kids;
  };
  auto shared = [&](const FusionUnit& a, const FusionUnit& b) {
    int k = 0;
    for (int x : a.atomset) k += b.atomset.count(x) > 0;
    return k;
  };
  auto group_of = [&](int start, const std::vector<char>& blocked) {
    std::vector<int> g = {start};
    std::vector<char> seen(nr, 0);
    seen[start] = 1;
    for (std::size_t h = 0; h < g.size(); ++h)
      for (int ei : radj[g[h]]) {
        const int o2 = edges[ei].r == g[h] ? edges[ei].s : edges[ei].r;
        if (!seen[o2] && !blocked[o2]) seen[o2] = 1, g.push_back(o2);
      }
    std::sort(g.begin(), g.end());
    return g;
  };
  auto cover = [&](const std::vector<int>& group, Attached& out) -> bool {
    std::vector<int> cands;
    for (std::size_t i = 0; i < units.size(); ++i) {
      const auto& ur = units[i].rings;
      if (units[i].orders.empty() || !std::includes(group.begin(), group.end(), ur.begin(), ur.end())) continue;
      if (shared(units[i], B) < 2) continue;
      cands.push_back(static_cast<int>(i));
    }
    std::sort(cands.begin(), cands.end(), [&](int x, int y) {
      if (units[x].rings.size() != units[y].rings.size()) return units[x].rings.size() > units[y].rings.size();
      return units[x].senior > units[y].senior;
    });
    for (int ci : cands) {
      const FusionUnit& u = units[static_cast<std::size_t>(ci)];
      Attached a;
      a.unit = ci;
      std::vector<char> blocked(nr, 1);
      for (int r : group) blocked[r] = 0;
      for (int r : u.rings) blocked[r] = 1;
      bool ok = true;
      std::vector<char> seen(nr, 0);
      for (int r : group) {
        if (blocked[r] || seen[r]) continue;
        const std::vector<int> sub = group_of(r, blocked);
        for (int x : sub) seen[x] = 1;
        int kid = -1;
        for (std::size_t i = 0; i < units.size(); ++i)
          if (units[i].rings == sub && !units[i].orders.empty()) kid = static_cast<int>(i);
        // a second-order unit sits on its first-order unit only
        if (kid < 0 || shared(units[static_cast<std::size_t>(kid)], u) < 2 || shared(units[static_cast<std::size_t>(kid)], B) >= 2) {
          ok = false;
          break;
        }
        Attached k;
        k.unit = kid;
        a.kids.push_back(k);
      }
      if (!ok) continue;
      out = a;
      return true;
    }
    return false;
  };
  std::vector<Attached> first;
  {
    std::vector<char> seen(nr, 0);
    for (int r = 0; r < nr; ++r) {
      if (in_base[r] || seen[r]) continue;
      const std::vector<int> group = group_of(r, in_base);
      for (int x : group) seen[x] = 1;
      Attached a;
      if (!cover(group, a)) return fail("fused ring system needing attached components beyond second order");
      first.push_back(a);
    }
  }

  // The atoms unit `u` shares with the periphery `op` of its parent, as a run of consecutive positions.
  auto shared_run = [&](const FusionUnit& u, const std::vector<int>& op, int& startpos, int& len) {
    const int n = static_cast<int>(op.size());
    std::vector<char> c(n, 0);
    len = 0;
    for (int i = 0; i < n; ++i) c[i] = u.atomset.count(op[i]) > 0, len += c[i];
    if (len < 2 || len >= n) return false;
    startpos = -1;
    for (int i = 0; i < n; ++i)
      if (c[i] && !c[(i + n - 1) % n]) {
        if (startpos >= 0) return false;  // two separate runs
        startpos = i;
      }
    return startpos >= 0;
  };
  struct Piece {
    std::string prefix, kids, locs, letters;
    int letter = 0;
    std::vector<double> key;
  };
  auto cite = [&](const std::vector<std::string>& ls, const std::string& mark) {  // "1,8" from 1, 8a, 8
    std::string out;
    for (const auto& l : ls) {
      const bool lettered = !l.empty() && l.back() >= 'a' && l.back() <= 'z';
      if (ls.size() > 2 && lettered) continue;
      out += (out.empty() ? "" : ",") + l + mark;
    }
    return out;
  };
  auto loc_keys = [&](const std::vector<std::string>& ls) {
    std::vector<double> k;
    for (const auto& l : ls) {
      const bool lettered = !l.empty() && l.back() >= 'a' && l.back() <= 'z';
      if (ls.size() > 2 && lettered) continue;
      k.push_back(locant_key(l));
    }
    std::vector<double> s2 = k;
    std::sort(s2.begin(), s2.end());
    s2.insert(s2.end(), k.begin(), k.end());
    return s2;
  };
  auto describe = [&](const Attached& a, const std::vector<int>& ob, Piece& out) -> bool {
    const FusionUnit& u = units[static_cast<std::size_t>(a.unit)];
    const int nb = static_cast<int>(ob.size());
    int sp0 = 0, len = 0;
    if (!shared_run(u, ob, sp0, len)) return false;
    if (sp0 + len - 2 > 25) return false;
    bool have = false;
    for (std::size_t oi = 0; oi < u.orders.size(); ++oi) {
      const auto& oa = u.orders[oi];
      const auto& la = u.locs[oi];
      Piece p;
      p.prefix = u.prefix;
      p.letter = sp0;
      for (int k = 0; k + 1 < len; ++k) p.letters += static_cast<char>('a' + (sp0 + k) % nb);
      std::vector<std::string> ls;
      bool ok = true;
      for (int k = 0; k < len && ok; ++k) {
        const auto it = std::find(oa.begin(), oa.end(), ob[(sp0 + k) % nb]);
        if (it == oa.end()) ok = false;
        else ls.push_back(la[static_cast<std::size_t>(it - oa.begin())]);
      }
      if (!ok) continue;
      p.key = loc_keys(ls);
      if (u.prefix_locants || !a.kids.empty()) p.locs = cite(ls, "");
      for (const Attached& k : a.kids) {
        const FusionUnit& ku = units[static_cast<std::size_t>(k.unit)];
        int s1 = 0, l1 = 0;
        if (!shared_run(ku, oa, s1, l1)) { ok = false; break; }
        const int na = static_cast<int>(oa.size());
        std::vector<int> run;
        for (int j = 0; j < l1; ++j) run.push_back(oa[(s1 + j) % na]);
        std::vector<std::string> ms;
        for (int x : run) ms.push_back(la[static_cast<std::size_t>(std::find(oa.begin(), oa.end(), x) - oa.begin())]);
        if (locant_key(ms.front()) > locant_key(ms.back())) {  // cite the first-order locants ascending
          std::reverse(run.begin(), run.end());
          std::reverse(ms.begin(), ms.end());
        }
        std::string bestk;
        std::vector<double> bestkk;
        for (std::size_t ki = 0; ki < ku.orders.size(); ++ki) {
          std::vector<std::string> ks;
          bool good = true;
          for (int x : run) {
            const auto it = std::find(ku.orders[ki].begin(), ku.orders[ki].end(), x);
            if (it == ku.orders[ki].end()) { good = false; break; }
            ks.push_back(ku.locs[ki][static_cast<std::size_t>(it - ku.orders[ki].begin())]);
          }
          if (!good) continue;
          const auto kk = loc_keys(ks);
          if (bestkk.empty() || kk < bestkk) bestkk = kk, bestk = cite(ks, "'");
        }
        if (bestkk.empty()) { ok = false; break; }
        p.kids += ku.prefix + "[" + (ku.prefix_locants ? bestk + ":" : "") + cite(ms, "") + "]";
        const auto mk = loc_keys(ms);
        p.key.insert(p.key.end(), mk.begin(), mk.end());
        p.key.insert(p.key.end(), bestkk.begin(), bestkk.end());
      }
      if (!ok) continue;
      if (!have || p.key < out.key) out = p, have = true;
    }
    return have;
  };
  std::string name;
  std::vector<double> name_key;
  for (const auto& ob : B.orders) {
    std::vector<Piece> ps;
    bool ok = true;
    for (const Attached& a : first) {
      Piece p;
      if (!describe(a, ob, p)) { ok = false; break; }
      ps.push_back(p);
    }
    if (!ok) continue;
    {
      // Peri-fusion is limited to one attached component sharing two bonds with the base: names with
      // several peri-fused components or longer shared paths could not be checked against a second
      // implementation, so they are refused rather than risked.
      int peri = 0;
      for (const auto& p : ps) {
        if (p.letters.size() > 2) return fail("peri-fused component sharing more than two bonds with its base");
        peri += p.letters.size() == 2;
      }
      if (peri > 1) return fail("fused ring system with several peri-fused components");
    }
    std::stable_sort(ps.begin(), ps.end(), [](const Piece& x, const Piece& y) {
      const std::string kx = alpha_key(x.prefix), ky = alpha_key(y.prefix);
      if (kx != ky) return kx < ky;
      if (x.kids != y.kids) return x.kids < y.kids;
      return x.letter < y.letter;
    });
    std::vector<double> key;
    {
      std::vector<double> letters;
      for (const auto& p : ps) letters.push_back(p.letter);
      std::sort(letters.begin(), letters.end());
      key = letters;
      for (const auto& p : ps) key.push_back(p.letter);
      for (const auto& p : ps) key.insert(key.end(), p.key.begin(), p.key.end());
    }
    std::string nm;
    for (std::size_t i = 0; i < ps.size();) {
      std::size_t j = i;
      while (j < ps.size() && ps[j].prefix == ps[i].prefix && ps[j].kids.empty() && ps[i].kids.empty()) ++j;
      if (j == i) j = i + 1;
      const int cnt = static_cast<int>(j - i);
      std::string desc;
      for (std::size_t q = i; q < j; ++q) {
        std::string locs = ps[q].locs;
        if (!locs.empty() && q > i) {  // primes on the locants of repeated components
          const std::string primes(q - i, '\'');
          std::string marked;
          std::size_t pos = 0;
          while (pos <= locs.size()) {
            std::size_t comma = locs.find(',', pos);
            if (comma == std::string::npos) comma = locs.size();
            marked += (marked.empty() ? "" : ",") + locs.substr(pos, comma - pos) + primes;
            pos = comma + 1;
          }
          locs = marked;
        }
        if (q > i) desc += ps[i].locs.empty() ? "," : ":";
        desc += locs.empty() ? ps[q].letters : locs + "-" + ps[q].letters;
      }
      nm += ps[i].kids + (cnt > 1 ? multiplier(cnt, false) : "") + ps[i].prefix + "[" + desc + "]";
      i = j;
    }
    nm += B.base;
    if (name.empty() || key < name_key) name = nm, name_key = key;
  }
  if (name.empty()) return fail("fusion descriptor could not be formed");
  {  // retained names of polycyclic hydrocarbons (P-25.1.1); their numbering follows the general rules
    static const std::map<std::string, std::string> retained = {
        {"benzo[b]anthracene", "tetracene"}, {"benzo[a]anthracene", "tetraphene"}, {"benzo[a]phenanthrene", "chrysene"},
        {"naphtho[2,3-b]anthracene", "pentacene"}, {"benzo[l]phenanthrene", "triphenylene"}, {"dibenzo[a,i]phenanthrene", "picene"}};
    if (auto it = retained.find(name); it != retained.end()) name = it->second;
  }
  for (const auto& f : numberings) {
    R.orders.push_back(f.o);
    R.locs.push_back(f.L);
    R.names.push_back(name);
  }
  return true;
}

// Polycyclic von Baeyer systems (P-23.2.5): main ring with the most atoms, the largest main bridge
// dividing it as evenly as possible, then secondary bridges (independent before dependent, larger
// first) with superscript locants as low as possible. Used for bridged systems and for fused
// systems that have no fusion name here.
bool Namer::von_baeyer_poly(RingSys& R) {
  const std::vector<int>& V = R.atoms;
  const int n = static_cast<int>(V.size());
  if (n > 40 || R.nrings > 8) return fail("ring system too large for von Baeyer nomenclature");
  std::map<int, int> idx;
  for (int i = 0; i < n; ++i) idx[V[i]] = i;
  std::vector<std::vector<int>> adj(n);
  int nbonds = 0;
  for (int i = 0; i < n; ++i) {
    const int a = V[i];
    if (at_[a].z != 6 && at_[a].z != 7 && at_[a].z != 8 && at_[a].z != 16 && at_[a].z != 15) return false;
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k) {
      auto it = idx.find(at_[a].nb[k]);
      if (ringbond_[a][k] && it != idx.end()) adj[i].push_back(it->second), ++nbonds;
    }
  }
  nbonds /= 2;
  // no spiro junctions: the ring graph must stay connected when any one atom is removed
  for (int cut = 0; cut < n; ++cut) {
    std::vector<char> seen(n, 0);
    std::vector<int> st = {cut == 0 ? 1 : 0};
    seen[st[0]] = 1;
    int cnt = 1;
    while (!st.empty()) {
      const int a = st.back();
      st.pop_back();
      for (int b : adj[a])
        if (b != cut && !seen[b]) seen[b] = 1, ++cnt, st.push_back(b);
    }
    if (cnt != n - 1) return fail("spiro junction in a polycyclic ring system at " + sym(V[cut]));
  }
  // all simple cycles (each once: smallest atom first, second smaller than last)
  std::vector<std::vector<int>> cycles;
  {
    std::vector<int> path;
    std::vector<char> on(n, 0);
    long steps = 0;
    auto dfs = [&](auto&& self, int start, int a) -> void {
      if (++steps > 2000000 || cycles.size() > 20000) return;
      for (int b : adj[a]) {
        if (b == start && path.size() >= 3 && path[1] < path.back()) cycles.push_back(path);
        if (b <= start || on[b]) continue;
        on[b] = 1;
        path.push_back(b);
        self(self, start, b);
        path.pop_back();
        on[b] = 0;
      }
    };
    for (int s0 = 0; s0 < n; ++s0) {
      path = {s0};
      on[s0] = 1;
      dfs(dfs, s0, s0);
      on[s0] = 0;
    }
    if (steps > 2000000 || cycles.size() > 20000) return fail("ring system too complex for von Baeyer nomenclature");
  }
  std::size_t maxring = 0;
  for (const auto& c : cycles) maxring = std::max(maxring, c.size());
  struct Main {
    std::vector<int> ring, bridge;  // bridge atoms from ring[pu] to ring[pv]
    int pu = 0, pv = 0;
  };
  std::vector<Main> mains;
  int best_bridge = -1, best_split = 1 << 20;
  for (const auto& C : cycles) {
    if (C.size() != maxring) continue;
    const int m = static_cast<int>(C.size());
    std::vector<int> pos(n, -1);
    for (int i = 0; i < m; ++i) pos[C[i]] = i;
    auto consider = [&](int pu, int pv, const std::vector<int>& bridge) {
      if (pu > pv) return;
      if (pu == pv) return;
      const int inner = pv - pu - 1, outer = m - 2 - inner;
      const int len = static_cast<int>(bridge.size()), split = std::abs(inner - outer);
      if (len < best_bridge || (len == best_bridge && split > best_split)) return;
      if (len > best_bridge || split < best_split) {
        best_bridge = len;
        best_split = split;
        mains.clear();
      }
      mains.push_back({C, bridge, pu, pv});
    };
    for (int i = 0; i < m; ++i) {
      const int u = C[i];
      for (int b : adj[u]) {
        if (pos[b] >= 0) {  // chord: a bridge of no atoms
          const int d = std::abs(pos[b] - i);
          if (d != 1 && d != m - 1) consider(i, pos[b], {});
          continue;
        }
        // paths through atoms outside the ring back to it
        std::vector<int> path = {b};
        std::vector<char> on(n, 0);
        on[b] = 1;
        auto walk = [&](auto&& self, int a) -> void {
          for (int x : adj[a]) {
            if (pos[x] >= 0) {
              if (x != u || path.size() > 1) {
                if (x != u) consider(i, pos[x], path);
              }
              continue;
            }
            if (on[x]) continue;
            on[x] = 1;
            path.push_back(x);
            self(self, x);
            path.pop_back();
            on[x] = 0;
          }
        };
        walk(walk, b);
      }
    }
  }
  if (mains.empty()) return fail("no main bridge in a polycyclic ring system");

  struct Sec {  // a secondary bridge: dependent?, atoms, superscript locants
    int dependent, len, lo, hi;
  };
  struct Result {
    std::vector<int> o;
    std::string desc;
    std::vector<int> supers;
  };
  std::vector<Result> results;
  long branches = 0;
  auto bond_id = [&](int a, int b) { return std::min(a, b) * n + std::max(a, b); };
  // Secondary bridges. First the remaining structure is cut into bridges (independent ones before
  // dependent, larger before smaller); then the bridge atoms are numbered, each bridge in turn from
  // the one on the highest-numbered bridgehead down, starting at the end nearer that bridgehead.
  auto extend = [&](std::vector<int> o, std::set<int> used, std::string desc, int main_size) -> void {
    if (++branches > 200000) return;
    struct Bridge {
      std::vector<int> path;
      int u, v;
      bool dependent;
    };
    std::vector<Bridge> bridges;
    std::vector<char> assigned(n, 0), in_main(n, 0);
    for (int a : o) assigned[a] = in_main[a] = 1;
    for (;;) {
      std::vector<Bridge> cands;
      for (int u = 0; u < n; ++u) {
        if (!assigned[u]) continue;
        for (int b : adj[u]) {
          if (used.count(bond_id(u, b))) continue;
          if (assigned[b]) {
            if (u < b) cands.push_back({{}, u, b, !(in_main[u] && in_main[b])});
            continue;
          }
          std::vector<int> path = {b};
          std::vector<char> on(n, 0);
          on[b] = 1;
          auto walk = [&](auto&& w, int a) -> void {
            for (int x : adj[a]) {
              if (assigned[x]) {
                if (x != u && u < x && !used.count(bond_id(a, x))) cands.push_back({path, u, x, !(in_main[u] && in_main[x])});
                continue;
              }
              if (on[x]) continue;
              on[x] = 1;
              path.push_back(x);
              w(w, x);
              path.pop_back();
              on[x] = 0;
            }
          };
          walk(walk, b);
        }
      }
      if (cands.empty()) break;
      const Bridge* pick = nullptr;
      for (const auto& c : cands)
        if (!pick || std::make_pair(!c.dependent, c.path.size()) > std::make_pair(!pick->dependent, pick->path.size())) pick = &c;
      Bridge b = *pick;
      used.insert(bond_id(b.u, b.path.empty() ? b.v : b.path.front()));
      for (std::size_t i = 0; i + 1 < b.path.size(); ++i) used.insert(bond_id(b.path[i], b.path[i + 1]));
      if (!b.path.empty()) used.insert(bond_id(b.path.back(), b.v));
      for (int a : b.path) assigned[a] = 1;
      bridges.push_back(b);
    }
    if (static_cast<int>(used.size()) != nbonds) return;
    std::vector<int> where(n, -1);
    for (std::size_t i = 0; i < o.size(); ++i) where[o[i]] = static_cast<int>(i);
    std::vector<char> done(bridges.size(), 0);
    std::vector<Sec> secs;
    for (std::size_t round = 0; round < bridges.size(); ++round) {
      int best = -1;
      std::pair<int, int> top = {-1, -1};
      for (std::size_t i = 0; i < bridges.size(); ++i) {
        const Bridge& b = bridges[i];
        if (done[i] || where[b.u] < 0 || where[b.v] < 0) continue;
        // independent bridges are numbered before dependent ones
        bool waits = false;
        if (b.dependent)
          for (std::size_t k = 0; k < bridges.size(); ++k) waits |= !done[k] && !bridges[k].dependent;
        if (waits) continue;
        const std::pair<int, int> key = {std::max(where[b.u], where[b.v]), std::min(where[b.u], where[b.v])};
        if (key > top) top = key, best = static_cast<int>(i);
      }
      if (best < 0) return;
      Bridge& b = bridges[static_cast<std::size_t>(best)];
      done[static_cast<std::size_t>(best)] = 1;
      std::vector<int> p = b.path;
      if (where[b.v] > where[b.u]) std::reverse(p.begin(), p.end());
      for (int a : p) {
        where[a] = static_cast<int>(o.size());
        o.push_back(a);
      }
      secs.push_back({b.dependent ? 1 : 0, static_cast<int>(b.path.size()), top.second + 1, top.first + 1});
    }
    if (static_cast<int>(o.size()) != n) return;
    (void)main_size;
    // cited: independent before dependent, larger first, then lower superscripts
    std::sort(secs.begin(), secs.end(), [](const Sec& a, const Sec& b) {
      return std::make_tuple(a.dependent, -a.len, a.lo, a.hi) < std::make_tuple(b.dependent, -b.len, b.lo, b.hi);
    });
    std::vector<int> supers;
    for (const Sec& x : secs) {
      desc += "." + std::to_string(x.len) + "^{" + std::to_string(x.lo) + "," + std::to_string(x.hi) + "}";
      supers.push_back(x.lo);
      supers.push_back(x.hi);
    }
    results.push_back({o, desc, supers});
  };
  for (const Main& M : mains) {
    const int m = static_cast<int>(M.ring.size());
    std::vector<int> seg1(M.ring.begin() + M.pu + 1, M.ring.begin() + M.pv), seg2;
    for (int i = M.pv + 1; i < M.pu + m; ++i) seg2.push_back(M.ring[i % m]);
    const int hu = M.ring[M.pu], hv = M.ring[M.pv];
    // seg1 runs hu -> hv, seg2 runs hv -> hu, the bridge hu -> hv
    for (int flip = 0; flip < 2; ++flip)
      for (int which = 0; which < 2; ++which) {
        std::vector<int> A = which ? seg2 : seg1, B = which ? seg1 : seg2, C = M.bridge;
        int h1 = hu, h2 = hv;
        if (which) std::swap(h1, h2);  // seg2 runs hv -> hu: start at hv
        if (A.size() < B.size()) continue;
        if (which) std::reverse(C.begin(), C.end());
        if (flip) {  // start from the other main bridgehead
          std::reverse(A.begin(), A.end());
          std::reverse(B.begin(), B.end());
          std::reverse(C.begin(), C.end());
          std::swap(h1, h2);
        }
        std::vector<int> o = {h1};
        o.insert(o.end(), A.begin(), A.end());
        o.push_back(h2);
        o.insert(o.end(), B.begin(), B.end());
        o.insert(o.end(), C.begin(), C.end());
        std::set<int> used;
        for (int i = 0; i < m; ++i) used.insert(bond_id(M.ring[i], M.ring[(i + 1) % m]));
        if (C.empty()) used.insert(bond_id(h1, h2));
        else {
          used.insert(bond_id(h1, C.front()));
          for (std::size_t i = 0; i + 1 < C.size(); ++i) used.insert(bond_id(C[i], C[i + 1]));
          used.insert(bond_id(C.back(), h2));
        }
        const std::string desc = std::to_string(A.size()) + "." + std::to_string(B.size()) + "." + std::to_string(C.size());
        extend(o, used, desc, static_cast<int>(o.size()));
      }
  }
  if (branches > 200000) return fail("ring system too complex for von Baeyer nomenclature");
  if (results.empty()) return fail("ring system that could not be described by von Baeyer nomenclature");
  auto skey = [](std::vector<int> s) {
    std::sort(s.begin(), s.end());
    return s;
  };
  std::vector<int> best = skey(results[0].supers);
  for (const auto& r : results) best = std::min(best, skey(r.supers));
  static const char* cyc[] = {"", "", "bicyclo", "tricyclo", "tetracyclo", "pentacyclo", "hexacyclo", "heptacyclo", "octacyclo"};
  std::set<std::vector<int>> seen;
  for (const auto& r : results) {
    if (skey(r.supers) != best) continue;
    std::vector<int> o;
    for (int i : r.o) o.push_back(V[i]);
    if (!seen.insert(o).second) continue;
    R.orders.push_back(o);
    R.names.push_back(std::string(cyc[R.nrings]) + "[" + r.desc + "]");
  }
  R.kind = PKind::VonBaeyer;
  R.unsat = 0;
  for (int x : V)
    for (std::size_t k = 0; k < at_[x].nb.size(); ++k)
      if (ringbond_[x][k] && at_[x].nb[k] > x && at_[x].bo[k] > 1 && idx.count(at_[x].nb[k])) ++R.unsat;
  return true;
}

// Names, numberings and saturation of one component of a spiro system (P-24.5): a saturated
// monocycle, or a ring system with a template or fusion name.
bool Namer::spiro_component(const std::vector<int>& atoms, int nrings, int centre, SpiroPart& out) {
  out = SpiroPart{};
  const std::set<int> in(atoms.begin(), atoms.end());
  auto ring_bond = [&](int a, std::size_t k) { return ringbond_[a][k] && in.count(at_[a].nb[k]) > 0; };
  RingSys T;
  T.atoms = atoms;
  T.nrings = nrings;
  // retained and fusion templates whose atoms are exactly this component
  {
    const auto& qs = template_queries();
    SmartsMatchOptions mo;
    mo.unique = false;
    mo.max_matches = 0;
    for (std::size_t t = 0; t < qs.size() && T.orders.empty(); ++t) {
      if (!qs[t].valid() || templates()[t].locants.size() != atoms.size() || pattern_rings(templates()[t].pattern) != nrings) continue;
      for (const auto& mt : match_smarts(st_, sp_, qs[t], mo)) {
        std::vector<int> o;
        for (std::uint32_t g : mt) {
          auto it = local_.find(g);
          if (it == local_.end() || !in.count(it->second)) break;
          o.push_back(it->second);
        }
        if (o.size() != mt.size()) continue;
        T.templ = static_cast<int>(t);
        T.orders.push_back(o);
        T.names.push_back(templates()[t].name);
        T.locs.push_back(templates()[t].locants);
      }
    }
  }
  bool mancude = T.templ >= 0 && templates()[T.templ].mancude;
  if (T.orders.empty() && nrings == 1) {
    // saturated monocycle: cycloalkane or Hantzsch-Widman name
    for (int a : atoms)
      for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
        if (ring_bond(a, k) && at_[a].bo[k] != 1) return fail("unsaturated monocycle in a spiro system");
    std::vector<int> cyc = {atoms[0]};
    for (int prev = -1, cur = atoms[0]; cyc.size() < atoms.size();) {
      int next = -1;
      for (std::size_t k = 0; k < at_[cur].nb.size(); ++k)
        if (ring_bond(cur, k) && at_[cur].nb[k] != prev && at_[cur].nb[k] != atoms[0]) next = at_[cur].nb[k];
      if (next < 0) return false;
      cyc.push_back(next);
      prev = cur;
      cur = next;
    }
    bool hetero = false;
    for (int a : atoms) hetero |= at_[a].z != 6;
    if (!hetero) {
      const int m = static_cast<int>(cyc.size());
      for (int st = 0; st < m; ++st)
        for (int dir : {1, -1}) {
          std::vector<int> o;
          std::vector<std::string> l;
          for (int i = 0; i < m; ++i) o.push_back(cyc[((st + dir * i) % m + m) % m]), l.push_back(std::to_string(i + 1));
          T.orders.push_back(o);
          T.locs.push_back(l);
          T.names.push_back("cyclo" + alkane_stem(m) + "ane");
        }
    } else {
      bool dummy = false;
      const std::string saved = err_;
      if (!hantzsch_widman(T, dummy) || T.orders.empty()) return false;
      err_ = saved;
    }
  }
  if (T.orders.empty() && nrings >= 2) {
    const std::string saved = err_;
    const bool ok = nrings == 2 ? fusion_bicycle(T) : fusion_polycycle(T);
    err_ = saved;
    if (!ok || T.orders.empty()) return fail("spiro component without a fusion name");
    mancude = true;
  }
  if (T.orders.empty()) return false;
  // lowest locants to heteroatoms, then to the spiro atom (P-24.5.3)
  std::vector<std::vector<double>> keys;
  for (std::size_t i = 0; i < T.orders.size(); ++i) {
    std::vector<std::string> het, sp;
    for (std::size_t j = 0; j < T.orders[i].size(); ++j) {
      if (at_[T.orders[i][j]].z != 6) het.push_back(T.locs[i][j]);
      if (T.orders[i][j] == centre) sp.push_back(T.locs[i][j]);
    }
    auto k = sorted_keys(het);
    k.push_back(-1);
    const auto k2 = sorted_keys(sp);
    k.insert(k.end(), k2.begin(), k2.end());
    keys.push_back(k);
  }
  const auto best = *std::min_element(keys.begin(), keys.end());
  for (std::size_t i = 0; i < T.orders.size(); ++i)
    if (keys[i] == best) {
      out.orders.push_back(T.orders[i]);
      out.locs.push_back(T.locs[i]);
      out.names.push_back(T.names[i]);
    }
  out.mancude = mancude;
  out.atoms = atoms;
  return true;
}

// Spiro systems with at least one polycyclic component (P-24.5): "spiro[cyclohexane-1,3'-indole]".
bool Namer::spiro_fused(RingSys& R) {
  const std::set<int> in(R.atoms.begin(), R.atoms.end());
  auto rdeg = [&](int a) {
    int d = 0;
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k) d += ringbond_[a][k] && in.count(at_[a].nb[k]);
    return d;
  };
  // the spiro atom: removing it splits the ring system in two
  int centre = -1;
  std::vector<std::vector<int>> parts;
  for (int c : R.atoms) {
    if (rdeg(c) != 4) continue;
    std::vector<std::vector<int>> comps;
    std::set<int> seen = {c};
    for (int s0 : R.atoms) {
      if (seen.count(s0)) continue;
      std::vector<int> comp = {s0};
      seen.insert(s0);
      for (std::size_t h = 0; h < comp.size(); ++h)
        for (std::size_t k = 0; k < at_[comp[h]].nb.size(); ++k) {
          const int b = at_[comp[h]].nb[k];
          if (ringbond_[comp[h]][k] && in.count(b) && !seen.count(b)) seen.insert(b), comp.push_back(b);
        }
      comps.push_back(comp);
    }
    if (comps.size() != 2) continue;
    if (centre >= 0) return fail("ring system with several spiro atoms");
    centre = c;
    parts = comps;
  }
  if (centre < 0) return false;
  SpiroPart P[2];
  for (int i = 0; i < 2; ++i) {
    std::vector<int> atoms = parts[i];
    atoms.push_back(centre);
    std::sort(atoms.begin(), atoms.end());
    std::set<int> as(atoms.begin(), atoms.end());
    int nb = 0;
    for (int a : atoms)
      for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
        if (ringbond_[a][k] && at_[a].nb[k] > a && as.count(at_[a].nb[k])) ++nb;
    const int nrings = nb - static_cast<int>(atoms.size()) + 1;
    if (!spiro_component(atoms, nrings, centre, P[i])) return false;
  }
  // components cited in alphabetical order; the second one takes primed locants
  if (alpha_key(P[1].names[0]) < alpha_key(P[0].names[0])) std::swap(P[0], P[1]);
  for (std::size_t i = 0; i < P[0].orders.size(); ++i)
    for (std::size_t j = 0; j < P[1].orders.size(); ++j) {
      std::vector<int> o = P[0].orders[i];
      std::vector<std::string> L = P[0].locs[i];
      std::string l0, l1;
      for (std::size_t k = 0; k < o.size(); ++k)
        if (o[k] == centre) l0 = L[k];
      for (std::size_t k = 0; k < P[1].orders[j].size(); ++k) {
        const int a = P[1].orders[j][k];
        if (a == centre) {
          l1 = P[1].locs[j][k] + "'";
          continue;
        }
        o.push_back(a);
        L.push_back(P[1].locs[j][k] + "'");
      }
      R.orders.push_back(o);
      R.locs.push_back(L);
      R.names.push_back("spiro[" + bracket_locants(P[0].names[i]) + "-" + l0 + "," + l1 + "-" + bracket_locants(P[1].names[j]) + "]");
    }
  R.kind = PKind::Template;
  R.templ = -1;
  R.mancude = P[0].mancude || P[1].mancude;
  R.spiro_mancude = true;
  for (int a : R.atoms)
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
      if (ringbond_[a][k] && at_[a].nb[k] > a && at_[a].bo[k] == 2) ++R.unsat;
  // saturation of the mancude components only; the spiro atom is saturated by the spiro union
  for (const auto& p : P) {
    if (!p.mancude) continue;
    const std::set<int> as(p.atoms.begin(), p.atoms.end());
    for (int a : p.atoms) {
      if (a == centre) continue;
      int rd = 0;
      bool dbl = false;
      for (std::size_t k = 0; k < at_[a].nb.size(); ++k)
        if (ringbond_[a][k] && as.count(at_[a].nb[k])) ++rd, dbl |= at_[a].bo[k] == 2;
      const bool inherent = ((at_[a].z == 8 || at_[a].z == 16 || at_[a].z == 15) && rd == 2) || (at_[a].z == 7 && rd == 3);
      if (inherent) continue;
      R.matchable.push_back(a);
      if (!dbl) R.X.push_back(a);
    }
  }
  // indicated hydrogen where a component needs it ("3',4'-dihydro-2'H-spiro[imidazolidine-4,1'-naphthalene]");
  // hydro() checks that a valid position exists for every numbering
  R.ih = static_cast<int>(R.matchable.size() % 2);
  if (!R.ih && !perfect_matching(R.matchable)) R.ih = 2;
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
  std::vector<int> picks;
  if (R.ih) {
    std::vector<int> xs = R.X;
    std::sort(xs.begin(), xs.end(), [&](int a, int b) { return locant_key(L[a]) < locant_key(L[b]); });
    auto feasible = [&](const std::vector<int>& out) {
      std::vector<int> rest;
      for (int a : R.matchable)
        if (std::find(out.begin(), out.end(), a) == out.end()) rest.push_back(a);
      return perfect_matching(rest);
    };
    if (R.ih == 1) {
      for (int x : xs)
        if (feasible({x})) {
          picks = {x};
          break;
        }
    } else {
      for (std::size_t i = 0; i < xs.size() && picks.empty(); ++i)
        for (std::size_t j = i + 1; j < xs.size() && picks.empty(); ++j)
          if (feasible({xs[i], xs[j]})) picks = {xs[i], xs[j]};
    }
    if (picks.empty()) {
      ind = "?";  // no indicated-hydrogen position gives a valid mancude parent
      return;
    }
    for (std::size_t i = 0; i < picks.size(); ++i) ind += (i ? "H," : "") + L[picks[i]];  // "1H,5" + "H-"
  }
  for (int x : R.X)
    if (std::find(picks.begin(), picks.end(), x) == picks.end()) hyd.push_back(L[x]);
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
  for (int a : o)  // cationic centres expressed by "-ium" next to another suffix (P-73.1.1.2)
    if (ringCat_[a] && !(sp.mode == Mode::Molecule && sp.principal == Cls::RingCation)) suf.push_back(L[a]);
  {
    std::string ind;
    std::vector<std::string> hyd;
    hydro(c, L, ind, hyd);
    if (ind == "?") valid = false;
    for (std::size_t p0 = 0; !ind.empty() && p0 <= ind.size();) {  // one key per indicated hydrogen
      std::size_t comma = ind.find(',', p0);
      if (comma == std::string::npos) comma = ind.size();
      s.ind.push_back(locant_key(ind.substr(p0, comma - p0)));
      p0 = comma + 1;
    }
    if (hyd.size() % 2) valid = false;  // hydro prefixes come in pairs; an odd set has no mancude reference
    un.insert(un.end(), hyd.begin(), hyd.end());
  }
  if (c.kind == PKind::Template && rs_[c.rsys].het_numbering) {
    // heteroatoms as a set, then in the order O, S, N (P-31.1.4.2.2)
    std::vector<std::string> all;
    for (int x : o)
      if (at_[x].z != 6) all.push_back(L[x]);
    s.het = sorted_keys(all);
    for (int z : {8, 16, 7, 15}) {
      std::vector<std::string> l;
      for (int x : o)
        if (at_[x].z == z) l.push_back(L[x]);
      const auto k = sorted_keys(l);
      s.het.insert(s.het.end(), k.begin(), k.end());
    }
  }
  if (sp.mode == Mode::Yl || sp.mode == Mode::Ylidene || sp.mode == Mode::Oyl) suf.push_back(L[sp.root]);
  if (sp.mode == Mode::Carbonyl) suf.push_back(L[sp.link]);
  if (c.kind == PKind::VonBaeyer)
    for (const auto& [loc, bo2] : ring_multiple_bonds(c, o)) {
      un.push_back(loc);
      if (bo2 == 2) ene.push_back(loc);
    }
  if (c.kind == PKind::Chain || c.kind == PKind::Carbocycle || c.kind == PKind::VonBaeyer) {
    const std::size_t m = o.size();
    const std::size_t lim = c.kind == PKind::VonBaeyer ? 0 : c.kind == PKind::Carbocycle ? m : m - 1;
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
      (void)counted;
      std::vector<std::string> het, by_kind;
      for (int z : {8, 16, 7, 15})
        for (int x : o)
          if (at_[x].z == z) het.push_back(L[x]);
      s.het = sorted_keys(het);
      for (int z : {8, 16, 7, 15}) {
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

// Lambda convention for ring sulfur oxides: "1,2-benzothiazole" -> "1λ6,2-benzothiazole",
// "thiolane" -> "1λ6-thiolane".
std::string Namer::lambda_name(std::string name, const Cand& c, const std::vector<int>& o) const {
  const auto L = locants(c, o);
  for (int a : o) {
    if (!lambdaS_[a]) continue;
    const std::string tag = L[a] + lambda_tag(a);
    std::size_t end = 0;
    // the leading locant list, including earlier lambda tags ("1,3,5,2λ5,4-")
    while (end < name.size() && ((name[end] >= '0' && name[end] <= '9') || name[end] == ',' || static_cast<unsigned char>(name[end]) >= 0x80 ||
                                 (name[end] >= 'a' && name[end] <= 'z' && end > 0 && name[end - 1] >= '0' && name[end - 1] <= '9')))
      ++end;
    bool done = false;
    if (end > 0 && end < name.size() && name[end] == '-') {
      std::string list = name.substr(0, end), out;
      std::size_t pos = 0;
      while (pos <= list.size()) {
        std::size_t comma = list.find(',', pos);
        if (comma == std::string::npos) comma = list.size();
        std::string tok = list.substr(pos, comma - pos);
        if (tok == L[a] && !done) tok = tag, done = true;
        out += (out.empty() ? "" : ",") + tok;
        pos = comma + 1;
      }
      if (done) name = out + name.substr(end);
    }
    if (!done) name = tag + "-" + name;
  }
  return name;
}

// Multiple bonds of a von Baeyer or spiro system: locant of the lower-numbered atom, with the other
// atom in parentheses when the two are not numbered consecutively ("2(7)", P-31.1.4.2.4).
std::vector<std::pair<std::string, int>> Namer::ring_multiple_bonds(const Cand& c, const std::vector<int>& o) const {
  std::vector<std::pair<std::string, int>> out;
  std::map<int, int> pos;
  for (std::size_t i = 0; i < o.size(); ++i) pos[o[i]] = static_cast<int>(i) + 1;
  (void)c;
  for (int a : o)
    for (std::size_t k = 0; k < at_[a].nb.size(); ++k) {
      const int b = at_[a].nb[k];
      auto it = pos.find(b);
      if (!ringbond_[a][k] || at_[a].bo[k] < 2 || it == pos.end() || pos[a] > it->second) continue;
      const int i = pos[a], j = it->second;
      out.push_back({j == i + 1 ? std::to_string(i) : std::to_string(i) + "(" + std::to_string(j) + ")", at_[a].bo[k]});
    }
  return out;
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
    append_piece(s, lambda_name(rs_[c.rsys].names[static_cast<std::size_t>(order_index(c, o))], c, o));
    return s;
  }
  const std::size_t m = o.size();
  const std::size_t lim = c.kind == PKind::VonBaeyer ? 0 : c.kind == PKind::Carbocycle ? m : m - 1;
  std::vector<std::string> enes, ynes;
  for (std::size_t i = 0; i < lim; ++i) {
    const int b = order(o[i], o[(i + 1) % m]);
    if (b == 2) enes.push_back(std::to_string(i + 1));
    if (b == 3) ynes.push_back(std::to_string(i + 1));
  }
  if (c.kind == PKind::VonBaeyer)
    for (const auto& [loc, bo2] : ring_multiple_bonds(c, o)) (bo2 == 2 ? enes : ynes).push_back(loc);
  std::string stem = (c.kind == PKind::Carbocycle ? "cyclo" : "") + alkane_stem(static_cast<int>(m));
  if (c.kind == PKind::VonBaeyer || (c.kind == PKind::Carbocycle && rs_[c.rsys].nhetero > 0)) {
    const auto L = locants(c, o);
    std::string rep_prefixes;
    for (const auto& [z, word] : std::vector<std::pair<int, std::string>>{{8, "oxa"}, {16, "thia"}, {7, "aza"}, {15, "phospha"}}) {
      std::vector<std::string> l;
      for (int x : o)
        if (at_[x].z == z) l.push_back(L[x]);
      if (l.empty()) continue;
      for (auto& x : l)
        for (int a2 : o)
          if (lambdaS_[a2] && L[a2] == x) x += lambda_tag(a2);
      rep_prefixes += (rep_prefixes.empty() ? "" : "-") + join_locants(l) + "-" + multiplier(static_cast<int>(l.size()), false) + word;
    }
    const RingSys& R = rs_[c.rsys];
    if (c.kind == PKind::VonBaeyer && R.spiro)
      stem = rep_prefixes + "spiro[" + std::to_string(R.vb[0]) + "." + std::to_string(R.vb[1]) + "]" + stem;
    else if (c.kind == PKind::VonBaeyer && !R.names.empty())
      stem = rep_prefixes + R.names[static_cast<std::size_t>(order_index(c, o))] + stem;
    else if (c.kind == PKind::VonBaeyer)
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
  std::string hyd = hydride(c, o, omit_unsat, saturated);
  {
    std::vector<std::string> ium;
    for (int a : o) {
      if (ringCat_[a] && !(mol && sp.principal == Cls::RingCation)) ium.push_back(L[a]);
      if (ringAn_[a] && !(mol && sp.principal == Cls::RingAnion)) return fail("ring anion outside the parent"), std::nullopt;
    }
    if (!ium.empty()) {
      if (c.kind == PKind::Chain) return fail("cationic chain atom"), std::nullopt;
      hyd = attach(hyd, join_locants(ium), multiplier(static_cast<int>(ium.size()), false) + "ium");
      saturated = false;  // "...-1-ium" stays whole before the next suffix
    }
  }
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
        case Cls::SulfonoHydrazide: word = "sulfonohydrazide"; break;
        case Cls::Ester: word = c.inchain ? "oate" : "carboxylate"; break;
        case Cls::Amide: word = c.inchain ? "amide" : "carboxamide"; break;
        case Cls::Nitrile: word = c.inchain ? "nitrile" : "carbonitrile"; break;
        case Cls::Aldehyde: word = c.inchain ? "al" : "carbaldehyde"; break;
        case Cls::Ketone: word = "one"; break;
        case Cls::Thione: word = "thione"; break;
        case Cls::Thioamide: word = c.inchain ? "thioamide" : "carbothioamide"; break;
        case Cls::Hydrazide: word = c.inchain ? "ohydrazide" : "carbohydrazide"; break;
        case Cls::Amidine: word = c.inchain ? "imidamide" : "carboximidamide"; break;
        case Cls::Alcohol: word = "ol"; break;
        case Cls::Thiol: word = "thiol"; break;
        case Cls::Peroxol: word = "peroxol"; break;
        case Cls::Diazonium: word = "diazonium"; break;
        case Cls::Amine: word = "amine"; break;
        case Cls::Carboxylate: word = c.inchain ? "oate" : "carboxylate"; break;
        case Cls::Sulfonate: word = "sulfonate"; break;
        case Cls::SulfEster: word = "sulfonate"; break;
        case Cls::SulfinicAcid: word = "sulfinic acid"; break;
        case Cls::ThioEster: word = c.inchain ? "thioate" : "carbothioate"; break;
        case Cls::Olate: word = "olate"; break;
        case Cls::Ammonium: word = "aminium"; break;
        case Cls::RingCation: word = "ium"; break;
        case Cls::RingAnion: word = "ide"; break;
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
      if ((k == Cls::Acid || k == Cls::Carboxylate || k == Cls::Anhydride || k == Cls::Nitrile || k == Cls::Carbonate || k == Cls::Carbamate || k == Cls::Urea || k == Cls::Thiourea ||
           k == Cls::Semicarbazide || k == Cls::ThioSemicarbazide || k == Cls::AcylSub || k == Cls::CarbonicAcid || k == Cls::CarbonicAnion) && !PA(a))
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
    if (acyl_[a] != Cls::None && static_cast<int>(acyl_[a]) < 500) senior(acyl_[a]);
    if (is_OH(a) && at_[at_[a].nb[0]].z == 6 && acyl_[at_[a].nb[0]] != Cls::Acid && !carbX_[a]) senior(Cls::Alcohol);
    if (is_SH(a) && !carbX_[a]) senior(Cls::Thiol);
    if (pseudo_[a] == Hydroperoxy) senior(Cls::Peroxol);
    if (ionC_[a] != Cls::None) senior(ionC_[a]);
    if (amineN_[a]) senior(Cls::Amine);
    if (cationN_[a]) senior(Cls::Ammonium);
    if (ringCat_[a]) senior(Cls::RingCation);
    if (ringAn_[a]) senior(Cls::RingAnion);
    if (anionO_[a] && at_[at_[a].nb[0]].z == 6 && acyl_[at_[a].nb[0]] != Cls::Carboxylate && !carbX_[a]) senior(Cls::Olate);
    if (imineC_[a]) senior(Cls::Imine);
    if (hzPartner_[a] >= 0) senior(Cls::Hydrazine);
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
  int silicon = -1, nsilicon = 0;
  for (int a = 0; a < n(); ++a)
    if (at_[a].z == 14) silicon = a, ++nsilicon;
  int azo = -1, nazo = 0;
  for (int a = 0; a < n(); ++a)
    if (azoPartner_[a] > a) azo = a, ++nazo;
  int phosphane = -1, nphosphane = 0;
  for (int a = 0; a < n(); ++a)
    if (acyl_[a] == Cls::Phosphane) phosphane = a, ++nphosphane;
  if (principal == Cls::Phosphane && !(nphosphane == 1 && nsilicon == 0 && nazo == 0)) principal = Cls::None;
  if (principal == Cls::Phosphane) {
    // phosphane as the senior parent hydride without a suffix: "triphenylphosphane"
    std::vector<int> xs(at_[phosphane].nb.begin(), at_[phosphane].nb.end());
    auto R = ligands(xs, phosphane, "phosphane");
    if (!R) return std::nullopt;
    result = R->text;
    cited = R->stereo;
  } else if (principal == Cls::None && nsilicon == 0 && nazo == 1) {
    // diazene as the senior parent hydride when there is no suffix (P-68.3.1.3): "diphenyldiazene"
    const int M = azoPartner_[azo];
    std::vector<Named> subs;
    std::vector<Key> extra;
    for (int X : {azo, M})
      for (int y : at_[X].nb) {
        if (y == azo || y == M) continue;
        auto R = sub(y, X, 1);
        if (!R) return std::nullopt;
        extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
        subs.push_back(*R);
      }
    std::map<std::string, std::pair<int, bool>> g;
    for (const auto& x : subs) ++g[x.text].first, g[x.text].second = x.compound;
    std::vector<std::pair<std::string, std::pair<int, bool>>> gv(g.begin(), g.end());
    std::sort(gv.begin(), gv.end(), [](const auto& a2, const auto& b2) { return alpha_key(a2.first) < alpha_key(b2.first); });
    std::string name;
    for (std::size_t i = 0; i < gv.size(); ++i) {
      const auto& [t, kc] = gv[i];
      std::string piece = kc.second ? (kc.first > 1 ? multiplier(kc.first, true) : "") + enclose(t) : multiplier(kc.first, false) + t;
      if (i > 0 && !is_enclosed(piece)) piece = enclose(piece);
      name += piece;
    }
    name += "diazene";
    if (auto it = stereo_.find(bond_key(azo, M)); it != stereo_.end()) {
      name = "(" + std::string(to_string(it->second)) + ")-" + name;
      extra.push_back(it->first);
    }
    result = name;
    cited = extra;
  } else if (principal == Cls::None && nsilicon == 1 && !at_[silicon].nb.empty()) {
    // silane as the senior parent hydride when there is no suffix (P-44.1.2.1): "trimethyl(phenyl)silane"
    std::vector<int> xs(at_[silicon].nb.begin(), at_[silicon].nb.end());
    auto R = ligands(xs, silicon, "silane");
    if (!R) return std::nullopt;
    result = R->text;
    cited = R->stereo;
  } else if (principal == Cls::PhosMono || principal == Cls::Phosphonic) {
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
  } else if (principal == Cls::Semicarbazide || principal == Cls::ThioSemicarbazide) {
    // semicarbazides and semicarbazones as hydrazine-1-carboxamides (P-68.3.1.2.2):
    // "2-(propan-2-ylidene)hydrazine-1-carboxamide", "N-phenylhydrazine-1-carbothioamide"
    int nsemi = 0;
    for (int C = 0; C < n(); ++C) nsemi += acyl_[C] == principal;
    for (int C = 0; C < n() && nsemi == 1; ++C) {
      if (acyl_[C] != principal) continue;
      int Na = -1, Nh = -1;
      for (int x : at_[C].nb)
        if (at_[x].z == 7) (hzPartner_[x] >= 0 ? Nh : Na) = x;
      if (Na < 0 || Nh < 0) break;
      const int Nz = hzPartner_[Nh];
      std::vector<std::pair<std::string, Named>> e;
      std::vector<Key> extra;
      bool ok = true;
      for (auto [X, loc] : {std::pair{Na, "N"}, std::pair{Nh, "1"}, std::pair{Nz, "2"}})
        for (std::size_t k = 0; k < at_[X].nb.size() && ok; ++k) {
          const int y = at_[X].nb[k];
          if (y == C || y == Nh || y == Nz) continue;
          auto R = sub(y, X, at_[X].bo[k]);
          if (!R) { ok = false; break; }
          extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
          e.push_back({loc, *R});
        }
      if (!ok) break;
      const std::string word = principal == Cls::Semicarbazide ? "carboxamide" : "carbothioamide";
      finish(Built{}, e.empty() ? "hydrazine" + word : format_prefixes(e, false) + "hydrazine-1-" + word, extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "semicarbazide could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::Hydrazine) {
    // hydrazine as parent (P-68.3.1.2): "phenylhydrazine", "1,1-dimethylhydrazine", "(propan-2-ylidene)hydrazine"
    for (int N = 0; N < n(); ++N) {
      const int M = hzPartner_[N];
      if (M < N) continue;
      std::vector<Named> sN, sM;
      std::vector<Key> extra;
      bool ok = true;
      for (int X : {N, M})
        for (std::size_t k = 0; k < at_[X].nb.size() && ok; ++k) {
          const int y = at_[X].nb[k];
          if (y == N || y == M) continue;
          auto R = sub(y, X, at_[X].bo[k]);
          if (!R) { ok = false; break; }
          extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
          (X == N ? sN : sM).push_back(*R);
        }
      if (!ok) continue;
      auto least = [](const std::vector<Named>& v) {
        std::string m = "~";
        for (const auto& x : v) m = std::min(m, alpha_key(x.text));
        return m;
      };
      if (sM.size() > sN.size() || (sM.size() == sN.size() && least(sM) < least(sN))) std::swap(sN, sM);
      std::vector<std::pair<std::string, Named>> e;
      for (const auto& x : sN) e.push_back({"1", x});
      for (const auto& x : sM) e.push_back({"2", x});
      finish(Built{}, format_prefixes(e, e.size() == 1) + "hydrazine", extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "hydrazine could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::Guanidine) {
    // "N-methyl-N'-phenylguanidine": amino nitrogens N and N', imino nitrogen N'' (P-66.4.1.2.1)
    for (int C = 0; C < n(); ++C) {
      if (acyl_[C] != Cls::Guanidine) continue;
      std::vector<int> amino;
      int imino = -1;
      for (int x : at_[C].nb)
        if (amidN_[x]) imino = x;
        else if (at_[x].z == 7) amino.push_back(x);
      std::vector<Key> extra;
      bool ok = imino >= 0 && amino.size() == 2;
      auto subs_of = [&](int N, std::vector<Named>& out) {
        for (int x : at_[N].nb) {
          if (x == C) continue;
          auto R = sub(x, N, 1);
          if (!R) return false;
          extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
          out.push_back(*R);
        }
        return true;
      };
      std::vector<Named> s0, s1, si;
      if (ok) ok = subs_of(amino[0], s0) && subs_of(amino[1], s1) && subs_of(imino, si);
      if (!ok) continue;
      if (s1.size() > s0.size()) std::swap(s0, s1);
      std::vector<std::pair<std::string, Named>> e;
      for (const auto& x : s0) e.push_back({"N", x});
      for (const auto& x : s1) e.push_back({"N'", x});
      for (const auto& x : si) e.push_back({"N''", x});
      finish(Built{}, format_prefixes(e, false) + "guanidine", extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "guanidine could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::Carbonate || principal == Cls::Carbamate || principal == Cls::Urea || principal == Cls::Thiourea ||
             principal == Cls::CarbonicAcid || principal == Cls::CarbonicAnion) {
    // functional class names of carbonic acid derivatives (P-65.2.1): "dimethyl carbonate",
    // "ethyl N-methylcarbamate", "N-methyl-N'-phenylurea"
    for (int C = 0; C < n(); ++C) {
      if (acyl_[C] != principal) continue;
      if (principal != Cls::Urea && principal != Cls::Thiourea) {
        if (auto r = carbonic_name(C)) finish(Built{}, r->first, r->second);
        continue;
      }
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
    }
    if (!result) {
      if (err_.empty()) err_ = "carbonic acid derivative could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::Anhydride) {
    // carboxylic anhydrides (P-65.7.1): "acetic anhydride", "acetic propanoic anhydride"
    auto acid_word = [](std::string t) {
      if (t == "acetyl") return std::string("acetic");
      if (t == "formyl") return std::string("formic");
      if (t == "benzoyl") return std::string("benzoic");
      if (ends_with(t, "oyl")) return t.substr(0, t.size() - 3) + "oic";
      if (ends_with(t, "carbonyl")) return t.substr(0, t.size() - 8) + "carboxylic";
      return std::string();
    };
    for (int a = 0; a < n(); ++a) {
      if (acyl_[a] != Cls::Anhydride) continue;
      const int O = ester_O(a);
      const int b2 = at_[O].nb[0] == a ? at_[O].nb[1] : at_[O].nb[0];
      if (b2 < a) continue;
      auto A1 = acyl_name(a, O), A2 = acyl_name(b2, O);
      if (!A1 || !A2) continue;
      std::string w1 = acid_word(A1->text), w2 = acid_word(A2->text);
      if (w1.empty() || w2.empty()) continue;
      if (alpha_key(w2) < alpha_key(w1)) std::swap(w1, w2);
      std::vector<Key> extra = A1->stereo;
      extra.insert(extra.end(), A2->stereo.begin(), A2->stereo.end());
      finish(Built{}, (w1 == w2 ? w1 : w1 + " " + w2) + " anhydride", extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "anhydride could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::PolyP) {
    // monoesters of di- and triphosphoric acid (P-67.2.1): "R trihydrogen diphosphate"
    std::vector<int> esterO;
    for (int a = 0; a < n(); ++a)
      if (at_[a].z == 8 && deg(a) == 2) {
        const int x = at_[a].nb[0], y = at_[a].nb[1];
        if ((acyl_[x] == Cls::PolyP && at_[y].z == 6) || (acyl_[y] == Cls::PolyP && at_[x].z == 6)) esterO.push_back(a);
      }
    // the phosphorus chain: P atoms joined by bridging oxygens
    std::vector<int> chainP;
    int hydrogens = 0;
    for (int a = 0; a < n(); ++a)
      if (acyl_[a] == Cls::PolyP) {
        chainP.push_back(a);
        for (int x : at_[a].nb) hydrogens += is_OH(x);
      }
    bool linear = esterO.size() == 1 && chainP.size() >= 2 && chainP.size() <= 3;
    if (linear) {  // the ester sits on a terminal phosphorus and the chain is unbranched
      for (int Pp : chainP) {
        int bridges = 0;
        for (int x : at_[Pp].nb)
          if (at_[x].z == 8 && deg(x) == 2)
            for (int y : at_[x].nb) bridges += y != Pp && at_[y].z == 15;
        if (bridges > 2) linear = false;
        if (bridges == 2 && std::count(at_[Pp].nb.begin(), at_[Pp].nb.end(), esterO[0])) linear = false;
      }
    }
    if (esterO.empty()) {  // diphosphoric and triphosphoric acid
      if (chainP.size() > 3 || hydrogens != static_cast<int>(chainP.size()) + 2) {
        err_ = "unsupported polyphosphoric acid";
        return std::nullopt;
      }
      result = chainP.size() == 2 ? "diphosphoric acid" : "triphosphoric acid";
      linear = false;
    } else if (!linear) {
      // diesters and branched esters: one phosphorus as the parent phosphate, the others as
      // phosphoryl groups ("[hydroxy(methoxy)phosphoryl] methyl hydrogen phosphate", as in
      // "acetyl dihydrogen phosphate"); the shortest such name
      for (int P : chainP) {
        std::map<std::string, int> g;
        std::set<std::string> wrapped;  // phosphoryl groups, already enclosed here
        std::vector<Key> extra;
        int hyd = 0;
        bool ok = true;
        for (std::size_t k = 0; k < at_[P].nb.size() && ok; ++k) {
          const int x = at_[P].nb[k];
          if (at_[P].bo[k] != 1) continue;
          if (is_OH(x)) { ++hyd; continue; }
          if (at_[x].z != 8 || deg(x) != 2) { ok = false; break; }
          const int r = at_[x].nb[0] == P ? at_[x].nb[1] : at_[x].nb[0];
          auto R = at_[r].z == 15 ? phosphoryl_name(r, x) : sub(r, x, 1);
          if (!R) { ok = false; break; }
          extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
          if (at_[r].z == 15) wrapped.insert(enclose(R->text));
          ++g[at_[r].z == 15 ? enclose(R->text) : organyl(*R)];
        }
        if (!ok) continue;
        std::vector<std::pair<std::string, int>> gv(g.begin(), g.end());
        std::sort(gv.begin(), gv.end(), [](const auto& a2, const auto& b2) { return alpha_key(a2.first) < alpha_key(b2.first); });
        std::string name;
        for (const auto& [t, k] : gv) {
          const bool comp = t[0] == '(' || t[0] == '[' || t[0] == '{' || has_digit(t);
          name += (name.empty() ? "" : " ") +
                  (k > 1 ? (comp ? multiplier(k, true) + (wrapped.count(t) ? t : enclose(t)) : multiplier(k, false) + t) : t);
        }
        if (hyd) name += " " + multiplier(hyd, false) + "hydrogen";
        finish(Built{}, name + " phosphate", extra);
      }
      if (!result) {
        if (err_.empty()) err_ = "polyphosphate ester could not be named";
        return std::nullopt;
      }
    }
    if (linear) {
    const int O = esterO[0];
    const int r = at_[at_[O].nb[0]].z == 6 ? at_[O].nb[0] : at_[O].nb[1];
    auto R = sub(r, O, 1);
    if (!R) return std::nullopt;
    result = organyl(*R) + " " + multiplier(hydrogens, false) + "hydrogen " + (chainP.size() == 2 ? "diphosphate" : "triphosphate");
    cited = R->stereo;
    }
  } else if (principal == Cls::AcylHalide || principal == Cls::SulfonylHalide) {
    // acid halides by functional class names (P-65.5.1): "benzoyl chloride", "methanesulfonyl chloride"
    for (int a = 0; a < n(); ++a) {
      if (acyl_[a] != principal) continue;
      int X = -1;
      for (int x : at_[a].nb)
        if (at_[x].z == 9 || at_[x].z == 17 || at_[x].z == 35 || at_[x].z == 53) X = x;
      auto R = principal == Cls::AcylHalide ? acyl_name(a, X) : sulfonyl_name(a, X);
      if (!R) continue;
      const int z = at_[X].z;
      finish(Built{}, R->text + (z == 9 ? " fluoride" : z == 17 ? " chloride" : z == 35 ? " bromide" : " iodide"), R->stereo);
    }
    if (!result) {
      if (err_.empty()) err_ = "acid halide could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::Sulfonium || principal == Cls::Phosphonium || principal == Cls::Borate || principal == Cls::Azanide ||
             principal == Cls::Nitrate) {
    // "trimethylsulfanium", "tetramethylphosphanium", "tetrafluoroborate", "(benzenesulfonyl)(phenyl)azanide",
    // "nitric acid", "nitrate" (P-73.1, P-72.2)
    int centre = -1, count = 0;
    for (int a = 0; a < n(); ++a)
      if (ionC_[a] == principal) centre = a, ++count;
    if (count != 1) {
      err_ = "several ionic centres of one kind";
      return std::nullopt;
    }
    if (principal == Cls::Nitrate) {
      int oh = 0;
      for (int x : at_[centre].nb) oh += is_OH(x);
      if (n() != 4) {
        err_ = "nitrate ester or nitric acid derivative";
        return std::nullopt;
      }
      result = oh ? "nitric acid" : "nitrate";
    } else {
      std::vector<int> xs(at_[centre].nb.begin(), at_[centre].nb.end());
      if (principal == Cls::Borate) {  // "(4-aminobutyl)trihydroxyborate": simple ligands unenclosed
        std::map<std::string, std::pair<int, bool>> g;
        for (int x : xs) {
          auto R = sub(x, centre, 1);
          if (!R) return std::nullopt;
          cited.insert(cited.end(), R->stereo.begin(), R->stereo.end());
          auto& e = g[R->text];
          ++e.first;
          e.second = R->compound;
        }
        std::vector<std::pair<std::string, std::pair<int, bool>>> v(g.begin(), g.end());
        std::sort(v.begin(), v.end(), [](const auto& x, const auto& y) { return alpha_key(x.first) < alpha_key(y.first); });
        std::string t;
        for (const auto& [txt, e] : v)
          t += e.second ? (e.first > 1 ? multiplier(e.first, true) : "") + enclose(txt) : multiplier(e.first, false) + txt;
        result = t + "borate";
      } else {
        const char* word = principal == Cls::Sulfonium ? "sulfanium" : principal == Cls::Phosphonium ? "phosphanium" : "azanide";
        auto R = ligands(xs, centre, word);
        if (!R) return std::nullopt;
        result = R->text;
        cited = R->stereo;
      }
    }
  } else if (principal == Cls::SulfuricDiamide || principal == Cls::SulfamateEster) {
    // "N,N-dimethylsulfuric diamide" (P-66.1.1.4.3), "methyl N,N-dimethylsulfamate" (P-67.1.3.2)
    for (int S = 0; S < n(); ++S) {
      if (acyl_[S] != principal) continue;
      std::vector<int> ns;
      int org = -1;
      for (int x : at_[S].nb) {
        if (at_[x].z == 7) ns.push_back(x);
        if (at_[x].z == 8 && deg(x) == 2) org = x;
      }
      std::vector<std::vector<Named>> nsubs;
      std::vector<Key> extra;
      bool ok = true;
      for (int N : ns) {
        nsubs.emplace_back();
        for (int x : at_[N].nb) {
          if (x == S) continue;
          auto R = sub(x, N, 1);
          if (!R) { ok = false; break; }
          extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
          nsubs.back().push_back(*R);
        }
      }
      if (!ok) continue;
      std::vector<std::pair<std::string, Named>> e;
      if (principal == Cls::SulfuricDiamide) {
        int first = 0;
        auto least = [](const std::vector<Named>& v) {
          std::string m = "~";
          for (const auto& x : v) m = std::min(m, alpha_key(x.text));
          return m;
        };
        if (nsubs[1].size() > nsubs[0].size() || (nsubs[1].size() == nsubs[0].size() && !nsubs[0].empty() && least(nsubs[1]) < least(nsubs[0])))
          first = 1;
        for (const auto& nm : nsubs[first]) e.push_back({"N", nm});
        for (const auto& nm : nsubs[1 - first]) e.push_back({"N'", nm});
        finish(Built{}, format_prefixes(e, false) + "sulfuric diamide", extra);
      } else {
        const int r = at_[org].nb[0] == S ? at_[org].nb[1] : at_[org].nb[0];
        auto R = sub(r, org, 1);
        if (!R) continue;
        extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
        for (const auto& nm : nsubs[0]) e.push_back({"N", nm});
        finish(Built{}, organyl(*R) + " " + format_prefixes(e, false) + "sulfamate", extra);
      }
    }
    if (!result) {
      if (err_.empty()) err_ = "sulfuric acid amide could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::Isothiourea || principal == Cls::IsothioAcid) {
    // isothioureas as carbamimidothioates (P-66.4.1.2): "2-aminoethyl carbamimidothioate",
    // "ethyl N'-phenylcarbamimidothioate", "N'-phenylcarbamimidothioic acid"; N amino, N' imino
    for (int C = 0; C < n(); ++C) {
      if (acyl_[C] != principal) continue;
      int Nam = -1, Nim = -1, Sx = -1;
      for (std::size_t k = 0; k < at_[C].nb.size(); ++k) {
        const int x = at_[C].nb[k];
        if (at_[x].z == 7) (at_[C].bo[k] == 2 ? Nim : Nam) = x;
        if (at_[x].z == 16) Sx = x;
      }
      std::vector<std::pair<std::string, Named>> e;
      std::vector<Key> extra;
      bool ok = true;
      for (auto [X, loc] : {std::pair{Nam, "N"}, std::pair{Nim, "N'"}})
        for (int y : at_[X].nb) {
          if (y == C) continue;
          auto R = sub(y, X, 1);
          if (!R) { ok = false; break; }
          extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
          e.push_back({loc, *R});
        }
      if (!ok) continue;
      std::string org;
      if (principal == Cls::Isothiourea) {
        const int r = at_[Sx].nb[0] == C ? at_[Sx].nb[1] : at_[Sx].nb[0];
        auto R = sub(r, Sx, 1);
        if (!R) continue;
        extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
        org = organyl(*R) + " ";
      }
      std::string word = format_prefixes(e, false) + "carbamimidothio" + (principal == Cls::Isothiourea ? "ate" : "ic acid");
      if (auto it = stereo_.find(bond_key(C, Nim)); it != stereo_.end()) {
        word = "(" + std::string(to_string(it->second)) + ")-" + word;
        extra.push_back(it->first);
      }
      finish(Built{}, org + word, extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "isothiourea could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::SulfateAnion || principal == Cls::SulfamicAcid || principal == Cls::SulfamateAnion ||
             principal == Cls::SulfuricInorg) {
    // "ethyl sulfate", "N-methylsulfamic acid", "N,N-dimethylsulfamate", "sulfuric acid", "hydrogen sulfate"
    for (int S = 0; S < n() && !result; ++S) {
      if (acyl_[S] != principal) continue;
      int nOH = 0, nOm = 0;
      for (int x : at_[S].nb) nOH += is_OH(x), nOm += is_Ominus(x);
      if (principal == Cls::SulfuricInorg) {
        finish(Built{}, nOH == 2 ? "sulfuric acid" : nOH == 1 ? "hydrogen sulfate" : "sulfate", {});
        continue;
      }
      if (principal == Cls::SulfateAnion) {
        for (int x : at_[S].nb) {
          if (at_[x].z != 8 || deg(x) != 2) continue;
          const int r = at_[x].nb[0] == S ? at_[x].nb[1] : at_[x].nb[0];
          if (auto R = sub(r, x, 1)) finish(Built{}, organyl(*R) + " sulfate", R->stereo);
        }
        continue;
      }
      int N = -1;
      for (int x : at_[S].nb)
        if (at_[x].z == 7) N = x;
      std::vector<std::pair<std::string, Named>> e;
      std::vector<Key> extra;
      bool ok = true;
      for (int x : at_[N].nb) {
        if (x == S) continue;
        auto R = sub(x, N, 1);
        if (!R) { ok = false; break; }
        extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
        e.push_back({"N", *R});
      }
      if (ok) finish(Built{}, format_prefixes(e, false) + (principal == Cls::SulfamicAcid ? "sulfamic acid" : "sulfamate"), extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "sulfuric acid derivative could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::SulfateMono) {
    // "R hydrogen sulfate" (P-67.1.3.2)
    for (int S = 0; S < n(); ++S) {
      if (acyl_[S] != Cls::SulfateMono) continue;
      for (int x : at_[S].nb) {
        if (at_[x].z != 8 || deg(x) != 2) continue;
        const int r = at_[x].nb[0] == S ? at_[x].nb[1] : at_[x].nb[0];
        auto R = sub(r, x, 1);
        if (R) finish(Built{}, organyl(*R) + " hydrogen sulfate", R->stereo);
      }
    }
    if (!result) {
      if (err_.empty()) err_ = "sulfate ester could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::PhosThio) {
    // "O,O-diethyl O-(4-nitrophenyl) phosphorothioate", "O,O-dimethyl S-alkyl phosphorodithioate" (P-67.1.3.2)
    for (int P = 0; P < n(); ++P) {
      if (acyl_[P] != Cls::PhosThio) continue;
      std::map<std::string, std::vector<std::string>> by;  // organyl -> letters
      std::map<std::string, bool> comp;
      std::vector<Key> extra;
      int nS = 0;
      bool ok = true;
      for (std::size_t k = 0; k < at_[P].nb.size(); ++k) {
        const int x = at_[P].nb[k];
        if (at_[P].bo[k] == 2) {
          nS += at_[x].z == 16;
          continue;
        }
        const int r = at_[x].nb[0] == P ? at_[x].nb[1] : at_[x].nb[0];
        auto R = sub(r, x, 1);
        if (!R) { ok = false; break; }
        nS += at_[x].z == 16;
        by[R->text].push_back(at_[x].z == 8 ? "O" : "S");
        comp[R->text] = R->compound || R->text[0] == '(';
        extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
      }
      if (!ok) continue;
      std::vector<std::pair<std::string, std::vector<std::string>>> v;
      for (auto& [t, ls] : by) {  // one group per organyl and letter
        std::map<std::string, int> cnt;
        for (const auto& l : ls) ++cnt[l];
        for (const auto& [l, c] : cnt) v.push_back({t, std::vector<std::string>(static_cast<std::size_t>(c), l)});
      }
      std::sort(v.begin(), v.end(), [](const auto& a2, const auto& b2) { return alpha_key(a2.first) < alpha_key(b2.first); });
      std::string name;
      for (const auto& [t, ls] : v) {
        std::string letters;
        for (const auto& l : ls) letters += (letters.empty() ? "" : ",") + l;
        const int c = static_cast<int>(ls.size());
        const std::string body = c > 1 ? (comp[t] ? multiplier(c, true) + enclose(t) : multiplier(c, false) + t) : (comp[t] ? enclose(t) : t);
        name += (name.empty() ? "" : " ") + letters + "-" + body;
      }
      name += nS >= 2 ? " phosphorodithioate" : " phosphorothioate";
      finish(Built{}, name, extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "thiophosphate could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::PhosOxide || principal == Cls::PhosphinicAcid) {
    // "triphenyl-λ5-phosphanone" (P-74.2), "methyl(phenyl)phosphinic acid" (P-67.1.1)
    for (int P = 0; P < n(); ++P) {
      if (acyl_[P] != principal) continue;
      std::vector<int> xs;
      for (std::size_t k = 0; k < at_[P].nb.size(); ++k)
        if (at_[P].bo[k] == 1 && at_[at_[P].nb[k]].z == 6) xs.push_back(at_[P].nb[k]);
      auto R = ligands(xs, P, "");
      if (!R) continue;
      finish(Built{}, principal == Cls::PhosOxide ? R->text + "-λ5-phosphanone" : R->text + "phosphinic acid", R->stereo);
    }
    if (!result) {
      if (err_.empty()) err_ = "phosphorus compound could not be named";
      return std::nullopt;
    }
  } else if (principal == Cls::PhosEster || principal == Cls::PhosAnion || principal == Cls::PhosAcidInorg) {
    // esters of phosphoric and phosphonic acids (P-67.1.3.2): "triethyl phosphate", "dimethyl methylphosphonate",
    // "diethyl hydrogen phosphate"
    for (int P = 0; P < n(); ++P) {
      if (acyl_[P] != principal) continue;
      if (principal == Cls::PhosAcidInorg) {
        finish(Built{}, "phosphoric acid", {});
        continue;
      }
      std::vector<std::string> org;
      std::vector<Key> extra;
      int hydrogens = 0;
      std::string cpart;
      bool ok = true;
      for (std::size_t k = 0; k < at_[P].nb.size() && ok; ++k) {
        const int x = at_[P].nb[k];
        if (at_[P].bo[k] != 1) continue;
        if (is_OH(x)) ++hydrogens;
        else if (is_Ominus(x)) continue;
        else if (at_[x].z == 8) {
          const int r = at_[x].nb[0] == P ? at_[x].nb[1] : at_[x].nb[0];
          auto R = sub(r, x, 1);
          if (!R) ok = false;
          else org.push_back(organyl(*R)), extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
        } else if (at_[x].z == 6) {
          auto R = sub(x, P, 1);
          if (!R) ok = false;
          else cpart = (R->compound ? enclose(R->text) : R->text), extra.insert(extra.end(), R->stereo.begin(), R->stereo.end());
        }
      }
      if (!ok) continue;
      std::map<std::string, int> g;
      for (const auto& o : org) ++g[o];
      std::vector<std::pair<std::string, int>> gv(g.begin(), g.end());
      std::sort(gv.begin(), gv.end(), [](const auto& a2, const auto& b2) { return alpha_key(a2.first) < alpha_key(b2.first); });
      std::string name;
      for (const auto& [t, k] : gv) {
        const bool comp = t[0] == '(' || t[0] == '[' || has_digit(t);
        name += (name.empty() ? "" : " ") + (k > 1 ? (comp ? multiplier(k, true) + enclose(t) : multiplier(k, false) + t) : t);
      }
      if (hydrogens) name += (name.empty() ? "" : " ") + multiplier(hydrogens, false) + "hydrogen";
      name += (name.empty() ? "" : " ") + (cpart.empty() ? std::string("phosphate") : cpart + "phosphonate");
      finish(Built{}, name, extra);
    }
    if (!result) {
      if (err_.empty()) err_ = "phosphorus ester could not be named";
      return std::nullopt;
    }
  } else if (principal != Cls::Ester && principal != Cls::ThioEster && principal != Cls::SulfEster) {
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
      if (acyl_[a] == principal) E.push_back(a);
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
      sp.principal = principal;
      auto b = build(sp);
      if (!b || b->ester_o.empty()) continue;
      std::map<std::string, std::pair<int, bool>> alk;
      std::map<std::string, std::vector<std::string>> alk_locs;
      std::vector<Key> extra;
      bool ok = true;
      for (const auto& [o, loc] : b->ester_o) {
        const int x0 = at_[o].nb[0], x1 = at_[o].nb[1];
        const int r = (acyl_[x0] == principal && sp.pok[x0]) ? x1 : x0;
        auto R = sub(r, o, 1);
        if (!R) {
          ok = false;
          break;
        }
        const std::string org = principal == Cls::ThioEster ? "S-" + (R->compound ? enclose(R->text) : R->text)  // S-methyl ethanethioate
                                                             : organyl(*R);
        auto& g = alk[org];
        ++g.first;
        g.second = R->compound;
        alk_locs[org].push_back(loc);
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
    if (neutrals.size() > 1 && cations.empty() && anions.empty() && adducts.empty()) {
      // addition compound (P-74.2): names joined by em dashes with the ratio, fewer components first
      std::map<std::string, int> cnt;
      for (const auto& x : neutrals) ++cnt[x.name];
      if (waters) cnt["water"] += waters;
      std::vector<std::pair<std::string, int>> v(cnt.begin(), cnt.end());
      std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
        return a.second != b.second ? a.second < b.second : alpha_key(a.first) < alpha_key(b.first);
      });
      std::string names, ratio;
      for (const auto& [nm, k] : v) {
        names += (names.empty() ? "" : "—") + nm;
        ratio += (ratio.empty() ? "" : "/") + std::to_string(k);
      }
      return names + " (" + ratio + ")";
    }
    if (neutrals.size() != 1 || !cations.empty() || !anions.empty() || (adducts.empty() && waters == 0)) {
      err = "mixture of neutral and ionic components";
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
