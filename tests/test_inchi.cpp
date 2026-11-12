// Standard InChI / InChIKey against reference strings.
//   inchi_curated.inc  hand-picked structures (reference: NCI CACTUS, InChI 1.0x)
//   inchi_pubchem.inc  500 random PubChem compounds that were never used to
//                      tune the generator (reference: PubChem's InChI)
#include <cstddef>
#include <string>

#include "chemlarp/inchi.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

struct Row {
  const char* smiles;
  const char* inchi;
  const char* key;
};

const Row kCurated[] = {
#include "data/inchi_curated.inc"
};
const Row kPubchem[] = {
#include "data/inchi_pubchem.inc"
};

template <std::size_t N>
std::size_t matches(const Row (&rows)[N], std::string* first_miss = nullptr) {
  std::size_t ok = 0;
  for (const Row& r : rows) {
    Diagnostics d;
    const InchiResult got = inchi_from_smiles(r.smiles, d);
    if (got.ok && got.inchi == r.inchi && got.key == r.key) ++ok;
    else if (first_miss && first_miss->empty()) *first_miss = std::string(r.smiles) + " -> " + got.inchi;
  }
  return ok;
}

std::string hex(const std::array<std::uint8_t, 32>& h) {
  static const char* d = "0123456789abcdef";
  std::string s;
  for (auto b : h) { s += d[b >> 4]; s += d[b & 15]; }
  return s;
}

}  // namespace

TEST(sha256_known_vectors) {
  CHECK_EQ(hex(sha256("")), std::string("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"));
  CHECK_EQ(hex(sha256("abc")), std::string("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"));
}

TEST(inchi_key_of_string) {
  CHECK_EQ(inchi_key("InChI=1S/CH4/h1H4"), std::string("VNWKTOKETHGBQD-UHFFFAOYSA-N"));
  CHECK_EQ(inchi_key("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"), std::string("QNAYBMKLOCPYGJ-REOHCLBHSA-N"));
  CHECK(inchi_key("CCO").empty());
}

TEST(inchi_independent_of_input_spelling) {
  Diagnostics d;
  const std::string a = inchi_from_smiles("OC(=O)[C@@H](N)Cc1ccccc1", d).inchi;
  CHECK_EQ(inchi_from_smiles("c1ccc(C[C@H](N)C(O)=O)cc1", d).inchi, a);
  CHECK_EQ(inchi_from_smiles("N[C@@H](Cc1ccccc1)C(=O)O", d).inchi, a);
  CHECK_EQ(inchi_from_smiles("[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]", d).inchi, a);  // zwitterion
  CHECK(inchi_from_smiles("N[C@H](Cc1ccccc1)C(=O)O", d).inchi != a);
  CHECK(!inchi_from_smiles("C(", d).ok);
}

TEST(inchi_curated_reference_table) {
  // One known miss: C[C@H]([2H])O needs the isotopic stereo layer (not implemented).
  std::string miss;
  const std::size_t ok = matches(kCurated, &miss);
  const std::size_t n = sizeof(kCurated) / sizeof(kCurated[0]);
  if (ok + 1 < n) throw cltest::Failure(std::to_string(ok) + "/" + std::to_string(n) + " match; first miss " + miss);
}

TEST(inchi_pubchem_sample) {
  // 492 of these 500 matched when the table was frozen; this guards against regressions.
  const std::size_t ok = matches(kPubchem);
  if (ok < 492) throw cltest::Failure("only " + std::to_string(ok) + "/500 PubChem InChIs match");
}
