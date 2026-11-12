#include <array>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "chemlarp/inchi.hpp"

namespace chemlarp {

std::array<std::uint8_t, 32> sha256(std::string_view data) {
  static constexpr std::uint32_t k[64] = {
      0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
      0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
      0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
      0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
      0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
      0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
      0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
      0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
  std::uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  std::vector<std::uint8_t> m(data.begin(), data.end());
  const std::uint64_t bits = static_cast<std::uint64_t>(data.size()) * 8u;
  m.push_back(0x80);
  while (m.size() % 64 != 56) m.push_back(0);
  for (int i = 7; i >= 0; --i) m.push_back(static_cast<std::uint8_t>(bits >> (8 * i)));
  auto rotr = [](std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); };
  for (std::size_t off = 0; off < m.size(); off += 64) {
    std::uint32_t w[64];
    for (std::size_t i = 0; i < 16; ++i)
      w[i] = (std::uint32_t{m[off + 4 * i]} << 24) | (std::uint32_t{m[off + 4 * i + 1]} << 16) |
             (std::uint32_t{m[off + 4 * i + 2]} << 8) | std::uint32_t{m[off + 4 * i + 3]};
    for (std::size_t i = 16; i < 64; ++i) {
      const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    std::uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (std::size_t i = 0; i < 64; ++i) {
      const std::uint32_t t1 = hh + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
      const std::uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
      hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
  }
  std::array<std::uint8_t, 32> out{};
  for (std::size_t i = 0; i < 8; ++i)
    for (std::size_t j = 0; j < 4; ++j) out[4 * i + j] = static_cast<std::uint8_t>(h[i] >> (24 - 8 * j));
  return out;
}

namespace {

// 14-bit value -> three letters. The 16384 triplets are AAA..ZZZ without those
// starting with 'E' and without TAA..TTV.
void put_triplet(std::string& s, std::uint32_t v) {
  std::uint32_t first = v / 676, rest = v % 676;
  // letters before 'T' in the reduced alphabet: A-D, F-S = 18 full blocks
  if (v >= 18u * 676u) {
    const std::uint32_t w = v - 18u * 676u;  // inside or after the truncated T block (160 entries: TTW..TZZ)
    if (w < 160u) {
      const std::uint32_t t = w + 516u;  // offset inside the T block
      s += 'T';
      s += static_cast<char>('A' + t / 26);
      s += static_cast<char>('A' + t % 26);
      return;
    }
    first = 19u + (w - 160u) / 676u;
    rest = (w - 160u) % 676u;
  }
  static constexpr char kFirst[] = "ABCDFGHIJKLMNOPQRSTUVWXYZ";  // no 'E'
  s += kFirst[first];
  s += static_cast<char>('A' + rest / 26);
  s += static_cast<char>('A' + rest % 26);
}

void put_doublet(std::string& s, std::uint32_t v) {  // 9 bits, plain AA..ZZ
  s += static_cast<char>('A' + v / 26);
  s += static_cast<char>('A' + v % 26);
}

}  // namespace

std::string inchi_key(std::string_view inchi) {
  constexpr std::string_view prefix = "InChI=1S/";
  if (inchi.substr(0, prefix.size()) != prefix) return {};
  inchi.remove_prefix(prefix.size());
  // major: formula, /c, /h, /q.  minor: every other layer except /p.
  std::string major, minor;
  int protons = 0;
  bool first = true;
  while (!inchi.empty() || first) {
    const std::size_t cut = inchi.find('/');
    const std::string_view layer = inchi.substr(0, cut);
    inchi = cut == std::string_view::npos ? std::string_view{} : inchi.substr(cut + 1);
    if (first) {
      major.append(layer);
      first = false;
      if (!layer.empty() && layer[0] == 'p') {  // "InChI=1S/p+1": no formula
        major.clear();
        protons = std::atoi(std::string(layer.substr(1)).c_str());
      }
    } else if (!layer.empty() && (layer[0] == 'c' || layer[0] == 'h' || layer[0] == 'q')) {
      major += '/';
      major.append(layer);
    } else if (!layer.empty() && layer[0] == 'p') {
      protons = std::atoi(std::string(layer.substr(1)).c_str());
    } else {
      minor += '/';
      minor.append(layer);
    }
    if (cut == std::string_view::npos) break;
  }
  if (minor.size() < 255) minor += minor;
  const auto a = sha256(major), b = sha256(minor);
  std::string key;
  put_triplet(key, a[0] | (std::uint32_t{a[1]} & 0x3fu) << 8);
  put_triplet(key, ((a[1] & 0xc0u) | std::uint32_t{a[2]} << 8 | (std::uint32_t{a[3]} & 0x0fu) << 16) >> 6);
  put_triplet(key, ((a[3] & 0xf0u) | std::uint32_t{a[4]} << 8 | (std::uint32_t{a[5]} & 0x03u) << 16) >> 4);
  put_triplet(key, ((a[5] & 0xfcu) | std::uint32_t{a[6]} << 8) >> 2);
  put_doublet(key, a[7] | (std::uint32_t{a[8]} & 0x01u) << 8);
  key += '-';
  put_triplet(key, b[0] | (std::uint32_t{b[1]} & 0x3fu) << 8);
  put_triplet(key, ((b[1] & 0xc0u) | std::uint32_t{b[2]} << 8 | (std::uint32_t{b[3]} & 0x0fu) << 16) >> 6);
  put_doublet(key, ((b[3] & 0xf0u) | (std::uint32_t{b[4]} & 0x1fu) << 8) >> 4);
  key += "SA-";
  if (protons == 0) key += 'N';
  else if (protons > 0) key += protons <= 12 ? "OPQRSTUVWXYZ"[protons - 1] : 'A';
  else key += -protons <= 12 ? "MLKJIHGFEDCB"[-protons - 1] : 'A';
  return key;
}

}  // namespace chemlarp
