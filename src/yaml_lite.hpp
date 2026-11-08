#pragma once
// Minimal YAML reader for mechanism files (internal; not a general YAML library).
//
// Supported: block mappings and sequences (indentation based), "- key: v" items,
// flow collections {...} and [...] (may span lines), plain / single / double
// quoted scalars, '#' comments, '|' and '>' block scalars (kept as text), '---'.
// Not supported: anchors/aliases, tags, multi-document streams, complex keys.

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace omnichem::yaml {

struct Node {
  enum class Kind { Null, Scalar, Map, Seq };
  Kind kind = Kind::Null;
  std::string scalar;
  std::vector<std::pair<std::string, Node>> map;
  std::vector<Node> seq;
  std::size_t line = 0;  // 1-based source line

  bool is_null() const noexcept { return kind == Kind::Null; }
  bool is_scalar() const noexcept { return kind == Kind::Scalar; }
  bool is_map() const noexcept { return kind == Kind::Map; }
  bool is_seq() const noexcept { return kind == Kind::Seq; }
  const Node* get(std::string_view key) const noexcept {
    if (kind != Kind::Map) return nullptr;
    for (const auto& [k, v] : map)
      if (k == key) return &v;
    return nullptr;
  }
};

struct ParseError {
  std::string message;
  std::size_t line = 0;
};

// Returns false and fills `error` on malformed input.
bool parse(std::string_view text, Node& out, ParseError& error);

}  // namespace omnichem::yaml
