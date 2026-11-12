#include "yaml_lite.hpp"

#include <cctype>

namespace chemlarp::yaml {
namespace {

struct Logical {
  int indent;
  std::string content;
  std::size_t line;
};

std::string trim(std::string_view s) {
  std::size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return std::string(s.substr(a, b - a));
}

// Removes a trailing comment ('#' at start or after whitespace, outside quotes).
std::string strip_comment(const std::string& s) {
  char quote = 0;
  for (std::size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (quote) {
      if (c == quote) quote = 0;
      continue;
    }
    if (c == '"' || c == '\'') quote = c;
    else if (c == '#' && (i == 0 || std::isspace(static_cast<unsigned char>(s[i - 1])))) return s.substr(0, i);
  }
  return s;
}

int bracket_balance(const std::string& s) {
  int depth = 0;
  char quote = 0;
  for (char c : s) {
    if (quote) {
      if (c == quote) quote = 0;
      continue;
    }
    if (c == '"' || c == '\'') quote = c;
    else if (c == '[' || c == '{') ++depth;
    else if (c == ']' || c == '}') --depth;
  }
  return depth;
}

// Position of the ':' separating key and value (followed by space or end), or npos.
std::size_t key_colon(const std::string& s) {
  char quote = 0;
  int depth = 0;
  for (std::size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (quote) {
      if (c == quote) quote = 0;
      continue;
    }
    if (c == '"' || c == '\'') quote = c;
    else if (c == '[' || c == '{') ++depth;
    else if (c == ']' || c == '}') --depth;
    else if (c == ':' && depth == 0 && (i + 1 == s.size() || s[i + 1] == ' ')) return i;
  }
  return std::string::npos;
}

std::string unquote(const std::string& s) {
  if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\''))) {
    std::string out;
    const bool dbl = s.front() == '"';
    for (std::size_t i = 1; i + 1 < s.size(); ++i) {
      if (dbl && s[i] == '\\' && i + 2 < s.size()) {
        const char e = s[++i];
        out += e == 'n' ? '\n' : e == 't' ? '\t' : e;
      } else if (!dbl && s[i] == '\'' && s[i + 1] == '\'') {
        out += '\'';
        ++i;
      } else {
        out += s[i];
      }
    }
    return out;
  }
  return s;
}

class Parser {
 public:
  Parser(std::vector<Logical> lines, ParseError& err) : l_(std::move(lines)), err_(err) {}

  bool run(Node& out) {
    if (l_.empty()) return true;
    i_ = 0;
    out = block(l_[0].indent);
    if (ok_ && i_ < l_.size()) fail("unexpected indentation", l_[i_].line);
    return ok_;
  }

 private:
  void fail(const std::string& m, std::size_t line) {
    if (!ok_) return;
    ok_ = false;
    err_.message = m;
    err_.line = line;
  }
  static bool is_item(const std::string& c) { return c == "-" || (c.size() >= 2 && c[0] == '-' && c[1] == ' '); }

  Node block(int indent) {
    if (i_ < l_.size() && is_item(l_[i_].content)) return sequence(indent);
    return mapping(indent);
  }

  Node sequence(int indent) {
    Node n;
    n.kind = Node::Kind::Seq;
    n.line = l_[i_].line;
    while (ok_ && i_ < l_.size() && l_[i_].indent == indent && is_item(l_[i_].content)) {
      const std::string& c = l_[i_].content;
      std::size_t off = 1;
      while (off < c.size() && c[off] == ' ') ++off;
      std::string rest = c.substr(off);
      if (rest.empty()) {
        ++i_;
        if (i_ < l_.size() && l_[i_].indent > indent) n.seq.push_back(block(l_[i_].indent));
        else n.seq.emplace_back();
      } else if ((rest[0] != '[' && rest[0] != '{' && rest[0] != '"' && rest[0] != '\'') &&
                 key_colon(rest) != std::string::npos) {
        // "- key: value" opens a mapping whose other keys sit at indent + off.
        l_[i_].indent = indent + static_cast<int>(off);
        l_[i_].content = rest;
        n.seq.push_back(mapping(l_[i_].indent));
      } else if (is_item(rest)) {
        l_[i_].indent = indent + static_cast<int>(off);
        l_[i_].content = rest;
        n.seq.push_back(sequence(l_[i_].indent));
      } else {
        const std::size_t line = l_[i_].line;
        ++i_;
        n.seq.push_back(value(rest, line));
      }
    }
    return n;
  }

  Node mapping(int indent) {
    Node n;
    n.kind = Node::Kind::Map;
    n.line = i_ < l_.size() ? l_[i_].line : 0;
    while (ok_ && i_ < l_.size() && l_[i_].indent == indent && !is_item(l_[i_].content)) {
      const std::string c = l_[i_].content;
      const std::size_t line = l_[i_].line;
      const std::size_t colon = key_colon(c);
      if (colon == std::string::npos) {
        fail("expected 'key: value'", line);
        return n;
      }
      std::string key = unquote(trim(c.substr(0, colon)));
      std::string val = trim(c.substr(colon + 1));
      ++i_;
      Node child;
      if (val.empty()) {
        if (i_ < l_.size() && (l_[i_].indent > indent || (l_[i_].indent == indent && is_item(l_[i_].content))))
          child = block(l_[i_].indent);
        child.line = child.line ? child.line : line;
      } else {
        child = value(val, line);
      }
      n.map.emplace_back(std::move(key), std::move(child));
    }
    if (ok_ && i_ < l_.size() && l_[i_].indent > indent) fail("unexpected indentation", l_[i_].line);
    return n;
  }

  Node value(const std::string& text, std::size_t line) {
    if (!text.empty() && (text[0] == '[' || text[0] == '{')) {
      std::size_t pos = 0;
      Node n = flow(text, pos, line);
      skip_ws(text, pos);
      if (ok_ && pos != text.size()) fail("trailing characters after flow collection", line);
      return n;
    }
    Node n;
    n.kind = Node::Kind::Scalar;
    n.scalar = unquote(text);
    n.line = line;
    if (n.scalar == "~" || n.scalar == "null") n.kind = Node::Kind::Null;
    return n;
  }

  static void skip_ws(const std::string& s, std::size_t& p) {
    while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p]))) ++p;
  }

  Node flow_scalar(const std::string& s, std::size_t& p, std::size_t line) {
    skip_ws(s, p);
    Node n;
    n.kind = Node::Kind::Scalar;
    n.line = line;
    if (p < s.size() && (s[p] == '"' || s[p] == '\'')) {
      const char q = s[p];
      std::size_t e = p + 1;
      while (e < s.size() && s[e] != q) e += (s[e] == '\\' && q == '"') ? 2 : 1;
      n.scalar = unquote(s.substr(p, e - p + 1));
      p = e + 1;
      return n;
    }
    const std::size_t start = p;
    while (p < s.size() && s[p] != ',' && s[p] != ']' && s[p] != '}' &&
           !(s[p] == ':' && (p + 1 == s.size() || s[p + 1] == ' ')))
      ++p;
    n.scalar = trim(s.substr(start, p - start));
    return n;
  }

  Node flow(const std::string& s, std::size_t& p, std::size_t line) {
    skip_ws(s, p);
    Node n;
    n.line = line;
    if (p >= s.size()) return n;
    if (s[p] == '[') {
      n.kind = Node::Kind::Seq;
      ++p;
      skip_ws(s, p);
      if (p < s.size() && s[p] == ']') { ++p; return n; }
      while (ok_ && p < s.size()) {
        skip_ws(s, p);
        n.seq.push_back((s[p] == '[' || s[p] == '{') ? flow(s, p, line) : flow_scalar(s, p, line));
        skip_ws(s, p);
        if (p < s.size() && s[p] == ',') { ++p; continue; }
        if (p < s.size() && s[p] == ']') { ++p; return n; }
        fail("expected ',' or ']' in flow sequence", line);
      }
      fail("unterminated '['", line);
      return n;
    }
    if (s[p] == '{') {
      n.kind = Node::Kind::Map;
      ++p;
      skip_ws(s, p);
      if (p < s.size() && s[p] == '}') { ++p; return n; }
      while (ok_ && p < s.size()) {
        Node key = flow_scalar(s, p, line);
        skip_ws(s, p);
        Node val;
        if (p < s.size() && s[p] == ':') {
          ++p;
          skip_ws(s, p);
          val = (p < s.size() && (s[p] == '[' || s[p] == '{')) ? flow(s, p, line) : flow_scalar(s, p, line);
        }
        n.map.emplace_back(key.scalar, std::move(val));
        skip_ws(s, p);
        if (p < s.size() && s[p] == ',') { ++p; continue; }
        if (p < s.size() && s[p] == '}') { ++p; return n; }
        fail("expected ',' or '}' in flow mapping", line);
      }
      fail("unterminated '{'", line);
      return n;
    }
    return flow_scalar(s, p, line);
  }

  std::vector<Logical> l_;
  ParseError& err_;
  std::size_t i_ = 0;
  bool ok_ = true;
};

}  // namespace

bool parse(std::string_view text, Node& out, ParseError& error) {
  // Physical lines -> logical lines (comments removed, block scalars folded,
  // multi-line flow collections joined).
  std::vector<std::pair<std::string, std::size_t>> raw;
  std::size_t start = 0, no = 0;
  while (start <= text.size()) {
    std::size_t e = text.find('\n', start);
    if (e == std::string_view::npos) e = text.size();
    std::string line(text.substr(start, e - start));
    if (!line.empty() && line.back() == '\r') line.pop_back();
    raw.emplace_back(std::move(line), ++no);
    if (e == text.size()) break;
    start = e + 1;
  }
  auto indent_of = [](const std::string& s) {
    int n = 0;
    while (n < static_cast<int>(s.size()) && s[static_cast<std::size_t>(n)] == ' ') ++n;
    return n;
  };
  std::vector<Logical> lines;
  for (std::size_t k = 0; k < raw.size(); ++k) {
    const std::string& r = raw[k].first;
    if (r.find('\t') != std::string::npos && trim(r).size() && r.find_first_not_of(" \t") < r.find('\t')) {
    } else if (r.find('\t') != std::string::npos && r.find_first_not_of(' ') != std::string::npos &&
               r[r.find_first_not_of(' ')] == '\t') {
      error = {"tab characters are not allowed for indentation", raw[k].second};
      return false;
    }
    std::string c = trim(strip_comment(r));
    if (c.empty() || c == "---" || c == "...") continue;
    const int ind = indent_of(r);
    // Block scalar: "key: |" / "key: >-" etc. Keep the following more-indented lines as text.
    const std::size_t colon = key_colon(c);
    if (colon != std::string::npos) {
      const std::string v = trim(c.substr(colon + 1));
      if (!v.empty() && (v[0] == '|' || v[0] == '>') && v.find_first_not_of("|>+-0123456789") == std::string::npos) {
        std::string body;
        while (k + 1 < raw.size() && (trim(raw[k + 1].first).empty() || indent_of(raw[k + 1].first) > ind)) {
          ++k;
          body += trim(raw[k].first) + (v[0] == '|' ? "\n" : " ");
        }
        std::string escaped = "\"";
        for (char ch : trim(body)) escaped += (ch == '"' || ch == '\\') ? std::string("\\") + ch : std::string(1, ch);
        lines.push_back({ind, c.substr(0, colon) + ": " + escaped + "\"", raw[k].second});
        continue;
      }
    }
    Logical lg{ind, c, raw[k].second};
    while (bracket_balance(lg.content) > 0 && k + 1 < raw.size()) {
      ++k;
      lg.content += " " + trim(strip_comment(raw[k].first));
    }
    if (bracket_balance(lg.content) != 0) {
      error = {"unbalanced brackets", lg.line};
      return false;
    }
    lines.push_back(std::move(lg));
  }
  Parser p(std::move(lines), error);
  return p.run(out);
}

}  // namespace chemlarp::yaml
