#include <memory>
#include <sstream>

#include "slmm_array_tree.hpp"

namespace slmm {
namespace tree {

namespace {

struct Node {
  typedef std::shared_ptr<Node> Ptr;
  std::vector<Ptr> kids;
  std::vector<Int> items;
  static Node::Ptr create () { return std::make_shared<Node>(); }
};

void read_tree (std::stringstream& ss, const Node::Ptr& node) {
  const auto get = [&] () {
    SIQK_THROW_IF(ss.eof(), "End of string; bad tree description.");
    return ss.get();
  };
  const auto eat_until_lparen = [&] () {
    while ( ! ss.eof() && ss.peek() != '(') get();
  };
  const auto eat_whitespace = [&] () {
    while ( ! ss.eof() && ss.peek() == ' ') get();
  };
  const auto read_int = [&] () -> int {
    int i;
    ss >> i;
    return i;
  };
  const auto read_item = [&] (Int& item) -> bool {
    eat_whitespace();
    const auto peek = ss.peek();
    if (peek == ')') return false;
    SIQK_THROW_IF(peek == '(', "parse::read_item: Expected integer.");
    item = read_int();
    return true;
  };
  eat_until_lparen();
  get();
  eat_whitespace();
  if (ss.peek() == '(') {
    for (;;) {
      const auto kid = Node::create();
      node->kids.push_back(kid);
      read_tree(ss, kid);
      eat_whitespace();
      const auto peek = ss.peek();
      if (peek == ')') break;
    }
  } else {
    Int item;
    while (read_item(item)) node->items.push_back(item);
  }
  get();
}

void write (std::stringstream& ss, const Node::Ptr& node) {
  ss << "(";
  if (node->kids.empty()) {
    bool first = true;
    for (const Int e : node->items) {
      if ( ! first) ss << " ";
      first = false;
      ss << e;
    }
  } else {
    for (size_t i = 0; i < node->kids.size(); ++i) {
      write(ss, node->kids[i]);
      if (i < node->kids.size()-1) ss << " ";
    }
  }
  ss << ")";
}

std::string tree_to_str (const Node::Ptr& node) {
  std::stringstream ss;
  write(ss, node);
  return ss.str();
}

Node::Ptr parse (const char* s) {
  std::stringstream ss(s);
  const auto tree = Node::create();
  read_tree(ss, tree);
  return tree;
}

Node::Ptr make_tensor2d_tree (const Int n, const Node::Ptr& x, const Node::Ptr& y) {
  const auto node = Node::create();
  const bool x_has_kids = ! x->kids.empty();
  const bool y_has_kids = ! y->kids.empty();
  if (x_has_kids) {
    for (const auto& xk : x->kids) {
      if (y_has_kids) {
        for (const auto& yk : y->kids)
          node->kids.push_back(make_tensor2d_tree(n, xk, yk));
      } else {
        node->kids.push_back(make_tensor2d_tree(n, xk, y));
      }
    }
  } else if (y_has_kids) {
    for (const auto& yk : y->kids)
      node->kids.push_back(make_tensor2d_tree(n, x, yk));
  } else {
    for (const auto yi : y->items)
      for (const auto xi : x->items)
        node->items.push_back(n*yi + xi);
  }
  return node;
}

Int get_n (const Node::Ptr& t) {
  Int n = 0;
  for (const auto& k : t->kids) n = std::max(n, get_n(k));
  for (const auto i : t->items) n = std::max(n, i);
  return n;
}

Node::Ptr make_tensor2d_tree (const Node::Ptr& t1d) {
  const Int n = get_n(t1d);
  return make_tensor2d_tree(n+1, t1d, t1d);
}

void fill_elems (const Node::Ptr& node, std::vector<Int>& elems) {
  for (const auto& k : node->kids) fill_elems(k, elems);
  for (const auto i : node->items) elems[i]++;
}

bool covered_once (const Node::Ptr& tree) {
  const Int n = get_n(tree) + 1;
  std::vector<Int> elems(n, 0);
  fill_elems(tree, elems);
  bool ok = true;
  for (const auto e : elems)
    if (e != 1) ok = false;
  return ok;
}

Int encode_tree (const Node::Ptr& n, std::vector<Int>& c, const Int nslots,
                 const Int parent) {
  assert(n->kids.empty() || n->items.empty());
  size_t pos = c.size();
  const Int ni = n->kids.size() + n->items.size();
  c.resize(c.size() + 2 + ni + nslots);
  c[pos] = ni;
  c[pos+1] = parent;
  if (n->kids.empty()) {
    for (Int i = 0; i < ni; ++i)
      c[pos+2+i] = n->items[i];
  } else {
    for (Int i = 0; i < ni; ++i)
      c[pos+2+i] = -encode_tree(n->kids[i], c, nslots, pos);
  }
  return pos;
}

// c follows the structure of the tree described in slmm_mesh.hpp.
void encode_tree (const Node::Ptr& tree, std::vector<Int>& c,
                  const Int nslots=0) {
  c.resize(0);
  encode_tree(tree, c, nslots, -1);
}

Int get_n (const Int* const t, const Int node) {
  const Int nkids = node_nkids(t, node);
  Int n = 0;
  if (node_has_cells(t, node)) {
    for (Int i = 0; i < nkids; ++i)
      n = std::max(n, node_kid(t, node, i));
  } else {
    for (Int i = 0; i < nkids; ++i)
      n = std::max(n, get_n(t, node_kid(t, node, i)));
  }
  return n;
}

void fill_elems (const Int* const t, const Int node,
                 std::vector<Int>& elems) {
  const Int nkids = node_nkids(t, node);
  if (node_has_cells(t, node)) {
    for (Int i = 0; i < nkids; ++i)
      elems[node_kid(t, node, i)]++;
  } else {
    for (Int i = 0; i < nkids; ++i)
      fill_elems(t, node_kid(t, node, i), elems);
  }
}

bool covered_once (const Int* const t) {
  const Int n = get_n(t, 0) + 1;
  std::vector<Int> elems(n, 0);
  fill_elems(t, 0, elems);
  bool ok = true;
  for (const auto e : elems)
    if (e != 1) ok = false;
  return ok;
}

} // namespace

void read_tree (const char* desc, std::vector<Int>& d, const Int nslots) {
  const auto tree = parse(desc);
  encode_tree(tree, d, nslots);
}

void read_tensor2d_tree (const char* desc, std::vector<Int>& d, const Int nslots) {
  const auto tree1d = parse(desc);
  const auto tree2d = make_tensor2d_tree(tree1d);
  encode_tree(tree2d, d, nslots);
}

Int get_max_tree_cell_idx (const Int* tree) {
  return get_n(tree, 0);
}

Int test_array_tree () {
  Int nerr = 0;
  {
    static const char* s0 = "(((0 1 2) (3 4)) ((5 6 7) ((8 9) ((10) (11 12)))))";
    const auto tree = parse(s0);
    const auto s = tree_to_str(tree);
    if (s != std::string(s0)) ++nerr;
    if (get_n(tree) != 12) ++nerr;
    const auto tree2d = make_tensor2d_tree(tree);
    if ( ! covered_once(tree2d)) ++nerr;
    std::vector<Int> c;
    encode_tree(tree2d, c);
    if (get_n(c.data(), 0) != get_n(tree2d)) ++nerr;
    if ( ! covered_once(c.data())) ++nerr;
  }
  {
    static const char* s0 = "(((0 1) (2 3)) ((4 5) (6 7)))";
    const auto tree = parse(s0);
    const auto tree2d = make_tensor2d_tree(tree);
    if ( ! covered_once(tree2d)) ++nerr;
  }
  return nerr;
}

} // namespace tree
} // namespace slmm
