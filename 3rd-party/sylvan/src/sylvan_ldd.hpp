/*
 * Copyright 2011-2016 Formal Methods and Tools, University of Twente
 * Copyright 2016 Tom van Dijk, Johannes Kepler University Linz
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
// Copyright 2020 Wieger Wesselink
//
/// \file mcrl2/3rdparty/sylvan/sylvan_mdd.hpp
/// \brief add your file description here.

#ifndef SYLVAN_LDD_HPP
#define SYLVAN_LDD_HPP

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include <lace.h>
#include <sylvan.h>
#include <sylvan_table.h>
#include <sylvan_ldd.h>
#include <sylvan_int.h>
#include <boost/iterator/iterator_facade.hpp>

static uint64_t cache_merge_id = 0;

namespace sylvan::ldds
{

namespace detail
{

typedef struct __attribute__((packed)) mddnode
{
  uint64_t a, b;
} * mddnode_t; // 16 bytes

inline
mddnode_t LDD_GETNODE(MDD x)
{
  return ((mddnode_t)llmsset_index_to_ptr(nodes, x));
};

inline
std::uint32_t mddnode_getvalue(mddnode_t n)
{
  return *(uint32_t*)((uint8_t*)n+6);
};

inline
MDD mddnode_getdown(mddnode_t n)
{
  return n->b >> 17;
};

/// \brief Forward iterator used for iterating over an ldd cube.
class ldd_iterator: public boost::iterator_facade<ldd_iterator, const std::uint32_t, boost::forward_traversal_tag, std::uint32_t>
{
  private:
    MDD mdd;
    mddnode_t node;

  public:
    explicit ldd_iterator(MDD mdd_)
      : mdd(mdd_)
    {
      assert(mdd != lddmc_false);
      if (mdd != lddmc_true)
      {
        node = LDD_GETNODE(mdd);
      }
    }

  private:
    friend class boost::iterator_core_access;

    void increment()
    {
      mdd = mddnode_getdown(node);
      if (mdd != lddmc_true)
      {
        node = LDD_GETNODE(mdd);
      }
    }

    bool equal(const ldd_iterator& other) const
    {
      return this->mdd == other.mdd;
    }

    std::uint32_t dereference() const
    {
      return mddnode_getvalue(node);
    }
};

} // namespace detail

// From the PhD thesis of Tom van Dijk:
//
// A list decision diagram is a rooted directed acyclic graph with leaves 0 and
// 1. Each internal node has a value v and two outgoing edges labeled > and =,
// also called the "right" and the "down" edge. Along the "right" edges,
// values v
// are encountered in ascending order. The "down" edge never points to leaf 0
// and the "right" edge never points to leaf 1. Duplicate nodes are forbidden.
// LDD nodes have a property called a level (and its dual, depth), which is
// defined as follows: the root node is at the first level, nodes along
// "right" edges
// stay in the same level, while "down" edges lead to the next level. The depth
// of an LDD node is the number of "down" edges to leaf 1. All maximal paths
// from an LDD node have the same depth.

class ldd
{
  friend class Sylvan;

  private:
    MDD mdd;

  public:
    ldd()
     : mdd(lddmc_false)
    {

      lddmc_protect(&mdd);
    }

    explicit ldd(const MDD other)
     : mdd(other)
    {
      lddmc_protect(&mdd);
    }

    ldd(const ldd& other)
     : mdd(other.mdd)
    {
      lddmc_protect(&mdd);
    }

    ldd& operator=(const ldd& other)
    {
      mdd = other.mdd;
      return *this;
    }

    ldd(ldd&& other)
     : mdd(other.mdd)
    {
      lddmc_protect(&mdd);
    }

    ldd& operator=(ldd&& other)
    {
      mdd = other.mdd;
      return *this;
    }

    ~ldd()
    {
      lddmc_unprotect(&mdd);
    }

    MDD get() const
    {
      return mdd;
    }

    // the value
    std::uint32_t value() const
    {
      return lddmc_getvalue(mdd);
    }

    // the down edge
    ldd down() const
    {
      return ldd(lddmc_getdown(mdd));
    }

    // the right edge
    // invariant: right() != true_()
    ldd right() const
    {
      return ldd(lddmc_getright(mdd));
    }

    bool operator==(const ldd& other) const
    {
      return mdd == other.mdd;
    }

    bool operator!=(const ldd& other) const
    {
      return !(*this == other);
    }

    // iterator for a cube in this ldd
    detail::ldd_iterator begin() const
    {
      assert(mdd != lddmc_false);
      return detail::ldd_iterator(mdd);
    }

    detail::ldd_iterator end() const
    {
      return detail::ldd_iterator(lddmc_true);
    }
};

// the leaf false
inline ldd false_()
{
  return ldd(lddmc_false);
}

inline ldd empty_set()
{
  return false_();
}

// the leaf true
inline ldd true_()
{
  return ldd(lddmc_true);
}

inline ldd empty_list()
{
  return true_();
}

// create a node with the given value and down and right arrows
//
// The semantics of this operation are defined as follows:
// [{ε}] = {ε}
// [∅] = ∅
// [node(a, n1, n2)] = {a w | w ∈ [n1]} ∪ [n2]
inline ldd node(std::uint32_t value, const ldd& down = true_(), const ldd& right = false_())
{
  return ldd(lddmc_makenode(value, down.get(), right.get()));
}

// The function extend inserts an edge
//
// extend(false, w, z) = false
// extend(true, w, z) = true
// extend(node(v, x, y), w, z) = v < w -> node(v, x, extend(y, w, z))
//                               v = w -> node(w, z, y)
//                               v > w -> node(w, z, node(v, x, y))
inline ldd extend(const ldd& A, std::uint32_t value, const ldd& B)
{
  return ldd(lddmc_extendnode(A.get(), value, B.get()));
}

// The function follow searches for an edge that has a given value as target
//
// follow(false, w) = false
// follow(true, w) = true
// follow(node(v, x, y), w) = v = w -> x
//                            v > w -> false
//                            v < w -> follow(y, w)
inline ldd follow(const ldd& A, std::uint32_t value)
{
  return ldd(lddmc_follow(A.get(), value));
}

// union(A,B) = A ∪ B
inline ldd union_(const ldd& A, const ldd& B)
{
  LACE_ME;
  return ldd(lddmc_union(A.get(), B.get()));
}

// minus(A,B) = A \ B
inline ldd minus(const ldd& A, const ldd& B)
{
  LACE_ME;
  return ldd(lddmc_minus(A.get(), B.get()));
}

// zip(A,B) = (X,Y) such that X = A U B, Y = B \ A
inline std::pair<ldd, ldd> zip(const ldd& A, const ldd& B)
{
  LACE_ME;
  MDD result[2];
  lddmc_zip(A.get(), B.get(), result);
  return {ldd(result[0]), ldd(result[1]) };
}

// project(A,p) = the projection of the set A according to the projection vector p
inline ldd project(const ldd& A, const ldd& p)
{
  LACE_ME;
  return ldd(lddmc_project(A.get(), p.get()));
}

// project_minus(A,p,B) = minus(project(A,p),B)
inline ldd project_minus(const ldd& A, const ldd& p, const ldd& B)
{
  LACE_ME;
  return ldd(lddmc_project_minus(A.get(), p.get(), B.get()));
}

// intersect(A,B) = A ∩ B
inline ldd intersect(const ldd& A, const ldd& B)
{
  LACE_ME;
  return ldd(lddmc_intersect(A.get(), B.get()));
}

// match(A,B,meta) = A ∩ B, with B defined on subset of the variables of A according to meta
inline ldd match(const ldd& A, const ldd& B, const ldd& meta)
{
  LACE_ME;
  return ldd(lddmc_match(A.get(), B.get(), meta.get()));
}

// join(A,B,pA,pB) = A ∩ B, but A and B are projections of the state vector, described by pA and pB
inline ldd join(const ldd& A, const ldd& B, const ldd& px, const ldd& py)
{
  LACE_ME;
  return ldd(lddmc_join(A.get(), B.get(), px.get(), py.get()));
}

// cube(v) = the singleton set containing the tuple with values in v
inline ldd cube(const std::uint32_t* v, std::size_t n)
{
  LACE_ME;
  return ldd(lddmc_cube(const_cast<std::uint32_t*>(v), n));
}

// cube(v) = the singleton set containing the tuple with values in v
inline ldd cube(const std::vector<std::uint32_t>& v)
{
  return cube(v.data(), v.size());
}

// member_cube(A,v) = check if cube(v) is in the set A
inline ldd member_cube(const ldd& A, const std::vector<std::uint32_t>& v)
{
  LACE_ME;
  return ldd(lddmc_member_cube(A.get(), const_cast<std::uint32_t*>(v.data()), v.size()));
}

// union_cube(A,v) = union_(A,cube(v))
inline ldd union_cube(const ldd& A, const std::uint32_t* v, std::size_t n)
{
  LACE_ME;
  return ldd(lddmc_union_cube(A.get(), const_cast<std::uint32_t*>(v), n));
}

// union_cube(A,v) = union_(A,cube(v))
inline ldd union_cube(const ldd& A, const std::vector<std::uint32_t>& v)
{
  return union_cube(A, v.data(), v.size());
}

// <undocumented>
inline ldd cube_copy(const std::vector<std::uint32_t>& v, const std::vector<int>& copy)
{
  LACE_ME;
  return ldd(lddmc_cube_copy(const_cast<std::uint32_t*>(v.data()), const_cast<int*>(copy.data()), v.size()));
}

// <undocumented>
inline int member_cube_copy(const ldd& A, const std::vector<std::uint32_t>& v, const std::vector<int>& copy)
{
  LACE_ME;
  return lddmc_member_cube_copy(A.get(), const_cast<std::uint32_t*>(v.data()), const_cast<int*>(copy.data()), v.size());
}

// union_cube(A,v) = union_(A,cube(v))
inline ldd union_cube_copy(const ldd& A, const std::uint32_t* v, const int* copy, std::size_t n)
{
  LACE_ME;
  return ldd(lddmc_union_cube_copy(A.get(), const_cast<std::uint32_t*>(v), const_cast<int*>(copy), n));
}

// <undocumented>
inline ldd union_cube_copy(const ldd& A, const std::vector<std::uint32_t>& v, const std::vector<int>& copy)
{
  LACE_ME;
  return union_cube_copy(A, v.data(), copy.data(), v.size());
}

// relprod(A,B,meta) = the successors of the states in A according to the transition relation B which is described by meta
// meta: -1 (end; rest not in rel), 0 (not in rel), 1 (read), 2 (write), 3 (only-read), 4 (only-write), 5 (action label)
inline ldd relprod(const ldd& A, const ldd& B, const ldd& meta)
{
  LACE_ME;
  return ldd(lddmc_relprod(A.get(), B.get(), meta.get()));
}

// relprod_union(A,B,meta) = union_(A,relprod(A,B,meta))
// TODO: the parameter U is undocumented
inline ldd relprod_union(const ldd& A, const ldd& B, const ldd& meta, const ldd& U)
{
  LACE_ME;
  return ldd(lddmc_relprod_union(A.get(), B.get(), meta.get(), U.get()));
}

// relprev(A,B,meta,U) = the predecessors of the states in A according to the transition relation B which is described by meta, restricted to states in U
inline ldd relprev(const ldd& A, const ldd& B, const ldd& meta, const ldd& U)
{
  LACE_ME;
  return ldd(lddmc_relprev(A.get(), B.get(), meta.get(), U.get()));
}

// satcount(A) = the size of the set A
inline double satcount(const ldd& A)
{
  LACE_ME;
  return lddmc_satcount(A.get());
}

// nodecount(A) = the number of nodes of the LDD A
inline
std::size_t nodecount(const ldd& A)
{
  return lddmc_nodecount(A.get());
}

// sat_one(A) = an arbitrary vector from the set A
// returns the value of the leaf node (0 or 1)
template <typename OutputIterator>
inline int sat_one(const ldd& A, OutputIterator to)
{
  LACE_ME;
  MDD mdd = A.get();
  while (mdd != lddmc_false && mdd != lddmc_true)
  {
    detail::mddnode_t n = detail::LDD_GETNODE(mdd);
    *to++ = mddnode_getvalue(n);
    mdd = mddnode_getdown(n);
  }
  return mdd == lddmc_false ? 0 : 1;
}

inline
std::vector<std::uint32_t> sat_one_vector(const ldd& A)
{
  std::vector<std::uint32_t> result;
  sat_one(A, std::back_inserter(result));
  return result;
}

// sat_all_nopar(A,cb) calls the callback function cb for each vector in A
inline
void sat_all(const ldd& A, lddmc_enum_cb cb, void* context = nullptr)
{
  LACE_ME;
  lddmc_sat_all_par(A.get(), cb, context);
}

// sat_all_par(A,cb) = sat_all(A,cb), but parallelized
inline
void sat_all_nopar(const ldd& A, lddmc_enum_cb cb, void* context = nullptr)
{
  LACE_ME;
  lddmc_sat_all_nopar(A.get(), cb, context);
}

// collect(A,cb) = sat_all_par(A,cb), but the callback cb returns some set encoded as an LDD, and all returned LDDs are combined using union
// TODO: this is probably wrong, because there are two more undocumented parameters
inline ldd collect(const ldd& A, lddmc_collect_cb cb, void* context = nullptr)
{
  LACE_ME;
  return ldd(lddmc_collect(A.get(), cb, context));
}

// match_sat(A,B,m,cb) = sat_all_par(match(A,B,m),cb)
inline
void match_sat(const ldd& A, const ldd& B, const ldd& meta, lddmc_enum_cb cb, void* context = nullptr)
{
  LACE_ME;
  lddmc_match_sat_par(A.get(), B.get(), meta.get(), cb, context);
}

inline ldd pick_cube(const ldd& A)
{
  return ldd(lddmc_pick_cube(A.get()));
}

inline
void print_dot(const std::string& filename, const ldd& A)
{
  FILE* out = fopen(filename.c_str(), "w");
  if (!out)
  {
    throw std::runtime_error("Could not open file " + filename);
  }
  lddmc_fprintdot(out, A.get());
}

inline
void ldd_solutions_callback(WorkerP*, Task*, std::uint32_t* v, std::size_t n, void* context = nullptr)
{
  std::vector<std::vector<std::uint32_t>>& V = *reinterpret_cast<std::vector<std::vector<std::uint32_t>>*>(context);
  V.push_back(std::vector<std::uint32_t>(v, v + n));
}

inline
std::vector<std::vector<std::uint32_t>> ldd_solutions(const ldds::ldd& x)
{
  std::vector<std::vector<std::uint32_t>> result;
  sat_all_nopar(x, ldd_solutions_callback, &result);
  return result;
}

inline
std::string print_ldd(const ldds::ldd& x)
{
  std::ostringstream out;
  auto solutions = ldd_solutions(x);
  bool multiline = solutions.size() > 1;
  std::string sep = multiline ? ",\n" : ", ";

  out << "{" << (multiline ? "\n" : " ");
  for (auto i = solutions.begin(); i != solutions.end(); ++i)
  {
    const std::vector<std::uint32_t>& solution = *i;
    if (i != solutions.begin())
    {
      out << sep;
    }
    out << "[ ";
    for (auto j = solution.begin(); j != solution.end(); ++j)
    {
      if (j != solution.begin())
      {
        out << ", ";
      }
      out << *j;
    }
    out << " ]";
  }
  out << (multiline ? "\n" : " ") << "}";

  return out.str();
}

inline
std::ostream& operator<<(std::ostream& out, const ldds::ldd& x)
{
  return out << ldds::print_ldd(x);
}

// Returns { x in X | x[0] = value }
inline
ldd fix_first_element(const ldd& X, std::uint32_t value)
{
  ldd x = X;
  while (x.value() != value)
  {
    x = x.right();
    if (x == false_())
    {
      return empty_set();
    }
  }
  return node(value, x.down());
}

// Initialise caches for custom LDD operations.
void initialise() {
  cache_merge_id = cache_next_opid();
}

ldd merge(ldd a, ldd b) {
  if (a == ldds::empty_list()) {
    return b;
  } else if (b == ldds::empty_list()) {
    return a;
  } else if (a == ldds::empty_set() || b == ldds::empty_set()) {
    return ldds::empty_set();
  } else {
    
    MDD result;
    if (cache_get(cache_merge_id, a.get(), b.get(), &result)) 
    {
      return ldd(result);
    }

    ldd new_result = node(a.value(), merge(b.down(), a), merge(a.right(), b));

    cache_put(cache_merge_id, a.get(), b.get(), new_result.get());
    return new_result;
  }

}

// Computes the meta for relprod and relprev
// read = the indices of read variables
// write = the indices of write variables
// has_action: if true, then a 5 is added to indicate there is an extra parameter that corresponds to an action
inline
ldd compute_meta(const std::vector<std::size_t>& read, const std::vector<std::size_t>& write, bool has_action = false)
{
  std::vector<std::uint32_t> meta;

  std::size_t r_k = read.size();
  std::size_t w_k = write.size();
  std::size_t r_i = 0;
  std::size_t w_i = 0;
  std::size_t i=0;

  for (;;)
  {
    // compute the type (0=read+write, 1=read, 2=write, 3=only-read, 4=only-write (???)
    std::size_t type = 0;
    if (r_i < r_k && read[r_i] == i)
    {
      r_i++;
      type += 1; // read
    }
    if (w_i < w_k && write[w_i] == i)
    {
      w_i++;
      type += 2; // write
    }

    // now that we have type, set meta */
    if (type == 0) meta.push_back(0);
    else if (type == 1) { meta.push_back(3); }
    else if (type == 2) { meta.push_back(4); }
    else if (type == 3) { meta.push_back(1); meta.push_back(2); }
    if (r_i == r_k && w_i == w_k)
    {
      if (has_action)
      {
        meta.push_back(5);
      }
      meta.push_back((std::uint32_t)-1);
      break;
    }
    i++;
  }
  return cube(meta);
}

// Returns true if V is a subset of U
inline
bool includes(const ldds::ldd& U, const ldds::ldd& V)
{
  return union_(U, V) == U;
}

} // namespace sylvan::ldds

#endif // SYLVAN_LDD_OBJ_H
