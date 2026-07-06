// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/alphabet_operations.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_ALPHABET_OPERATIONS_H
#define MCRL2_PROCESS_ALPHABET_OPERATIONS_H

#include "mcrl2/core/identifier_string.h"
#include "mcrl2/process/communication_expression.h"
#include "mcrl2/process/multi_action_name.h"
#include "mcrl2/process/action_names.h"
#include "mcrl2/process/process_expression.h"
#include "mcrl2/process/rename_expression.h"
#include "mcrl2/utilities/detail/atomic_wrapper.h"
#include "mcrl2/utilities/sequence.h"
#include <algorithm>
#include <iostream>
#include <iterator>

namespace mcrl2::process::alphabet_operations
{

//-----------------------------------------------------//
// block operations
//-----------------------------------------------------//

inline
multi_action_name_set block(const core::identifier_string_list& B, const multi_action_name_set& A, bool A_includes_subsets = false)
{
  multi_action_name_set result;
  multi_action_name beta(B.begin(), B.end());

  if (A_includes_subsets)
  {
    for (multi_action_name alpha: A)
    {
      for (const core::identifier_string& b: B)
      {
        alpha.erase(b);
      }
      if (!alpha.empty())
      {
        result.insert(alpha);
      }
    }
  }
  else
  {
    for (const multi_action_name& alpha: A)
    {
      if (has_empty_intersection(beta, alpha))
      {
        result.insert(alpha);
      }
    }
  }
  return result;
}

//-----------------------------------------------------//
// hide operations
//-----------------------------------------------------//

// Hides elements in I from alpha
inline
multi_action_name hide(const multi_action_name& alpha, const std::set<core::identifier_string>& I)
{
  using utilities::detail::contains;
  multi_action_name result;
  for (const core::identifier_string& a: alpha)
  {
    if (!contains(I, a))
    {
      result.insert(a);
    }
  }
  return result;
}

// Hides (or: removes) elements in I from C
inline
std::set<core::identifier_string> hide(const core::identifier_string_list& I, const std::set<core::identifier_string>& J)
{
  using utilities::detail::contains;
  std::set<core::identifier_string> result;
  for (const core::identifier_string& j: J)
  {
    if (!contains(I, j))
    {
      result.insert(j);
    }
  }
  return result;
}

inline
multi_action_name hide(const std::set<core::identifier_string>& I, const multi_action_name& alpha)
{
  using utilities::detail::contains;
  multi_action_name result;
  for (const core::identifier_string& i: alpha)
  {
    if (!contains(I, i))
    {
      result.insert(i);
    }
  }
  return result;
}

template <typename IdentifierContainer>
multi_action_name_set hide(const IdentifierContainer& I, const multi_action_name_set& A, bool /* A_includes_subsets */ = false)
{
  multi_action_name m(I.begin(), I.end());
  multi_action_name_set result;
  for (multi_action_name alpha: A)
  {
    for (const core::identifier_string& i: I)
    {
      alpha.erase(i);
    }
    result.insert(alpha);
  }
  return result;
}

inline
multi_action_name_set hide_inverse(const core::identifier_string_list& I, const multi_action_name_set& A, bool A_includes_subsets = false)
{
  return block(I, A, A_includes_subsets);
}

//-----------------------------------------------------//
// merge/left_merge/sync operations
//-----------------------------------------------------//

inline
multi_action_name_set concat(const multi_action_name_set& A1, const multi_action_name_set& A2)
{
  multi_action_name_set result;
  for (const multi_action_name& i: A1)
  {
    for (const multi_action_name& j: A2)
    {
       result.insert(multiset_union(i, j));
    }
  }
  return result;
}

// Returns true if alpha in concat(A1, A2).
// TODO: Make the implementation more inefficient!
inline
bool concat_includes(const multi_action_name_set& A1, const multi_action_name_set& A2, const multi_action_name& alpha)
{
  for (const multi_action_name& alpha1: A1)
  {
    for (const multi_action_name& alpha2: A2)
    {
      if (alpha == multiset_union(alpha1, alpha2))
      {
        return true;
      }
    }
  }
  return false;
}

// Returns the intersection of concat(A1, A2) and subsets(A)
inline
multi_action_name_set bounded_concat(const multi_action_name_set& A1, const multi_action_name_set& A2, const multi_action_name_set& A)
{
  multi_action_name_set result;
  for (const multi_action_name& i: A1)
  {
    for (const multi_action_name& j: A2)
    {
      multi_action_name alpha = multiset_union(i, j);
      if (subset_includes(A, alpha))
      {
        result.insert(alpha);
      }
    }
  }
  return result;
}

inline
multi_action_name_set left_arrow1(const multi_action_name_set& A1, const multi_action_name_set& A2)
{
  multi_action_name_set result = A1; // needed because tau is not explicitly stored
  for (const multi_action_name& beta: A2)
  {
    for (const multi_action_name& gamma: A1)
    {
      if (includes(gamma, beta))
      {
        multi_action_name alpha = multiset_difference(gamma, beta);
        if (!alpha.empty())
        {
          result.insert(alpha);
        }
      }
    }
  }
  return result;
}

inline
multi_action_name_set left_arrow(const multi_action_name_set& A1, bool A1_includes_subsets, const multi_action_name_set& A2)
{
  multi_action_name_set result;
  if (A1_includes_subsets)
  {
    result = A1;
  }
  else
  {
    result = set_union(A1, left_arrow1(A1, A2));
  }
  return result;
}

// returns left_arrow(A I*, A2)
inline
multi_action_name_set left_arrow2(const multi_action_name_set& A, const std::set<core::identifier_string>& I, const multi_action_name_set& A2)
{
  multi_action_name_set result = A; // needed because tau is not explicitly stored
  for (const multi_action_name& alpha2: A2)
  {
    multi_action_name beta = hide(I, alpha2);
    for (const multi_action_name& alpha: A)
    {
      if (includes(alpha, beta))
      {
        multi_action_name gamma = multiset_difference(alpha, beta);
        if (!gamma.empty())
        {
          result.insert(hide(I, gamma));
        }
      }
    }
  }
  return result;
}

inline
multi_action_name_set merge(const multi_action_name_set& A1, const multi_action_name_set& A2)
{
  return set_union(set_union(A1, A2), concat(A1, A2));
}

// Returns the intersection of merge(A1, A2) and subsets(A)
inline
multi_action_name_set bounded_merge(const multi_action_name_set& A1, const multi_action_name_set& A2, const multi_action_name_set& A)
{
  return set_union(set_union(A1, A2), bounded_concat(A1, A2, A));
}

inline
multi_action_name_set left_merge(const multi_action_name_set& A1, const multi_action_name_set& A2)
{
  return merge(A1, A2);
}

inline
multi_action_name_set sync(const multi_action_name_set& A1, const multi_action_name_set& A2)
{
  return concat(A1, A2);
}

// Returns the intersection of sync(A1, A2) and subsets(A)
inline
multi_action_name_set bounded_sync(const multi_action_name_set& A1, const multi_action_name_set& A2, const multi_action_name_set& A)
{
  return bounded_concat(A1, A2, A);
}

//-----------------------------------------------------//
// comm operations
//-----------------------------------------------------//

/// Returns true is all action a in alpha are in action_names.
inline
bool includes(const action_name_set& action_names, const multi_action_name& alpha)
{
  return std::all_of(alpha.begin(), alpha.end(), [&action_names](const core::identifier_string& a) {
    return utilities::detail::contains(action_names, a);
  });
}

/// Cache communication expressions by their name, to speed up the implementation of comm_inverse.
using comm_inverse_cache = std::unordered_map<core::identifier_string, std::vector<communication_expression>>;

/// Calculate the inverse of a communication expression list C. The result is cached
inline
comm_inverse_cache calculate_comm_inverse_cache(const communication_expression_list& C)
{
  comm_inverse_cache cache;
  for (const communication_expression& c: C)
  {
    cache[c.name()].push_back(c);
  }
  return cache;
}

/// Calculates the inverse of a communication expression on alpha1, and accumulates
/// expressions to which the inverse of communication have been applied in alpha2.
/// \param C_inverse a cache of communication expressions, indexed by their name.
/// \param action_names the set of action names that are allowed to be used in the result.
/// \param alpha1 the multi action name to which the inverse of communication expressions are applied.
/// \param alpha2 the multi action name to which the inverse of communication expressions have been applied.
/// \param A_includes_subsets if true, in case not all action names in \f$ \alpha = \alpha1 \cup \alpha2\f$ are in action_names, the result will include the maximum subset of \f$ \alpha \f$ of actions in action_names
/// \param result the set of multi action names that are the result of applying the inverse of communication expressions to alpha1 and alpha2.
inline
void comm_inverse(const comm_inverse_cache& C_inverse, const action_name_set& action_names, const multi_action_name& alpha1, const multi_action_name& alpha2, const bool A_includes_subsets, multi_action_name_set& result)
{
  mCRL2log(log::trace) << "comm_inverse: action_names = " << core::pp(action_names) << ", alpha1 = " << pp(alpha1) << ", alpha2 = " << pp(alpha2) << "\n";

  if (includes(action_names, alpha1) && includes(action_names, alpha2))
  {
    mCRL2log(log::trace) << "comm_inverse: result += {" << multiset_union(alpha1, alpha2) << "}\n";
    result.insert(multiset_union(alpha1, alpha2));
  }
  else if (A_includes_subsets)
  {
    // keep the largest subset that is included in action_names
    multi_action_name alpha;
    for (const core::identifier_string& a: alpha1)
    {
      if (utilities::detail::contains(action_names, a))
      {
        alpha.insert(a);
      }
    }
    for (const core::identifier_string& a: alpha2)
    {
      if (utilities::detail::contains(action_names, a))
      {
        alpha.insert(a);
      }
    }
    result.insert(alpha);
    mCRL2log(log::trace) << "comm_inverse: result += {" << alpha << "}\n";
  }

  if (!alpha1.empty())
  {
    multi_action_name beta1(alpha1);
    multi_action_name::const_iterator beta1_last_it = beta1.nth(beta1.size() - 1);
    core::identifier_string beta1_last = *beta1_last_it;
    beta1.erase(beta1_last_it);

    /// Case where action name is not the result of communication
    if (action_names.find(beta1_last) != action_names.end())
    {
      multi_action_name beta2(alpha2);
      beta2.insert(beta1_last);
      comm_inverse(C_inverse, action_names, beta1, beta2, A_includes_subsets, result);
    }

    comm_inverse_cache::const_iterator i = C_inverse.find(beta1_last);
    if (i != C_inverse.end())
    {
      for (const communication_expression& c: i->second)
      {
        multi_action_name beta2(alpha2);
        beta2.insert(c.action_name().names().begin(), c.action_name().names().end());
        comm_inverse(C_inverse, action_names, beta1, beta2, A_includes_subsets, result);
      }
    }
  }
}

/// Returns the set of multiactions that, with application of communication expressions in C, can result in a multiaction in A.
/// \param C a list of communication expressions
/// \param action_names the set of action names that are allowed to be used in the result
/// \param A a set of multi action names
inline multi_action_name_set comm_inverse1(
  const communication_expression_list& C,
  const action_name_set& action_names,
  const multi_action_name_set& A,
  const bool A_includes_subsets = false)
{
  mCRL2log(log::trace) << "comm_inverse1: C = " << core::pp(C) << ", action_names = " << core::pp(action_names) << ", A = " << pp(A) << "\n";
  multi_action_name_set result;
  multi_action_name empty;
  const comm_inverse_cache C_inverse = calculate_comm_inverse_cache(C);
  for (multi_action_name alpha: A)
  {
    comm_inverse(C_inverse, action_names, alpha, empty, A_includes_subsets, result);
  }
  mCRL2log(log::trace) << "comm_inverse1: result = " << pp(result) << "\n";
  return result;
}

// Note that the result is flattened.
inline
std::set<core::identifier_string> comm_inverse(const communication_expression_list& C, const action_name_set& action_names, const std::set<core::identifier_string>& I)
{
  const comm_inverse_cache C_inverse = calculate_comm_inverse_cache(C);
  std::set<core::identifier_string> result = I;
  for (const core::identifier_string& i: I)
  {
    comm_inverse_cache::const_iterator C_inverse_it = C_inverse.find(i);
    if (C_inverse_it != C_inverse.end())
    {
      for (const communication_expression& j: C_inverse_it->second)
      {
        std::set<core::identifier_string> lhs(j.action_name().names().begin(), j.action_name().names().end());
        if (utilities::detail::set_includes(action_names, lhs))
        {
          result.insert(lhs.begin(), lhs.end());
        }
      }
    }
  }
  return result;
}

inline
void apply_comm(const communication_expression& c, multi_action_name_set& A)
{
  core::identifier_string_list names = c.action_name().names();
  const core::identifier_string& a = c.name();
  multi_action_name alpha(names.begin(), names.end());
  // c == alpha -> a

  multi_action_name_set to_be_added;
  for (multi_action_name beta: A)
  {
    while (includes(beta, alpha))
    {
      for (const core::identifier_string& a: alpha)
      {
        beta.erase(beta.find(a));
      }
      beta.insert(a);
      to_be_added.insert(beta);
    }
  }
  A.insert(to_be_added.begin(), to_be_added.end());
}

inline
multi_action_name_set comm(const communication_expression_list& C, const multi_action_name_set& A, bool /* A_includes_subsets */ = false)
{
  multi_action_name_set result = A;

  // sequentially apply the communication rules to result
  for (const communication_expression& c: C)
  {
    apply_comm(c, result);
  }

  return result;
}

inline
communication_expression_list filter_comm_set(const communication_expression_list& C, const multi_action_name_set& alphabet)
{
  std::vector<communication_expression> result;
  for (const communication_expression& c: C)
  {
    core::identifier_string_list lhs = c.action_name().names();
    multi_action_name alpha(lhs.begin(), lhs.end());
    if (includes(alphabet, alpha))
    {
      result.push_back(c);
    }
  }
  return communication_expression_list(result.begin(), result.end());
}

//-----------------------------------------------------//
// rename operations
//-----------------------------------------------------//

// returns true if x is a source of one of the rename rules in R
inline
bool is_source(const rename_expression_list& R, const core::identifier_string& x)
{
  for (const rename_expression& r: R)
  {
    if (r.source() == x)
    {
      return true;
    }
  }
  return false;
}

inline
core::identifier_string apply_rename(const rename_expression_list& R, const core::identifier_string& x)
{
  for (const rename_expression& r: R)
  {
    if (x == r.source())
    {
      return r.target();
    }
  }
  return x;
}

inline
multi_action_name apply_rename(const rename_expression_list& R, const multi_action_name& a)
{
  multi_action_name result;
  for (const core::identifier_string& i: a)
  {
    result.insert(apply_rename(R, i));
  }
  return result;
}

/// \brief Computes R(A)
inline
multi_action_name_set rename(const rename_expression_list& R, const multi_action_name_set& A, bool /* A_includes_subsets */ = false)
{
  multi_action_name_set result;
  for (const multi_action_name& alpha: A)
  {
    result.insert(apply_rename(R, alpha));
  }
  return result;
}

using rename_inverse_map = std::map<core::identifier_string, std::vector<core::identifier_string>>;

// Example: R = {b -> c}, then rename_inverse(R) = {b -> [], c -> [b, c]}
inline
rename_inverse_map rename_inverse(const rename_expression_list& R)
{
  rename_inverse_map Rinverse;
  for (auto i = R.begin(); i != R.end(); ++i)
  {
    Rinverse[i->target()].push_back(i->source());
    if (!is_source(R, i->target()))
    {
      Rinverse[i->target()].push_back(i->target());
    }
    Rinverse[i->source()]; // this is to make sure that i->source() is in the map
  }
  return Rinverse;
}

inline
std::string print_rename_inverse_map(const rename_inverse_map& m)
{
  std::ostringstream out;
  out << "{ ";
  for (auto i = m.begin(); i != m.end(); ++i)
  {
    if (i != m.begin())
    {
      out << ", ";
    }
    out << core::pp(i->first) << " -> [";
    const std::vector<core::identifier_string>& v = i->second;
    for (auto j = v.begin(); j != v.end(); ++j)
    {
      if (j != v.begin())
      {
        out << ", ";
      }
      out << *j;
    }
    out << "]";
  }
  out << " }";
  return out.str();
}

struct rename_inverse_apply
{
  const multi_action_name& alpha;
  const std::vector<core::identifier_string>& beta;
  multi_action_name_set& A;

  rename_inverse_apply(const multi_action_name& alpha_, const std::vector<core::identifier_string>& beta_, multi_action_name_set& A_)
    : alpha(alpha_), beta(beta_), A(A_)
  {}

  void operator()()
  {
    multi_action_name gamma = alpha;
    gamma.insert(beta.begin(), beta.end());
    A.insert(gamma);
  }
};

inline
std::set<core::identifier_string> rename_inverse(const rename_expression_list& R, const std::set<core::identifier_string>& I)
{
  alphabet_operations::rename_inverse_map Rinverse = alphabet_operations::rename_inverse(R);

  std::set<core::identifier_string> result;
  for (const core::identifier_string& i: I)
  {
    auto j = Rinverse.find(i);
    if (j != Rinverse.end())
    {
      result.insert(j->second.begin(), j->second.end());
    }
    else
    {
      result.insert(i);
    }
  }
  return result;
}

inline
void rename_inverse(const rename_inverse_map& Rinverse, const multi_action_name& x, bool x_includes_subsets, multi_action_name_set& result)
{
  std::vector<std::vector<core::identifier_string> > V;

  multi_action_name alpha = x;

  // remove elements that appear in Rinverse, and put the replacements in V
  for (auto i = alpha.begin(); i != alpha.end(); )
  {
    auto j = Rinverse.find(*i);
    if (j != Rinverse.end())
    {
      i = alpha.erase(i);
      if (!j->second.empty() || !x_includes_subsets)
      {
        V.push_back(j->second);
      }
    }
    else
    {
      ++i;
    }
  }

  // v will hold a replacement
  std::vector<core::identifier_string> v(V.size());

  // generate all permutations of the replacements in V, and put them in result
  utilities::foreach_sequence(V, v.begin(), rename_inverse_apply(alpha, v, result));
}

/// \brief Computes R^[-1}(A)
inline
multi_action_name_set rename_inverse(const rename_expression_list& R, const multi_action_name_set& A, bool A_includes_subsets = false)
{
  rename_inverse_map Rinverse = rename_inverse(R);

  multi_action_name_set result;
  for (const multi_action_name& alpha: A)
  {
    rename_inverse(Rinverse, alpha, A_includes_subsets, result);
  }
  return result;
}

//-----------------------------------------------------//
// allow operations
//-----------------------------------------------------//

inline
multi_action_name_set allow(const action_name_multiset_list& V, const multi_action_name_set& A, bool A_includes_subsets = false)
{
  multi_action_name_set result;

  for (const action_name_multiset& s: V)
  {
    const core::identifier_string_list& names = s.names();
    multi_action_name v(names.begin(), names.end());
    for (const multi_action_name& alpha: A)
    {
      bool keep = A_includes_subsets ? includes(alpha, v) : alpha == v;
      if (keep)
      {
        result.insert(v);
      }
    }
  }
  return result;
}

} // namespace mcrl2::process::alphabet_operations

#endif // MCRL2_PROCESS_ALPHABET_OPERATIONS_H
