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
/// \file mcrl2/3rdparty/sylvan/sylvan_bdd.hpp
/// \brief add your file description here.

#ifndef SYLVAN_BDD_HPP
#define SYLVAN_BDD_HPP

#include <string>
#include <vector>

#include <lace.h>
#include <sylvan.h>
#include <sylvan_bdd.h>

namespace sylvan::bdds
{

class bdd
{
  friend class bdd_substitution;

  private:
    BDD m_bdd;

  public:
    bdd()
     : m_bdd(sylvan_false)
    {
      sylvan_protect(&m_bdd);
    }

    explicit bdd(const BDD from)
     : m_bdd(from)
    {
      sylvan_protect(&m_bdd);
    }

    bdd(const bdd& from)
     : m_bdd(from.m_bdd)
    {
      sylvan_protect(&m_bdd);
    }

    explicit bdd(const std::uint32_t var)
     : m_bdd(sylvan_ithvar(var))
    {
      sylvan_protect(&m_bdd);
    }

    bdd& operator=(const bdd& other) = default;

    ~bdd()
    {
      sylvan_unprotect(&m_bdd);
    }

    BDD get() const
    {
      return m_bdd;
    }
    
    int is_true() const
    {
      return m_bdd == sylvan_true;
    }

    bool is_false() const
    {
      return m_bdd == sylvan_false;
    }

    /// \brief Returns true if this bdd has the value true or false
    bool is_constant() const
    {
      return is_true() || is_false();
    }

    /// \brief Returns true if this bdd has the value true or false
    bool is_terminal() const
    {
      return is_constant();
    }

    /// \brief Returns the top variable index of this bdd (the variable in the root node)
    std::uint32_t top() const
    {
      return sylvan_var(m_bdd);
    }

    /// \brief Follows the high edge ("then") of the root node of this bdd
    bdd then() const
    {
      return bdd(sylvan_high(m_bdd));
    }

    /// \brief Follows the low edge ("else") of the root node of this bdd
    bdd else_() const
    {
      return bdd(sylvan_low(m_bdd));
    }

    /// \brief Returns the number of nodes in this bdd. Not thread-safe!
    std::size_t node_count() const
    {
      return sylvan_nodecount(m_bdd);
    }

    bool operator==(const bdd& other) const
    {
      return m_bdd == other.m_bdd;
    }

    bool operator!=(const bdd& other) const
    {
      return m_bdd != other.m_bdd;
    }

    bool operator<=(const bdd& other) const
    {
      // TODO: better implementation, since we are not interested in the BDD result
      return sylvan_ite(m_bdd, sylvan_not(other.m_bdd), sylvan_false) == sylvan_false;
    }

    bool operator>=(const bdd& other) const
    {
      // TODO: better implementation, since we are not interested in the BDD result
      return other <= *this;
    }

    bool operator<(const bdd& other) const
    {
      return m_bdd != other.m_bdd && *this <= other;
    }

    bool operator>(const bdd& other) const
    {
      return m_bdd != other.m_bdd && *this >= other;
    }

    bdd operator!() const
    {
      return bdd(sylvan_not(m_bdd));
    }

    bdd operator~() const
    {
      return bdd(sylvan_not(m_bdd));
    }

    bdd operator*(const bdd& other) const
    {
      return bdd(sylvan_and(m_bdd, other.m_bdd));
    }

    bdd& operator*=(const bdd& other)
    {
      m_bdd = sylvan_and(m_bdd, other.m_bdd);
      return *this;
    }

    bdd operator&(const bdd& other) const
    {
      return bdd(sylvan_and(m_bdd, other.m_bdd));
    }

    bdd& operator&=(const bdd& other)
    {
      m_bdd = sylvan_and(m_bdd, other.m_bdd);
      return *this;
    }

    bdd operator+(const bdd& other) const
    {
      return bdd(sylvan_or(m_bdd, other.m_bdd));
    }

    bdd& operator+=(const bdd& other)
    {
      m_bdd = sylvan_or(m_bdd, other.m_bdd);
      return *this;
    }

    bdd operator|(const bdd& other) const
    {
      return bdd(sylvan_or(m_bdd, other.m_bdd));
    }

    bdd& operator|=(const bdd& other)
    {
      m_bdd = sylvan_or(m_bdd, other.m_bdd);
      return *this;
    }

    bdd operator^(const bdd& other) const
    {
      return bdd(sylvan_xor(m_bdd, other.m_bdd));
    }

    bdd& operator^=(const bdd& other)
    {
      m_bdd = sylvan_xor(m_bdd, other.m_bdd);
      return *this;
    }

    bdd operator-(const bdd& other) const
    {
      return bdd(sylvan_and(m_bdd, sylvan_not(other.m_bdd)));
    }

    bdd& operator-=(const bdd& other)
    {
      m_bdd = sylvan_and(m_bdd, sylvan_not(other.m_bdd));
      return *this;
    }
};

/// \brief Creates a bdd representing just the variable index in its positive form
/// The variable index must be a 0<=index<=2^23 (we use 24 bits internally)
inline bdd make_variable(std::uint32_t index)
{
  return bdd(sylvan_ithvar(index));
}

/// \brief Returns the bdd representing "True"
inline bdd true_()
{
  return bdd(sylvan_true);
}

/// \brief Returns the bdd representing "False"
inline bdd false_()
{
  return bdd(sylvan_false);
}

/// \brief Computes \exists v: x
/// The set of variables v is represented by a conjunction of variables
inline bdd exists(const bdd& v, const bdd& x)
{
  return bdd(sylvan_exists(x.get(), v.get()));
}

/// \brief  Computes \forall v: x
/// The set of variables v is represented by a conjunction of variables
inline bdd forall(const bdd& v, const bdd& x)
{
  return bdd(sylvan_forall(x.get(), v.get()));
}

/// \brief  Computes if b then x else y
inline bdd ite(const bdd& b, const bdd& x, const bdd& y)
{
  return bdd(sylvan_ite(b.get(), x.get(), y.get()));
}

/// \brief Computes x <=> y
inline
bdd equiv(const bdd& x, const bdd& y)
{
  return bdd(sylvan_equiv(x.get(), y.get()));
}

/// \brief Computes x => y
inline
bdd implies(const bdd& x, const bdd& y)
{
  return (!x) | y;
}

inline double satcount(const bdd& x, const bdd &variables)
{
  return sylvan_satcount(x.get(), variables.get());
}

class bdd_substitution
{
  protected:
    BDD m_bdd;

  public:
    bdd_substitution()
     : m_bdd(sylvan_map_empty())
    {
      sylvan_protect(&m_bdd);
    }

    // TODO: copying is disabled, because it leads to crashes It's not clear
    // why, because this class is very similar to the class bdd.
    bdd_substitution(const bdd_substitution& other) = delete;
    ~bdd_substitution()
    {
      sylvan_unprotect(&m_bdd);
    }

    BDD get() const
    {
      return m_bdd;
    }

    /// \brief Adds a key-value pair to the substitution
    /// \param index The index of a bdd variable
    void put(std::uint32_t index, const bdd& value)
    {
      m_bdd = sylvan_map_add(m_bdd, index, value.get());
    }
};

/// \brief Applies the substitution sigma to x.
inline bdd let(const bdd_substitution& sigma, const bdd& x)
{
  return bdd(sylvan_compose(x.get(), sigma.get()));
}

/// \brief  Computes the reverse application of a transition relation to this set.
/// @param relation the transition relation to apply
/// @param v the variables that are in the transition relation
/// This function assumes that s,t are interleaved with s even and t odd (s+1).
/// Other variables in the relation are ignored (by existential quantification)
/// Set v to "false" (illegal v) to assume all encountered variables are in s,t
///
/// Use this function to concatenate two relations   --> -->
/// or to take the 'previous' of a set               -->  S
inline bdd relprev(const bdd& relation, const bdd& x, const bdd& variables)
{
  return bdd(sylvan_relprev(relation.get(), x.get(), variables.get()));
}

/// \brief  Computes the application of a transition relation to this set.
/// @param relation the transition relation to apply
/// @param variables the variables that are in the transition relation
/// This function assumes that s,t are interleaved with s even and t odd (s+1).
/// Other variables in the relation are ignored (by existential quantification)
/// Set v to "false" (illegal v) to assume all encountered variables are in s,t
///
/// Use this function to take the 'next' of a set     S  -->
inline bdd relnext(const bdd& relation, const bdd& x, const bdd& variables)
{
  return bdd(sylvan_relnext(relation.get(), x.get(), variables.get()));
}

// Applies relation R to x
bdd relation_forward(const bdd& R, const bdd& x, const bdd& variables, const bdd_substitution& prev_substitution, bool optimized = false)
{
  if (optimized)
  {
    // Note that the parameter prev_substitution is not used. This is because Sylvan makes assumptions about the variable order
    return relnext(R, x, variables);
  }
  return let(prev_substitution, exists(variables, x & R));
};

// Applies the inverse of the relation R to x
bdd relation_backward(const bdd& R, const bdd& x, const bdd& next_variables, const bdd_substitution& next_substitution, bool optimized = false)
{
  if (optimized)
  {
    // Note that the parameter next_substitution is not used. This is because Sylvan makes assumptions about the variable order
    return relprev(R, x, next_variables);
  }
  return exists(next_variables, let(next_substitution, x) & R);
};

inline bdd all(const std::vector<bdd>& v)
{
  if (v.empty())
  {
    return true_();
  }
  bdd result = v[0];
  for (auto i = ++v.begin(); i != v.end(); ++i)
  {
    result &= *i;
  }
  return result;
}

inline bdd any(const std::vector<bdd>& v)
{
  if (v.empty())
  {
    return false_();
  }
  bdd result = v[0];
  for (auto i = ++v.begin(); i != v.end(); ++i)
  {
    result |= *i;
  }
  return result;
}

} // namespace sylvan::bdds

#endif // SYLVAN_LDD_OBJ_H
