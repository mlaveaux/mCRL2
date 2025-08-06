// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/substitutions/mutable_indexed_substitution.h
/// \brief add your file description here.

// The code below contains an experiment to replace the classical mutable_indexed_substitution
// by a std::unordered_map from variables to expressions. This is up to 1.5 times 
// slower in most state space generations when there are not too many variables, which is
// typical for state space generation without complex sum operations or quantifiers.
// When a large number of variables exist, generally generated as fresh variables,
// std::unordered_map can perform much better, leading to the time to generate a 
// state space with a factor 2. 


#ifndef MCRL2_DATA_SUBSTITUTIONS_MUTABLE_INDEXED_SUBSTITUTION_H
#define MCRL2_DATA_SUBSTITUTIONS_MUTABLE_INDEXED_SUBSTITUTION_H

#include "mcrl2/data/data_expression.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/standard_containers/unordered_map.cxx"
  #include "mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.cxx"
#else
  import atermpp;
#endif

namespace mcrl2::data {

/// \brief Generic substitution function.
/// \details This substitution assumes a function variable -> std::size_t, that, for
///          each variable gives a unique index. The substitutions are stored
///          internally as a vector, mapping std::size_t to expression.
///          Provided that, given a variable, its index can be computed in O(1)
///          time, insertion is O(1) amortized, and lookup is O(1).
///          Memory required is O(n) where n is the largest index used.
///            This substitution can maintain the variables occurring in the rhs of
///          expressions. If it is requested that whether a variable occurs in a rhs, 
///          the substitution automatically maintains these variables. This requires
///          time O(n log n) per rhs, where n is the size of the rhs. 
template <typename VariableType = data::variable, typename ExpressionType = data_expression >
class mutable_indexed_substitution
{
protected:
  using substitution_type = atermpp::utilities::unordered_map<VariableType, ExpressionType>;

  /// \brief Internal storage for substitutions.
  /// Required to be a container with random access through [] operator.
  /// It is essential to store the variable also in the container, as it might be that
  /// this variable is not used anywhere although it has a valid assignment. This happens
  /// for instance when the assignment is already parsed, while the expression to which it
  /// needs to be applied must still be parsed. 
  substitution_type m_substitution;
  mutable bool m_variables_in_rhs_set_is_defined;
  mutable std::multiset<variable> m_variables_in_rhs;

public:

  /// \brief Type of variables
  using variable_type = VariableType;

  /// \brief Type of expressions
  using expression_type = ExpressionType;

  using argument_type = variable_type;
  using result_type = expression_type;

  /// \brief Default constructor
  mutable_indexed_substitution()
    : m_variables_in_rhs_set_is_defined(false)
  {
  }

  mutable_indexed_substitution(const substitution_type& substitution,
                               const bool variables_in_rhs_set_is_defined,
                               const std::multiset<variable_type>& variables_in_rhs)
      : m_substitution(substitution), 
        m_variables_in_rhs_set_is_defined(variables_in_rhs_set_is_defined),
        m_variables_in_rhs(variables_in_rhs)
  {
  }

  /// \brief Wrapper class for internal storage and substitution updates using operator()
  struct assignment
  {
    const variable_type& m_variable;
    mutable_indexed_substitution < VariableType, ExpressionType >& m_super;

    /// \brief Constructor.
    /// \param[in] v a variable.
    /// \param[in] super A reference to the surrounding indexed substitution.
    assignment(const variable_type& v, 
               mutable_indexed_substitution < VariableType, ExpressionType >& super)
     : m_variable(v),
       m_super(super)
    { }

    /// \brief Actual assignment
    void operator=(const expression_type& e)
    {
      assert(e.defined());
      const typename substitution_type::iterator i = m_super.m_substitution.find(m_variable);
      if (i!=m_super.m_substitution.end())
      {
        // Found.
        assert(i->first==m_variable);
        if (e==i->second)  // No change in the substitution is required. 
        {
          return;
        }
        if (m_super.m_variables_in_rhs_set_is_defined)
        {
          std::set<variable_type> s=find_free_variables(i->second);
          for(const variable& v: s) 
          { 
            // Remove one occurrence of v. 
            const typename std::multiset<variable_type>::const_iterator j=m_super.m_variables_in_rhs.find(v);
            if (j!=m_super.m_variables_in_rhs.end())
            { 
              m_super.m_variables_in_rhs.erase(j);
            }
          }
        }
        if (e != m_variable)
        {
            i->second=e;
        }
        else
        {
          m_super.m_substitution.erase(i);
        }
      }
      else
      {
        // Not found.
        if (e!=m_variable)
        {
          m_super.m_substitution.emplace(m_variable,e);
        }
      }
     
      if (e != m_variable && m_super.m_variables_in_rhs_set_is_defined)
      {
        std::set<variable_type> s1=find_free_variables(e);
        m_super.m_variables_in_rhs.insert(s1.begin(),s1.end());
      }
    }
  };

  /// \brief Application operator; applies substitution to v.
  /// \param   v The variable to which the subsitution is applied.
  /// \result  The value to which v is mapped, or v itself if v is not
  ///          mapped to a value. 
  const expression_type& operator()(const variable_type& v) const
  {
    typename substitution_type::const_iterator i=m_substitution.find(v);
    if (i==m_substitution.end())  // not found.
    {
      return v;
    }
    // no value assigned to v;
    assert(i->first==v); 
    return i->second; 
  }

  /// \brief Application operator; applies substitution to v.
  /// \details This must deliver an expression, and not a reference
  ///          to an expression, as the expressions are stored in 
  ///          a vector that can be resized and moved. 
  /// \param   v The variable to which the subsitution is applied.
  /// \param   target The target into which the substitution is stored. 
  template <class ResultType>
  void apply(const variable_type& v, ResultType& target)
  {
    static_assert(
        std::is_same_v<ResultType&, expression_type&> || std::is_same_v<ResultType&, atermpp::unprotected_aterm_core&>);
    const typename substitution_type::const_iterator i=m_substitution.find(v);
    if (i==m_substitution.end()) // not found.
    {
      target=v;
      return;
    }
    target=i->second;
  }

  /// \brief Application operator; applies substitution to v.
  /// \details This must deliver an expression, and not a reference
  ///          to an expression, as the expressions are stored in 
  ///          a vector that can be resized and moved. 
  /// \param   v The variable to which the subsitution is applied.
  /// \param   target The target into which the substitution is stored. 
  void apply(const variable_type& v, 
             expression_type& target,
             atermpp::detail::thread_aterm_pool& thread_aterm_pool)
  {
    const typename substitution_type::iterator i=m_substitution.find(v);
    if (i==m_substitution.end()) // not found.
    {
      target.assign(v, thread_aterm_pool);
      return;
    }
    target.assign(i->second, thread_aterm_pool);
  }

  /// \brief Index operator.
  assignment operator[](variable_type const& v)
  {
    return assignment(v, *this);
  }

  /// \brief Clear substitutions.
  void clear()
  {
    m_substitution.clear();
    m_variables_in_rhs_set_is_defined=false;
    m_variables_in_rhs.clear();
  }

  /// \brief Create a clone of the rewriter in which the underlying rewriter is
  /// copied, and not passed as a shared pointer.
  /// \details This is useful when the rewriter is used in different parallel
  /// processes. One rewriter can only be used sequentially. \return A rewriter,
  /// with a copy of the underlying jitty, jittyc or jittyp rewriting engine.
  mutable_indexed_substitution<VariableType, ExpressionType> clone()
  {
    return mutable_indexed_substitution<VariableType, ExpressionType>(
                       m_substitution, m_variables_in_rhs_set_is_defined, m_variables_in_rhs);
  }

  /// \brief Compare substitutions
  template <typename Substitution>
  bool operator==(const Substitution&) const
  {
    return false;
  }

  /// \brief Provides a multiset containing the variables in the rhs.
  /// \return A multiset with variables in the right hand side. 
  const std::multiset<variable_type>& variables_occurring_in_right_hand_sides() const
  {
    if (!m_variables_in_rhs_set_is_defined)
    {
      for(const auto& p: m_substitution)
      {
        std::set<variable_type> s=find_free_variables(p.second);
        m_variables_in_rhs.insert(s.begin(),s.end());
      }
      m_variables_in_rhs_set_is_defined=true;
    }
    return m_variables_in_rhs;
  }

  /// \brief Checks whether a variable v occurs in one of the rhs's of the current substitution.
  /// \return A boolean indicating whether v occurs in sigma(x) for some variable x. 
  bool variable_occurs_in_a_rhs(const variable& v)
  {
    const std::multiset<variable>& variables_in_rhs=variables_occurring_in_right_hand_sides();
    return variables_in_rhs.find(v)!=variables_in_rhs.end();
  }

  /// \brief Returns the number of assigned variables in the substitution.
  bool size()
  {
    return m_substitution.size();
  }

  /// \brief Returns true if the substitution is empty.
  bool empty()
  {
    return m_substitution.empty();
  }

public:
  /// \brief string representation of the substitution. N.B. This is an expensive operation!
  std::string to_string() const
  {
    std::stringstream result;
    bool first = true;
    result << "[";
    for (const auto& p: m_substitution)
    {
      if (first)
      {
        first = false;
      }
      else
      {
        result << "; ";
      }
        
      result << p.first << " := " << p.second;
    }
    result << "]";
    return result.str();
  }

};


template <typename VariableType, typename ExpressionType>
std::ostream& operator<<(std::ostream& out, const mutable_indexed_substitution<VariableType, ExpressionType>& sigma)
{
  return out << sigma.to_string();
}

template <typename VariableType = variable, typename ExpressionType = data_expression>
std::multiset<variable> substitution_variables(const mutable_indexed_substitution<VariableType, ExpressionType>& sigma)
{
  return sigma.variables_occurring_in_right_hand_sides();
}

} // namespace mcrl2::data



#endif // MCRL2_DATA_SUBSTITUTIONS_MUTABLE_INDEXED_SUBSTITUTION_H
