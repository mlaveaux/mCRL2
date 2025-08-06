// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/add_binding.h
/// \brief add your file description here.

#ifndef MCRL2_CORE_ADD_BINDING_H
#define MCRL2_CORE_ADD_BINDING_H

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/type_traits.cxx"
#else
  import atermpp;
#endif

namespace mcrl2::core
{

/// \brief Traverser that defines functions for maintaining bound variables.
template <template <class> class TraverserOrBuilder, typename Derived, typename Variable>
class add_binding : public TraverserOrBuilder<Derived>
{
  public:
    using super = TraverserOrBuilder<Derived>;
    using variable_type = Variable;

  protected:
    std::multiset<variable_type> m_bound_variables;

    /// \brief Add a variable to the multiset of bound variables.
    void increase_bind_count(const variable_type& var)
    {
      m_bound_variables.insert(var);
    }

    /// \brief Add a sequence of variables to the multiset of bound variables.
    template <typename Container>
    void increase_bind_count(const Container& variables, typename atermpp::enable_if_container<Container, variable_type>::type* = nullptr)
    {
      for (typename Container::const_iterator i = variables.begin(); i != variables.end(); ++i)
      {
        m_bound_variables.insert(*i);
      }
    }

    /// \brief Remove a variable from the multiset of bound variables.
    void decrease_bind_count(const variable_type& var)
    {
      m_bound_variables.erase(m_bound_variables.find(var));
    }

    /// \brief Remove a sequence of variables from the multiset of bound variables.
    template <typename Container>
    void decrease_bind_count(const Container& variables, typename atermpp::enable_if_container<Container, variable_type>::type* = nullptr)
    {
      for (typename Container::const_iterator i = variables.begin(); i != variables.end(); ++i)
      {
        m_bound_variables.erase(m_bound_variables.find(*i));
      }
    }

  public:
    /// \brief Returns true if the variable v is bound.
    bool is_bound(variable_type const& v) const
    {
      return m_bound_variables.find(v) != m_bound_variables.end();
    }

    /// \brief Returns the bound variables.
    const std::multiset<variable_type>& bound_variables() const
    {
      return m_bound_variables;
    }

    /// \brief Returns the bind count of the variable v.
    typename std::multiset<variable_type>::size_type bind_count(const variable_type& v)
    {
      return m_bound_variables.count(v);
    }
};

} // namespace mcrl2::core

#endif // MCRL2_CORE_ADD_BINDING_H
