// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_REWRITE_UTILITY_H
#define MCRL2_DATA_DETAIL_REWRITE_UTILITY_H

#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/variable.h"
#include "mcrl2/utilities/stack_array.h"

#include <assert.h>

namespace mcrl2
{
namespace data
{
namespace detail
{

/// \brief A capture-avoiding substitution of sigma applied to the given term.
template<typename Substitution>
static inline data_expression capture_avoiding_substitution(const data_expression& term, const Substitution& sigma)
{
  // C(x, sigma, V) = sigma(x), where x is a variable
  if (is_variable(term))
  {
    const auto& var = static_cast<const variable&>(term);
    return sigma.at(var);
  }
  // C(f, sigma, V) = f, where f is a function symbol.
  else if (is_function_symbol(term))
  {
    return term;
  }
  // C(lambda x . t, sigma, V) = lambda y . C(t, sigma[x := y], V), where x and y are variables.
  else if (is_abstraction(term))
  {
    //const auto& abstraction = static_cast<const class abstraction&>(term);
    assert(false);
    return term;
  }
  // C(t(t_1, ..., t_n, sigma, V) = C(t, sigma, V) ( C(t_1, sigma, V), ..., C(t_n, sigma, V) )
  else
  {
    assert(is_application(term));
    const auto& appl = static_cast<const application&>(term);

    // Substitution of all arguments.
    MCRL2_DECLARE_STACK_ARRAY(arguments, data_expression, appl.size());
    for (std::size_t index = 0; index < appl.size(); ++index)
    {
      arguments[index] = capture_avoiding_substitution(appl[index], sigma);
    }

    // Construct the application, also subsituting the head.
    return application(capture_avoiding_substitution(appl.head(), sigma), arguments.begin(), arguments.end());
  }
}

/// \brief Matches a single left-hand side with the given term and creates the substitution.
template<typename Substitution>
static inline bool match_lhs(const data_expression& term,  const data_expression& lhs, Substitution& sigma)
{
  if (term == lhs)
  {
    // If both sides are equivalent, then they match under any substitution. If one contains
    // a variable then that variable is also not bound in the other term.
    return true;
  }
  else if (is_function_symbol(lhs))
  {
    return term == lhs;
  }
  else if (is_variable(lhs))
  {
    const auto& var = static_cast<const variable&>(lhs);

    if (sigma.count(var) != 0)
    {
      // If the variable was already assigned it must match the previously assigned value.
      return sigma.at(var) == term;
    }
    else
    {
      // Else substitute the current term for this variable.
      sigma[var] = term;
      return true;
    }
  }
  else
  {
    // The term and lhs are applications.
    const application& lhs_appl  = static_cast<const application&>(lhs);
    const application& term_appl = static_cast<const application&>(term);

    // Both must have the same arity, the head symbol must match and their arguments must match.
    if (lhs_appl.size() != term_appl.size())
    {
      return false;
    }

    if (!match_lhs(term_appl.head(), lhs_appl.head(), sigma))
    {
      return false;
    }

    for (std::size_t i = 0; i < term_appl.size(); i++)
    {
      if (!match_lhs(term_appl[i], lhs_appl[i], sigma))
      {
        return false;
      }
    }

    return true;
  }
}

/// \brief A construction stack is a data structure that can be initialized from a given term to quickly reconstruct this term
///        bottom-up while applied a given substitution. A typical usecase is the right-hand side of a rewrite rule. Here, the structure
///        of the term is known beforehand.
class ConstructionStack
{
private:
  /// \brief A function symbol with arity (data_expression and an arity)
  class function_symbol_arity : public atermpp::aterm_appl
  {
  public:
    static atermpp::function_symbol& g_function_symbol_arity()
    {
      static atermpp::function_symbol constant("@@function_symbol_arity@@", 2);
      return constant;
    }

    function_symbol_arity(const data_expression& expression, std::size_t arity)
      : atermpp::aterm_appl(g_function_symbol_arity(), expression, atermpp::aterm_int(arity))
    {}

    const data_expression& head() const { return static_cast<const data_expression&>(this->operator[](0)); }

    std::size_t arity() const { return static_cast<const atermpp::aterm_int&>(this->operator[](1)).value(); }
  };

  /// \returns True when the given term is of type function_symbol_arity.
  static inline bool is_function_symbol_arity(const atermpp::aterm& term)
  {
    return term.function() == function_symbol_arity::g_function_symbol_arity();
  }

public:
  /// \brief Constructs a stack that can be used to reconstruct the given term without traversing it top-down.
  ConstructionStack(const data_expression& term)
  {
    build_construction_stack_impl(term);
  }

  /// \returns The evaluation of the construction stack, defined by e(Q, S, sigma), equivalent to capture_avoiding_substitution(term, sigma); where
  ///          term is the term that was used for construction of this object.
  template<typename Substitution>
  data_expression construct_term(const Substitution& sigma) const
  {
    std::vector<data_expression> argument_stack;
    return construct_term(sigma, argument_stack);
  }

  /// \brief The same function as above, but reuses an existing argument stack to prevent reallocations. This "stack" is a vector
  ///        to be able to iterate over the arity number of arguments directly.
  template<typename Substitution>
  data_expression construct_term(const Substitution& sigma, std::vector<data_expression>& argument_stack) const
  {
    // Define an iterative version to prevent a large number of nested calls before returning the value.
    for (const auto& term : m_stack)
    {
      if (is_variable(term))
      {
        // e(x |> Q, S, sigma) = e(Q, S |> sigma(x), sigma)
        const auto& var = static_cast<const variable&>(term);
        argument_stack.push_back(sigma.at(var));
      }
      else if (is_function_symbol_arity(term))
      {
        // e(f |> Q, S |> t_0 |> ... |> t_{arity(f)}, sigma) = e(Q, S |> MATCH_APPLY(f(t_0, ..., t__{arity(f)}), sigma)
        const auto& symbol = static_cast<const function_symbol_arity&>(term);
        data_expression result = application(symbol.head(), argument_stack.end() - static_cast<long>(symbol.arity()), argument_stack.end());

        // Remove arity(f) number of arguments from the stack
        argument_stack.erase(argument_stack.end() - static_cast<long>(symbol.arity()), argument_stack.end());
        argument_stack.push_back(result);
      }
      else
      {
        // e(t |> Q, S, sigma) = e(Q, S |> t, sigma)
        argument_stack.push_back(static_cast<const data_expression>(term));
      }
    }

    // e([], t |> [], sigma) = t
    assert(argument_stack.size() == 1);

    // Clear the argument stack to be able to reuse it.
    data_expression term = argument_stack.back();
    argument_stack.clear();
    return term;
  }

private:

  /// \brief Updates the stack directly, to prevent unnecessary copying and stack concatenation.
  void build_construction_stack_impl(const data_expression& term)
  {
    if (is_variable(term))
    {
      // c(x) = x
      m_stack.push_back(term);
    }
    else if (is_function_symbol(term))
    {
      // c(f) = f   as FV(f) = empty (normally handled by the f(t_0, ..., t_n) case)
      m_stack.push_back(term);
    }
    else if (is_abstraction(term))
    {
      // Ignored for now.
    }
    else
    {
      assert(is_application(term));
      const auto& appl = static_cast<const application&>(term);

      if (find_free_variables(term).empty())
      {
        // c(f(t_0, ..., t_n)) = f(t_0, ..., t_n)  if FV(f(t_0, ..., t_n)) = empty
        m_stack.push_back(term);
      }
      else
      {
        // c(f(t_0, ..., t_n)) = c(t_0) |> ... |> c(t_n) |> f  if FV(f(t_0, ..., t_n)) != empty
        for (auto& argument : appl)
        {
          build_construction_stack_impl(argument);
        }

        m_stack.push_back(function_symbol_arity(appl.head(), appl.size()));
      }
    }
  }

  std::vector<atermpp::aterm_appl> m_stack;
};

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_REWRITE_UTILITY_H
