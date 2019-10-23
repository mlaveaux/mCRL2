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

#include "mcrl2/data/abstraction.h"
#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/variable.h"
#include "mcrl2/data/detail/rewrite/substitute.h"
#include "mcrl2/utilities/stack_array.h"

#include <assert.h>

namespace mcrl2
{
namespace data
{
namespace detail
{

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
  /// \brief Default constructor.
  ConstructionStack() {}

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
  template<typename Substitution,
    typename Generator>
  data_expression construct_term(Substitution& sigma, Generator& generator, std::vector<data_expression>& argument_stack) const
  {
    // Define an iterative version to prevent a large number of nested calls before returning the value.
    for (const atermpp::aterm_appl& term : m_stack)
    {
      if (is_variable(term))
      {
        // e(x |> Q, S, sigma) = e(Q, S |> sigma(x), sigma)
        const auto& var = static_cast<const variable&>(term);
        argument_stack.push_back(sigma(var));
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
      else if (is_abstraction(term))
      {
        // e(lambda x . t |> Q, S, sigma) = e(Q, S |> (lambda x . t)^sigma, sigma)
        argument_stack.push_back(capture_avoiding_substitution(static_cast<const data_expression>(term), sigma, generator));
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
      m_stack.push_back(term);
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
