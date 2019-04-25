// Author(s): Muck van Weerdenburg
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_REWRITE_JITTY_H
#define MCRL2_DATA_DETAIL_REWRITE_JITTY_H

#include "mcrl2/data/detail/rewrite.h"
#include "mcrl2/data/detail/rewrite/strategy_rule.h"

namespace mcrl2
{
namespace data
{
namespace detail
{

class RewriterJitty: public Rewriter
{
  public:
    typedef Rewriter::substitution_type substitution_type;

    RewriterJitty(const data_specification& data_spec, const used_data_equation_selector &);
    virtual ~RewriterJitty();

    rewrite_strategy getStrategy();

    data_expression rewrite(const data_expression &term, substitution_type &sigma);

  private:
    std::map< function_symbol, data_equation_list > jitty_eqns;
    std::vector<strategy> jitty_strat;

    data_expression rewrite_aux(const data_expression& term, substitution_type& sigma);

    data_expression rewrite_aux_function_symbol(
                      const function_symbol& op,
                      const data_expression& term,
                      substitution_type& sigma);

    data_expression rewrite_aux_const_function_symbol(
                      const function_symbol& op,
                      substitution_type& sigma);

    /// \brief Auxiliary function to take care that the array jitty_strat is sufficiently large
    ///        to access element i.
    void make_jitty_strat_sufficiently_larger(const std::size_t i);
    void rebuild_strategy();
};

/// \brief removes auxiliary expressions this_term_is_in_normal_form from data_expressions that are being rewritten.
/// \details The function below is intended to remove the auxiliary function this_term_is_in_normal_form from a term
///          such that it can for instance be pretty printed. This auxiliary function is used internally in terms
///          when rewriting to avoid to rewrite too often.
data_expression remove_normal_form_function(const data_expression& t);

}
}
}

#endif
