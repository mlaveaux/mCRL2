// Author(s): Jan Friso Groote. Based on pbes/rewriters/enumerate_quantifiers_rewriter.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/rewriters/enumerate_quantifiers_rewriter.h
/// \brief A rewriter that enumerates infimum, supremum and sum operators in pres expressions. 

#ifndef MCRL2_PRES_REWRITERS_ENUMERATE_QUANTIFIERS_REWRITER_H
#define MCRL2_PRES_REWRITERS_ENUMERATE_QUANTIFIERS_REWRITER_H

#include <numeric>
#include "mcrl2/atermpp/standard_containers/vector.h"
#include "mcrl2/data/detail/split_finite_variables.h"
#include "mcrl2/data/real_utilities.h"
#include "mcrl2/data/enumerator.h"
#include "mcrl2/data/rewriters/enumerate_quantifiers_rewriter.h"
#include "mcrl2/pbes/rewriters/quantifier_expansion_mode.h"
#include "mcrl2/pres/rewriters/simplify_rewriter.h"



namespace mcrl2::pres_system {

namespace detail {

// Simplifying PRES rewriter that eliminates quantifiers using enumeration.
/// \param MutableSubstitution This must be a MapSubstitution.
template <typename Derived, typename DataRewriter, typename MutableSubstitution>
struct enumerate_quantifiers_builder: public simplify_data_rewriter_builder<Derived, DataRewriter, MutableSubstitution>
{
  using super = simplify_data_rewriter_builder<Derived, DataRewriter, MutableSubstitution>;
  using self = enumerate_quantifiers_builder<Derived, DataRewriter, MutableSubstitution>;
  using enumerator_element = data::enumerator_list_element<pres_expression>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;
  using super::sigma;
  using super::R;

  const data::data_specification& m_dataspec;

  /// If true, quantifier variables of infinite sort are enumerated.
  const pbes_system::enumerate_quantifiers_mode m_enum_mode; 

  /// identifier generator
  data::enumerator_identifier_generator& m_id_generator;

  /// The enumerator
  data::enumerator_algorithm<self> E;

  /// \brief Constructor.
  /// \param R A data rewriter.
  /// \param sigma A mutable substitution.
  /// \param dataspec A data specification.
  /// \param id_generator A generator to generate fresh variable names.
  /// \param enum_mode Indicates how quantifiers, including inf and sup are expanded. 
  enumerate_quantifiers_builder(const DataRewriter& R,
                                MutableSubstitution& sigma,
                                const data::data_specification& dataspec,
                                data::enumerator_identifier_generator& id_generator,
                                const pbes_system::enumerate_quantifiers_mode enum_mode = pbes_system::expand_infinite_sorts_and_use_data_rewriter)
    : super(dataspec, R, sigma), 
      m_dataspec(dataspec), 
      m_enum_mode(enum_mode), 
      m_id_generator(id_generator),
      E(*this, m_dataspec, R, id_generator, (std::numeric_limits<std::size_t>::max)())
  { }

  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }

  atermpp::vector<data::data_expression> undo_substitution(const data::variable_list& variables)
  {
    atermpp::vector<data::data_expression> result;
    for (const data::variable& v : variables)
    {
      result.push_back(sigma(v));
      sigma[v] = v;
    }
    return result;
  }

  void redo_substitution(const data::variable_list& v, const atermpp::vector<data::data_expression>& undo)
  {
    assert(v.size() == undo.size());
    auto i = v.begin();
    auto j = undo.begin();
    while (i != v.end())
    {
      sigma[*i++] = *j++;
    }
  }

  // We assume that phi is already rewritten. 
  void enumerate_infimum(pres_expression& result, const data::variable_list& v, const pres_expression& phi)
  {   
    assert(!v.empty());
    assert(&result!=&phi);
    atermpp::vector<data::data_expression> undo = undo_substitution(v);
    result = true_();
#ifndef NDEBUG
    pres_expression phi_;
    derived().apply(phi_, phi);
    assert(phi_==phi); // phi is assumed to be rewritten. 
#endif
    E.enumerate(enumerator_element(v, phi),
                sigma,
                [&](const enumerator_element& p)
                {
                  pres_system::make_optimized_and(result, result, p.expression());
                  return is_false(result);
                },
                is_true,
                is_false
               );
      
    redo_substitution(v, undo);
  }                        

  void enumerate_supremum(pres_expression& result, const data::variable_list& v, const pres_expression& phi)
  {
    assert(!v.empty());
    assert(&result!=&phi);
    atermpp::vector<data::data_expression> undo = undo_substitution(v);
    result = false_();
#ifndef NDEBUG
    pres_expression phi_; 
    derived().apply(phi_, phi);
    assert(phi_==phi); // phi is assumed to be rewritten. 
#endif
    E.enumerate(enumerator_element(v, phi),
                sigma,
                [&](const enumerator_element& p)
                {
                  pres_system::make_optimized_or(result, result, p.expression());
                  return is_true(result);
                },
                is_false,
                is_true
    );

    redo_substitution(v, undo);
  }

  void enumerate_sum(pres_expression& result, const data::variable_list& v, const pres_expression& phi)
  {
    assert(!v.empty());
    assert(&result!=&phi);
    atermpp::vector<data::data_expression> undo = undo_substitution(v);
    result = atermpp::down_cast<pres_expression>(data::sort_real::real_zero());
#ifndef NDEBUG
    pres_expression phi_; 
    derived().apply(phi_, phi);
    assert(phi_==phi); // phi is assumed to be rewritten. 
#endif
    E.enumerate(enumerator_element(v, phi),
                sigma,
                [&](const enumerator_element& p)
                {
                  make_optimized_plus(result, result, p.expression());
                  return is_true(result);
                },
                data::sort_real::is_zero,  
                is_true
    );

    redo_substitution(v, undo);
  }

  template <class T>
  void apply(T& result, const infimum& x)
  {
    atermpp::vector<data::data_expression> undo = undo_substitution(x.variables());
    derived().apply(result, x.body());
    std::set<data::variable> free_variables = find_free_variables(result);
    if (m_enum_mode==pbes_system::expand_finite_sorts)
    {
      data::variable_list finite;
      data::variable_list infinite;
      data::variable_list unused;
      data::detail::split_finite_variables(x.variables(), m_dataspec,
                                           finite, infinite, unused,
                                           [&free_variables](const data::variable& v){ return free_variables.count(v)>0; });
      if (finite.empty())
      { 
        make_optimized_infimum(result, infinite, result);
      }
      else
      {
        pres_expression phi_;
        enumerate_infimum(phi_, finite, result);
        make_optimized_infimum(result, infinite, phi_);
      }
    }
    else
    {
      data::variable_list enumerable;
      data::variable_list non_enumerable;
      data::variable_list unused;
      data::detail::split_enumerable_variables(x.variables(), m_dataspec, super::R,
                                               enumerable, non_enumerable, unused,
                                               [&free_variables](const data::variable& v){ return free_variables.count(v)>0; });
      if (enumerable.empty())
      {
        make_optimized_infimum(result, non_enumerable, result);
      }
      else
      {
        pres_expression phi_;
        enumerate_infimum(phi_, enumerable, result);
        make_optimized_infimum(result, non_enumerable, phi_);
      }
    }
    redo_substitution(x.variables(), undo);
  }

  template <class T>
  void apply(T& result, const supremum& x)
  {
    atermpp::vector<data::data_expression> undo = undo_substitution(x.variables());
    derived().apply(result, x.body());
    std::set<data::variable> free_variables = find_free_variables(result);
    if (m_enum_mode==pbes_system::expand_finite_sorts)
    {
      data::variable_list finite;
      data::variable_list infinite;
      data::variable_list unused;
      data::detail::split_finite_variables(x.variables(), m_dataspec,
                                           finite, infinite, unused,
                                           [&free_variables](const data::variable& v){ return free_variables.count(v)>0; });
      if (finite.empty())
      { 
        make_optimized_supremum(result, infinite, result);
      }
      else
      {
        pres_expression phi_;
        enumerate_supremum(phi_, finite, result);
        make_optimized_supremum(result, infinite, phi_);
      }
    }
    else
    {
      data::variable_list enumerable;
      data::variable_list non_enumerable;
      data::variable_list unused;
      data::detail::split_enumerable_variables(x.variables(), m_dataspec, super::R,
                                               enumerable, non_enumerable, unused,
                                               [&free_variables](const data::variable& v){ return free_variables.count(v)>0; });
      if (enumerable.empty())
      {
        make_optimized_supremum(result, non_enumerable, result);
      }
      else
      {
        pres_expression phi_;
        enumerate_supremum(phi_, enumerable, result);
        make_optimized_supremum(result, non_enumerable, phi_);
      }
    }
    redo_substitution(x.variables(), undo);
  }

  template <class T>
  void apply(T& result, const sum& x)
  {
    atermpp::vector<data::data_expression> undo = undo_substitution(x.variables());
    derived().apply(result, x.body());
    std::set<data::variable> free_variables = find_free_variables(result);
    if (m_enum_mode==pbes_system::expand_finite_sorts)
    {
      data::variable_list finite;
      data::variable_list infinite;
      data::variable_list unused;
      data::detail::split_finite_variables(x.variables(), m_dataspec,
                                           finite, infinite, unused,
                                           [&free_variables](const data::variable& v){ return free_variables.count(v)>0; });
      if (finite.empty())
      { 
        make_optimized_sum(result, infinite + unused, result, m_dataspec, super::R);
      }
      else
      {
        pres_expression phi_;
        enumerate_sum(phi_, finite, result);
        make_optimized_sum(result, infinite + unused, phi_, m_dataspec, super::R);
      }
    }
    else
    {
      data::variable_list enumerable;
      data::variable_list non_enumerable;
      data::variable_list unused;
      data::detail::split_enumerable_variables(x.variables(), m_dataspec, super::R,
                                               enumerable, non_enumerable, unused,
                                               [&free_variables](const data::variable& v){ return free_variables.count(v)>0; });
      if (enumerable.empty())
      {
        make_optimized_sum(result, non_enumerable + unused, result, m_dataspec, super::R);
      }
      else
      {
        pres_expression phi_;
        enumerate_sum(phi_, enumerable, result);
        make_optimized_sum(result, non_enumerable + unused, phi_, m_dataspec, super::R);
      }
    }
    redo_substitution(x.variables(), undo);
  }

  template <atermpp::IsATerm T>
  void apply(T& result, const data::data_expression& x)
  { 
    if (m_enum_mode==pbes_system::expand_infinite_sorts_and_use_data_rewriter)
    {
      R(atermpp::assign_cast<data::data_expression>(result),x,sigma);
    }
    else
    {
      data::data_expression temporary_result;
      data::enumerate_quantifiers_rewriter(R,m_dataspec,m_id_generator,m_enum_mode==pbes_system::expand_infinite_sorts)(temporary_result, x, sigma);
      R(atermpp::assign_cast<data::data_expression>(result),temporary_result,sigma);
    }
  }
 
  // N.B. This function has been added to make this class operate well with the enumerator.
  pres_expression operator()(const pres_expression& x, MutableSubstitution& )
  {
    pres_expression result;
    derived().apply(result, x);
    return result;
  }

  // N.B. This function has been added to make this class operate well with the enumerator.
  //      that employs the "make_rewrite" variant with an explicit result. 
  void operator()(pres_expression& result, const pres_expression& x, MutableSubstitution& )
  {
    derived().apply(result, x);
  }
};

template <template <class, class, class> class Builder, class DataRewriter, class MutableSubstitution>
struct apply_enumerate_builder: public Builder<apply_enumerate_builder<Builder, DataRewriter, MutableSubstitution>, DataRewriter, MutableSubstitution>
{
  using super
      = Builder<apply_enumerate_builder<Builder, DataRewriter, MutableSubstitution>, DataRewriter, MutableSubstitution>;
  using super::enter;
  using super::leave;

  apply_enumerate_builder(const DataRewriter& R, MutableSubstitution& sigma, const data::data_specification& dataspec, data::enumerator_identifier_generator& id_generator, const pbes_system::enumerate_quantifiers_mode enum_mode)
    : super(R, sigma, dataspec, id_generator, enum_mode)
  {}
}; 

template <template <class, class, class> class Builder, class DataRewriter, class MutableSubstitution>
apply_enumerate_builder<Builder, DataRewriter, MutableSubstitution>
make_apply_enumerate_builder(const DataRewriter& R, MutableSubstitution& sigma, const data::data_specification& dataspec, data::enumerator_identifier_generator& id_generator, const pbes_system::enumerate_quantifiers_mode enum_mode)
{
  return apply_enumerate_builder<Builder, DataRewriter, MutableSubstitution>(R, sigma, dataspec, id_generator, enum_mode);
}

} // namespace detail

/// \brief An attempt for improving the efficiency.
struct enumerate_quantifiers_rewriter
{
  protected:
    /// \brief A data rewriter
    data::rewriter m_rewriter;

    /// \brief A data specification
    data::data_specification m_dataspec;

    /// \brief If true, quantifier variables of infinite sort are enumerated.
    const pbes_system::enumerate_quantifiers_mode m_enum_mode;

    mutable data::enumerator_identifier_generator m_id_generator;

  public:
    using term_type = pres_expression;
    using variable_type = data::variable;

    enumerate_quantifiers_rewriter(const data::rewriter& R, const data::data_specification& dataspec, const pbes_system::enumerate_quantifiers_mode enum_mode = pbes_system::expand_infinite_sorts_and_use_data_rewriter)
      : m_rewriter(R), m_dataspec(dataspec), m_enum_mode(enum_mode)
    {}

    pres_expression operator()(const pres_expression& x) const
    {
      data::rewriter::substitution_type sigma;
      pres_expression result;
      detail::apply_enumerate_builder<detail::enumerate_quantifiers_builder, data::rewriter, data::rewriter::substitution_type>(m_rewriter, sigma, m_dataspec, m_id_generator, m_enum_mode).apply(result, x);
      return result;
    }

    template <typename MutableSubstitution>
    pres_expression operator()(const pres_expression& x, MutableSubstitution& sigma) const
    {
      pres_expression result;
      detail::apply_enumerate_builder<detail::enumerate_quantifiers_builder, data::rewriter, MutableSubstitution>(m_rewriter, sigma, m_dataspec, m_id_generator, m_enum_mode).apply(result, x);
      return result;
    }

    template <typename MutableSubstitution>
    void operator()(pres_expression& result, const pres_expression& x, MutableSubstitution& sigma) const
    {
      detail::apply_enumerate_builder<detail::enumerate_quantifiers_builder, data::rewriter, MutableSubstitution>(m_rewriter, sigma, m_dataspec, m_id_generator, m_enum_mode).apply(result, x);
    }

    void clear_identifier_generator()
    {
      m_id_generator.clear();
    }

    /// \brief Create a clone of the rewriter in which the underlying rewriter is copied, and not passed as a shared pointer. 
    /// \details This is useful when the rewriter is used in different parallel processes. One rewriter can only be used sequentially. 
    /// \return A rewriter, with a copy of the underlying jitty, jittyc or jittyp rewriting engine. 
    enumerate_quantifiers_rewriter clone()
    {
      return enumerate_quantifiers_rewriter(m_rewriter.clone(), m_dataspec, m_enum_mode);
    }

    /// \brief Initialises this rewriter with thread dependent information. 
    /// \details This function sets a pointer to the m_busy_flag and m_forbidden_flag of this.
    ///          process, such that rewriter can access these faster than via the general thread.
    ///          local mechanism. It is expected that this is not needed when compilers become.
    ///          faster, and should be removed in due time. 
    void thread_initialise()
    {
      m_rewriter.thread_initialise();
    }
};

} // namespace mcrl2::pres_system



#endif // MCRL2_PRES_REWRITERS_ENUMERATE_QUANTIFIERS_REWRITER_H

