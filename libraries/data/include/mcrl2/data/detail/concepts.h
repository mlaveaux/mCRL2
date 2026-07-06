#ifndef MCRL2_DATA_DETAIL_CONCEPTS_H
#define MCRL2_DATA_DETAIL_CONCEPTS_H

#include <boost/concept_check.hpp>
#include <boost/iterator/iterator_concepts.hpp>

#include <boost/concept/detail/concept_def.hpp>

namespace mcrl2::data::concepts
{

/// \cond INTERNAL_DOCS

BOOST_concept(Substitution,(S)) :
  boost::UnaryFunction< S, typename S::expression_type, typename S::expression_type >
{
  using variable_type = typename S::variable_type;
  using expression_type = typename S::expression_type;

  variable_type    v;
  expression_type  e;

  S s;

  BOOST_CONCEPT_USAGE(Substitution)
  {
    static_assert(std::is_convertible_v<variable_type, expression_type>,
        "variable_type and expression_type must be convertible");

    s(static_cast< expression_type const& >(e)) == e;
    s(static_cast< variable_type const& >(v)) == v;
  }
};

BOOST_concept(MutableSubstitution,(S)) :
  boost::DefaultConstructible< S >,
  Substitution< S >
{
  using variable_type = typename S::variable_type;
  using expression_type = typename S::expression_type;

  variable_type    v;
  expression_type  e;

  S s;

  BOOST_CONCEPT_USAGE(MutableSubstitution)
  {
    s[static_cast< variable_type const& >(v)] = e;
    s[static_cast< variable_type const& >(v)] = v;
  }
};


BOOST_concept(Evaluator,(C)(S))
{

  BOOST_CONCEPT_ASSERT((Substitution< S >));

  C                           c;
  S                           s;
  typename C::expression_type e;

  BOOST_CONCEPT_USAGE(Evaluator)
  {
    e = c(e);
    e = c(e, s);
  }
};

BOOST_concept(Enumerator,(E)) : boost_concepts::ForwardTraversal< E >
{

  BOOST_CONCEPT_ASSERT((Substitution< typename E::value_type >));
  BOOST_CONCEPT_ASSERT((Evaluator< typename E::evaluator_type, typename E::substitution_type >));

  typename E::value_type   s;
  E                        e;

  BOOST_CONCEPT_USAGE(Enumerator)
  {

    ++e;
    s = *e;
  }
};
} // namespace mcrl2::data::concepts

/// \endcond



#include <boost/concept/detail/concept_undef.hpp>

#endif
