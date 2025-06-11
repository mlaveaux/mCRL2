// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/data_expression.h
/// \brief The class data_expression.

#ifndef MCRL2_DATA_DATA_EXPRESSION_H
#define MCRL2_DATA_DATA_EXPRESSION_H

#include "mcrl2/data/container_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/untyped_sort.h"

namespace mcrl2
{

namespace data
{

/// \brief Returns true if the term t is an abstraction
inline bool is_abstraction(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Binder;
}

/// \brief Returns true if the term t is a lambda abstraction
inline bool is_lambda(const atermpp::aterm& x)
{
  return is_abstraction(x) && x[0].function() == core::detail::function_symbols::Lambda;
}

/// \brief Returns true if the term t is a universal quantification
inline bool is_forall(const atermpp::aterm& x)
{
  return is_abstraction(x) && x[0].function() == core::detail::function_symbols::Forall;
}

/// \brief Returns true if the term t is an existential quantification
inline bool is_exists(const atermpp::aterm& x)
{
  return is_abstraction(x) && x[0].function() == core::detail::function_symbols::Exists;
}

/// \brief Returns true if the term t is a set comprehension
inline bool is_set_comprehension(const atermpp::aterm& x)
{
  return is_abstraction(x) && x[0].function() == core::detail::function_symbols::SetComp;
}

/// \brief Returns true if the term t is a bag comprehension
inline bool is_bag_comprehension(const atermpp::aterm& x)
{
  return is_abstraction(x) && x[0].function() == core::detail::function_symbols::BagComp;
}

/// \brief Returns true if the term t is a set/bag comprehension.
inline bool is_untyped_set_or_bag_comprehension(const atermpp::aterm& x)
{
  return is_abstraction(x) && x[0].function() == core::detail::function_symbols::UntypedSetBagComp;
}

/// \brief Returns true if the term t is a function symbol.
inline bool is_function_symbol(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::OpId;
}

/// \brief Returns true if the term t is a variable.
inline bool is_variable(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::DataVarId;
}

/// \brief Returns true if the term t is a machine_number.
inline bool is_machine_number(const atermpp::aterm& x)
{
  return x.type_is_int();
}

/// \brief Returns true if the term t is an application
/// \details This function is inefficient as the arity of a term must 
///          be determined and an inspection must take place in an 
///          array of function symbols. Therefore, there is an more efficient overload
///          is_application(const data_expression& x).
inline bool is_application(const atermpp::aterm& x)
{
  return core::detail::gsIsDataAppl(x);
}

/// \brief Returns true if the term t is an application, but it does not check
///        whether an application symbol of sufficient arity exists, assuming
///        this is ok.
inline bool is_application_no_check(const atermpp::aterm& x)
{
  return core::detail::gsIsDataAppl_no_check(x);
}

/// \brief Returns true if the term t is a where clause
inline bool is_where_clause(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Whr;
}

/// \brief Returns true if the term t is an identifier
inline bool is_untyped_identifier(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::UntypedIdentifier;
}

class application; // prototype

/// \brief data expression.
///
/// A data expression can be any of:
/// - variable
/// - function symbol
/// - application
/// - abstraction
/// - where clause
/// - set enumeration
/// - bag enumeration

//--- start generated class data_expression ---//
/// \\brief A data expression
class data_expression: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    data_expression()
      : atermpp::aterm(core::detail::default_values::DataExpr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit data_expression(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_DataExpr(*this));
    }

    /// Move semantics
    data_expression(const data_expression&) noexcept = default;
    data_expression(data_expression&&) noexcept = default;
    data_expression& operator=(const data_expression&) noexcept = default;
    data_expression& operator=(data_expression&&) noexcept = default;
//--- start user section data_expression ---//
    /// \brief A function to efficiently determine whether a data expression is 
    ///        made by the default constructor.
    bool is_default_data_expression() const
    {
      return *this==core::detail::default_values::DataExpr;
    }

    application operator()(const data_expression& e) const;

    application operator()(const data_expression& e1,
                           const data_expression& e2) const;

    application operator()(const data_expression& e1,
                           const data_expression& e2,
                           const data_expression& e3) const;

    application operator()(const data_expression& e1,
                           const data_expression& e2,
                           const data_expression& e3,
                           const data_expression& e4) const;

    application operator()(const data_expression& e1,
                           const data_expression& e2,
                           const data_expression& e3,
                           const data_expression& e4,
                           const data_expression& e5) const;

    application operator()(const data_expression& e1,
                           const data_expression& e2,
                           const data_expression& e3,
                           const data_expression& e4,
                           const data_expression& e5,
                           const data_expression& e6) const;

    /// \brief Returns the sort of the data expression
    sort_expression sort() const;

    private:
      // The following methods are not to be used in a data_expression.
      // They are restricted to a data_appl. 
      const_iterator begin() const;
      const_iterator end() const;

//--- end user section data_expression ---//
};

/// \\brief list of data_expressions
typedef atermpp::term_list<data_expression> data_expression_list;

/// \\brief vector of data_expressions
typedef std::vector<data_expression>    data_expression_vector;

// prototype declaration
std::string pp(const data_expression& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const data_expression& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(data_expression& t1, data_expression& t2)
{
  t1.swap(t2);
}
//--- end generated class data_expression ---//

inline void make_data_expression(data_expression& result)
{
  static_cast<atermpp::aterm&>(result)=core::detail::default_values::DataExpr;
}

/// \brief Test for a data_expression expression
/// \param x A term
/// \return True if it is a data_expression expression
inline
bool is_data_expression(const atermpp::aterm& x)
{
  return is_machine_number(x)                   ||
         is_lambda(x)                           ||
         is_forall(x)                           ||
         is_exists(x)                           ||
         is_set_comprehension(x)                ||
         is_bag_comprehension(x)                ||
         is_untyped_set_or_bag_comprehension(x) ||
         is_function_symbol(x)                  ||
         is_variable(x)                         ||
         is_application(x)                      ||
         is_where_clause(x)                     ||
         is_untyped_identifier(x);
}

/// \brief Converts an container with data expressions to data_expression_list
/// \param r A range of data expressions.
/// \note This function uses implementation details of the iterator type
/// and hence is sometimes efficient than copying all elements of the list.
template < typename Container >
inline data_expression_list make_data_expression_list(Container const& r, typename atermpp::enable_if_container< Container, data_expression >::type* = nullptr)
{
  return data_expression_list(r.begin(),r.end());
}

class variable;

// template function overloads
std::string pp(const data_expression_list& x, bool precedence_aware = true);
std::string pp(const data_expression_vector& x, bool precedence_aware = true);
data::data_expression translate_user_notation(const data::data_expression& x);
std::set<data::sort_expression> find_sort_expressions(const data::data_expression& x);
std::set<data::variable> find_all_variables(const data::data_expression& x);
std::set<data::variable> find_all_variables(const data::data_expression_list& x);
std::set<data::variable> find_free_variables(const data::data_expression& x);
std::set<data::variable> find_free_variables(const data::data_expression_list& x);
bool search_variable(const data::data_expression& x, const data::variable& v);

inline
bool is_constant(const data_expression& x)
{
  return find_free_variables(x).empty();
}

typedef atermpp::term_list<variable> variable_list;
variable_list free_variables(const data_expression& x);

} // namespace data

} // namespace mcrl2

// The trick of including application.h only at this point is needed to break
// the circular dependencies between application and data_expression. This
// dependency was introduced to allow the application operator to be defined for
// all data expression types.
#ifndef MCRL2_DATA_APPLICATION_H
#include "mcrl2/data/application.h"
#endif

namespace mcrl2
{
namespace data
{
/// \brief Returns true if the term t is an application.
/// \param t The variable that is checked. 
inline bool is_application(const data_expression& t)
{
  return !(is_machine_number(t) ||
           is_function_symbol(t) ||
           is_variable(t) ||
           is_where_clause(t) ||
           is_abstraction(t) ||
           is_untyped_identifier(t));
} 

/// \brief Transform a variable_list into a data_expression_list
inline 
const data_expression_list& variable_list_to_data_expression_list(const variable_list& l)
{
  return atermpp::down_cast<data_expression_list>(static_cast<const atermpp::aterm&>(l));
}

/// \brief Apply a data expression to a data expression.
inline
application data_expression::operator()(const data_expression& e) const
{
  return application(*this, e);
}

/// \brief Apply a data expression to two data expressions
inline
application data_expression::operator()(const data_expression& e1, const data_expression& e2) const
{
  return application(*this, e1, e2);
}

/// \brief Apply a data expression to three data expressions
inline
application data_expression::operator()(const data_expression& e1, const data_expression& e2, const data_expression& e3) const
{
  return application(*this, e1, e2, e3);
}

/// \brief Apply a data expression to four data expressions
inline
application data_expression::operator()(const data_expression& e1, const data_expression& e2, const data_expression& e3, const data_expression& e4) const
{
  return application(*this, e1, e2, e3, e4);
}

/// \brief Apply a data expression to five data expressions
inline
application data_expression::operator()(const data_expression& e1, const data_expression& e2, const data_expression& e3, const data_expression& e4, const data_expression& e5) const
{
  return application(*this, e1, e2, e3, e4, e5);
}

/// \brief Apply a data expression to six data expressions
inline
application data_expression::operator()(const data_expression& e1, const data_expression& e2, const data_expression& e3, 
                                        const data_expression& e4, const data_expression& e5, const data_expression& e6) const
{
  return application(*this, e1, e2, e3, e4, e5, e6);
}

} // namespace data
} // namespace mcrl2


namespace std
{

template<>
struct hash<mcrl2::data::data_expression>
{
    std::size_t operator()(const mcrl2::data::data_expression& v) const
    {
      const hash<atermpp::aterm> hasher;
      return hasher(v);
    }
};

} // namespace std



#endif // MCRL2_DATA_DATA_EXPRESSION_H

