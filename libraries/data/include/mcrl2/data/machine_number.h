// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/machine_number.h
/// \brief The class machine_number, which is a subclass of data_expression.

#ifndef MCRL2_DATA_MACHINE_NUMBER_H
#define MCRL2_DATA_MACHINE_NUMBER_H

#include "mcrl2/atermpp/aterm_list.h"
#include "mcrl2/data/data_expression.h"

namespace mcrl2
{

namespace data
{

typedef std::pair<atermpp::aterm, atermpp::aterm> machine_number_key_type;


//--- start generated class machine_number ---//
/// \\brief A machine number
class machine_number: public data_expression
{
  public:
    /// \\brief Default constructor X2.
    machine_number()
      : data_expression(atermpp::aterm_int(std::size_t(0)))
    {}

    /// \\brief Constructor based on an aterm.
    /// \param term A term
    explicit machine_number(const atermpp::aterm& term)
      : data_expression(term)
    {
      assert(this->type_is_int());
    }

    /// \\brief Constructor Z13.
    machine_number(std::size_t value)
      : data_expression(atermpp::aterm_int(value))
    {}

    /// Move semantics
    machine_number(const machine_number&) noexcept = default;
    machine_number(machine_number&&) noexcept = default;
    machine_number& operator=(const machine_number&) noexcept = default;
    machine_number& operator=(machine_number&&) noexcept = default;

    std::size_t value() const
    {
      return atermpp::down_cast<atermpp::aterm_int>(static_cast<const atermpp::aterm&>(*this)).value();
    }
};

/// \\brief Make_machine_number constructs a new term into a given address.
/// \\ \param t The reference into which the new machine_number is constructed. 
template <class... ARGUMENTS>
inline void make_machine_number(atermpp::aterm& t, size_t n)
{
  atermpp::make_aterm_int(reinterpret_cast<atermpp::aterm_int&>(t), n);
}

/// \\brief list of machine_numbers
typedef atermpp::term_list<machine_number> machine_number_list;

/// \\brief vector of machine_numbers
typedef std::vector<machine_number>    machine_number_vector;

// prototype declaration
std::string pp(const machine_number& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const machine_number& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(machine_number& t1, machine_number& t2)
{
  t1.swap(t2);
}
//--- end generated class machine_number ---//

/// \brief A string representation indicating the maximal machine number + 1. 
/// \return A string containing the number of available machine numbers. 
inline std::string max_machine_number_string()
{
  std::size_t n=std::numeric_limits<std::size_t>::max();
  assert(n % 100 != 99); // This is not possible, as n has the shape 2^n-1, but we check it anyhow. 
  return std::to_string(n/100) + std::to_string(1+n % 100);
}



// template function overloads
std::string pp(const machine_number_list& x, bool precedence_aware = true);
std::string pp(const machine_number_vector& x, bool precedence_aware = true);
std::set<data::variable> find_all_variables(const data::machine_number& x);

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_MACHINE_NUMBER_H

