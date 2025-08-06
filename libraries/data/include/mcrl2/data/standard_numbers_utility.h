// Author(s): Jeroen van der Wulp
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/standard_numbers_utility.h
/// \brief Provides utilities for working with data expressions of standard sorts

#ifndef MCRL2_DATA_STANDARD_NUMBERS_UTILITY_H
#define MCRL2_DATA_STANDARD_NUMBERS_UTILITY_H

#include <cmath>

#include "mcrl2/data/standard_utility.h"

namespace mcrl2::data
{

/// \cond INTERNAL_DOCS
namespace detail
{
// Convert to number represented as character array where each character
// represents a decimal digit
inline std::vector< char > string_to_vector_number(const std::string& s)
{
  assert(s.size() > 0);
  std::vector< char > result;

  result.reserve(s.size());

  for (char i: s)
  {
    assert('0' <= i && i <= '9');

    result.push_back(i - '0');
  }

  return result;
}

// Convert from number represented as character array where each character
// represents a decimal digit
inline std::string vector_number_to_string(const std::vector< char >& v)
{
  assert(v.size() > 0);
  std::string result;

  result.reserve(v.size());

  for (char i: v)
  {
    result.push_back(i + '0');
  }

  return result;
}

/// \brief Multiplies a number in decimal notation represented by an array by two
/// \param[in,out] number the number
/// A number d0 d1 ... dn is represented as s[0] s[1] ... s[n]
inline void decimal_number_multiply_by_two(std::vector< char >& number)
{
  assert(0 < number.size());

  std::vector< char >           result(number.size() + 2, 0);
  std::vector< char >::iterator j(result.begin());

  if (5 <= number[0])
  {
    *(j++) = number[0] / 5;
  }

  for (std::vector< char >::const_iterator i = number.begin(); i < number.end(); ++i, ++j)
  {
    // result[a] = 2*(number[b] mod 5) + number[b+1] div 5   where result[a] = *j and number[b] = *(i)
    if (i == number.end() - 1)
    {
      *j = 2 * (*i % 5);
    }
    else
    {
      *j = 2 * (*i % 5) + *(i+1) / 5;
    }
  }

  result.resize(j - result.begin());

  number.swap(result);
}


/// \brief Adds one to a number in decimal notation represented by an array
/// \param[in,out] number the number
/// A number d0 d1 ... dn is represented as s[0] s[1] ... s[n]
inline void decimal_number_increment(std::vector< char >& number)
{
  assert(0 < number.size());

  for (std::vector< char >::reverse_iterator i = number.rbegin(); i != number.rend(); ++i)
  {
    if (*i < 9)
    {
      ++(*i);

      return;
    }
    else
    {
      *i = 0;
    }
  }

  number.insert(number.begin(), 1);
}


/// Type T is an unsigned type
template< typename T >
inline std::string as_decimal_string(T t)
{
  if (t != 0)
  {
    std::string result;

    while (0 < t)
    {
      result.append(1, '0' + static_cast< char >(t % 10));

      t /= 10;
    }

    return std::string(result.rbegin(), result.rend());
  }

  return "0";
}

/// \brief Divides a number in decimal notation represented by an array by two
/// \param[in,out] s the number
/// A number d0 d1 ... dn is represented as s[0] s[1] ... s[n]
inline void decimal_number_divide_by_two(std::vector< char >& number)
{
  assert(0 < number.size());

  std::vector< char >           result(number.size(), 0);
  std::vector< char >::iterator j(result.begin());

  if (2 <= number[0])
  {
    *(j++) = number[0] / 2;
  }

  for (std::vector< char >::const_iterator i = number.begin() + 1; i != number.end(); ++i, ++j)
  {
    // result[a] = 5*(number[b - 1] mod 2) + number[b] div 2   where result[a] = *j, number[b - 1] = *(i - 1)
    *j = 5 * (*(i - 1) % 2) + *i / 2;
  }

  result.resize(j - result.begin());

  number.swap(result);
}

#ifdef MCRL2_ENABLE_MACHINENUMBERS
// Check whether a number vector is equal to zero.
inline bool is_zero_number_vector(const std::vector< std::size_t >& v)
{
  for(const std::size_t n: v)
  {
    if (n>0)
    {
      return false;
    }
  }
  // All digits are zero. Thus the number is zero.
  assert(v.size()<2);  // Zero can only occur as a vector of length 0 or 1. 
  return true;
}


// Add the two number vectors. The least significant digit is at position 0.
inline std::vector<std::size_t> add_number_vectors(const std::vector< std::size_t >& v1,
                                                   const std::vector< std::size_t >& v2)
{
  std::vector<std::size_t> result;
  std::size_t carry=0;
  for(std::size_t i=0; i<std::max(v1.size(),v2.size()); i++)
  {
    std::size_t n1 = (i<v1.size()? v1[i]: 0);
    std::size_t n2 = (i<v2.size()? v2[i]: 0);

    std::size_t sum = n1+n2+carry;
    if (sum<n1)  // There is a carry overflow. 
    {
      assert(sum<n2);
      carry=1;
    }
    else
    {
      carry=0;
    } 
    result.push_back(sum);
  } 
  if (carry==1)
  {    
    result.push_back(1);
  } 
  assert(result.size()<=1 || result.back()!=0);
  return result;
}   


inline std::vector<size_t> multiply_by10_and_add(const std::vector< size_t >& v, size_t val)
{ 
  std::vector<size_t> val_vector(1,val);
  assert(val_vector.size()==1);
  std::vector<size_t> result=add_number_vectors(v,v);   // result is 2 v.
  result = add_number_vectors(result, result);          // result is 4 v.
  result = add_number_vectors(result, v);               // result is 5 v.
  result = add_number_vectors(result, result);          // result is 10 v.
  result = add_number_vectors(result, val_vector);      // result is 10 v + val.
  return result;
}

// Transform a number vector to a string representation of the same number.
inline std::string number_vector_as_string(const std::vector< std::size_t >& v)
{ 
  assert(!detail::is_zero_number_vector(v));
  
  std::vector< char > result = detail::string_to_vector_number("0");
  bool is_zero=true; // Avoid doing many multiplications for small numbers. 

  for(size_t i=v.size(); i>0 ; --i)
  {
    std::size_t n= v.at(i-1);
    for(std::size_t mask = std::size_t(1)<<(8*sizeof(std::size_t)-1);  mask>0; mask=mask>>1)
    {     
      if (!is_zero) 
      {  
        detail::decimal_number_multiply_by_two(result);
      }
      if ((n & mask) >0)
      {
        detail::decimal_number_increment(result);
        is_zero=false;
      }
    }
  }

  return detail::vector_number_to_string(result);
}


// Multiply a vector consisting of size_t's with 10 and add the extra value.
// Convert to number represented as character array where each character
// represents a decimal digit
inline std::vector< size_t > number_string_to_vector_number(const std::string& s)
{
  assert(s.size() > 0);
  std::vector< size_t > result;

  result.reserve(s.size()/18);  // approximately 18 digits fit in one size_t.
  
  for (char i: s)
  {
    if ('0' <= i && i <= '9')
    {
       result = multiply_by10_and_add(result, i - '0');
    }
    else
    {
      throw mcrl2::runtime_error("The string " + s + " is expected to only consist of digits. ");
    }
  }
  
  return result;
} 
#endif

} // namespace detail
/// \endcond

namespace sort_pos
{
/// \brief Constructs expression of type Bool from an integral type
/// Type T is an unsigned integral type.
template <typename T>
inline data_expression pos(const T t)
  requires std::is_integral_v<T>
{
  assert(t>0);

#ifdef MCRL2_ENABLE_MACHINENUMBERS
  static_assert(sizeof(T)<=sizeof(std::size_t),"Can only convert numbers up till a size_t.");
  return sort_pos::most_significant_digit(machine_number(t));
#else
  std::vector< bool > bits;
  bits.reserve(8 * sizeof(T));

  for (T u = t; 1 < u; u /= 2)
  {
    bits.push_back(u % 2 != 0);
  }

  data_expression result(sort_pos::c1());

  for (std::vector< bool >::reverse_iterator i = bits.rbegin(); i != bits.rend(); ++i)
  {
    result = sort_pos::cdub(sort_bool::bool_(*i), result);
  }

  return result;
#endif
}

/// \brief Constructs expression of type Pos from a string
/// \param n A string
inline data_expression pos(const std::string& n)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  std::vector<std::size_t> number_vector= detail::number_string_to_vector_number(n);

  assert(!detail::is_zero_number_vector(number_vector));

  data_expression result=most_significant_digit(machine_number(number_vector.back()));
  for(std::size_t i=number_vector.size()-1; i>0; --i)
  {
    result = sort_pos::concat_digit(result,machine_number(number_vector.at(i-1)));
  }
  return result;

#else
  std::vector< char > number_as_vector(detail::string_to_vector_number(n));

  std::vector< bool > bits;

  while (0 < number_as_vector.size() && !((number_as_vector.size() == 1) && number_as_vector[0] == 1))   // number != 1
  {
    bits.push_back((static_cast< int >(*number_as_vector.rbegin()) % 2 != 0));

    detail::decimal_number_divide_by_two(number_as_vector);
  }
  data_expression result(sort_pos::c1());

  for (std::vector< bool >::reverse_iterator i = bits.rbegin(); i != bits.rend(); ++i)
  {
    result = sort_pos::cdub(sort_bool::bool_(*i), result);
  }

  return result;
#endif
}

/// \brief Determines whether n is a positive constant
/// \param n A data expression
inline bool is_positive_constant(const data_expression& n)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  return (sort_pos::is_most_significant_digit_application(n) &&
          is_machine_number(sort_pos::arg(n)))   ||
         (sort_pos::is_concat_digit_application(n) &&
          sort_pos::is_positive_constant(sort_pos::arg1(n))
         );
#else
  return sort_pos::is_c1_function_symbol(n) ||
         (sort_pos::is_cdub_application(n) &&
          sort_bool::is_boolean_constant(sort_pos::left(n)) &&
          sort_pos::is_positive_constant(sort_pos::right(n))
         );
#endif
}

/// \brief Return the string representation of a positive number
/// \param n_in A data expression
/// \pre is_positive_constant(n)
/// \return String representation of n
/// Transforms a positive constant n into a character array containing
/// the decimal representation of n.
inline
std::string positive_constant_as_string(const data_expression& n_in)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  std::vector<std::size_t> number_vector;
  data_expression n=n_in;

  while (sort_pos::is_concat_digit_application(n))
  {
    number_vector.push_back(atermpp::down_cast<machine_number>(sort_pos::arg2(n)).value());
    n = sort_pos::arg1(n);
  }

  assert(sort_pos::is_most_significant_digit_application(n));
  number_vector.push_back(atermpp::down_cast<machine_number>(sort_pos::arg(n)).value());

  return detail::number_vector_as_string(number_vector);
#else
  std::vector<bool> bits;
  data_expression n=n_in;

  while (sort_pos::is_cdub_application(n))
  {
    bits.push_back(sort_bool::is_true_function_symbol(sort_pos::left(n)));
    n = sort_pos::right(n);
  }

  assert(sort_pos::is_c1_function_symbol(n));

  std::vector< char > result = detail::string_to_vector_number("1");

  for (std::vector<bool>::reverse_iterator i = bits.rbegin(); i != bits.rend(); ++i)
  {
    detail::decimal_number_multiply_by_two(result);
    if (*i)
    {
      detail::decimal_number_increment(result);
    }
  }

  return detail::vector_number_to_string(result);
#endif
}

/// \brief Returns the NUMERIC_TYPE representation of a positive number
/// \param n_in A data expression
/// \pre is_positive_constant(n)
/// \return Representation of n as sort NUMERIC_TYPE
/// Transforms a positive constant n into number of sort NUMERIC_TYPE. Throws an exception when it does not fit. 
/// the decimal representation of n.
template <class NUMERIC_TYPE>
inline
NUMERIC_TYPE positive_constant_to_value(const data_expression& n)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  if constexpr (std::is_integral_v<NUMERIC_TYPE>)
  {
    if (is_concat_digit_application(n))
    {
      throw mcrl2::runtime_error("Number " + pp(n) + " is too large to transform to a machine number.");
    }
    assert(is_most_significant_digit_application(n));
    assert(atermpp::down_cast<machine_number>(sort_pos::arg(n)).value()>0);
    return atermpp::down_cast<machine_number>(sort_pos::arg(n)).value();
  }
  else
  {
    if (is_concat_digit_application(n))
    {
      return positive_constant_to_value<NUMERIC_TYPE>(sort_pos::arg1(n))*pow(2.0,64) +
             atermpp::down_cast<machine_number>(sort_pos::arg2(n)).value();
    }
    assert(is_most_significant_digit_application(n));
    return atermpp::down_cast<machine_number>(sort_pos::arg(n)).value();
  }
#else
  if (sort_pos::is_cdub_application(n))
  {
    NUMERIC_TYPE result = 2*positive_constant_to_value<NUMERIC_TYPE>(sort_pos::right(n));
    if (sort_bool::is_true_function_symbol(sort_pos::left(n)))
    {
      result=result+1;
    }
    return result;
  }

  assert(sort_pos::is_c1_function_symbol(n));

  return static_cast<NUMERIC_TYPE>(1);
#endif
}

} // namespace sort_pos

namespace sort_nat
{

/// \brief Constructs expression of type pos from an integral type
template <typename T>
inline data_expression nat(T t)
  requires std::is_integral_v<T>
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  static_assert(sizeof(T)<=sizeof(std::size_t),"Can only convert numbers up till a size_t.");
  return sort_nat::most_significant_digit_nat(machine_number(t));
#else
  if (t == 0)
  {
    return sort_nat::c0();
  }
  return sort_nat::cnat(sort_pos::pos(t));
#endif
}

/// \brief Constructs expression of type Nat from a string
/// \param n A string
inline data_expression nat(const std::string& n)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  std::vector<std::size_t> number_vector= detail::number_string_to_vector_number(n);

  if (number_vector.empty())
  {
    return most_significant_digit_nat(machine_number(0));
  }
  data_expression result=most_significant_digit_nat(machine_number(number_vector.back()));
  for(std::size_t i=number_vector.size()-1; i>0; --i)
  {
    result = sort_nat::concat_digit(result,machine_number(number_vector.at(i-1)));
  }
  return result;
#else
  std::vector< char > number_as_vector(detail::string_to_vector_number(n));
    
  std::vector< bool > bits;

  while (0 < number_as_vector.size() && !((number_as_vector.size() == 1) && number_as_vector[0] == 1))   // number != 1
  {
    bits.push_back((static_cast< int >(*number_as_vector.rbegin()) % 2 != 0));

    detail::decimal_number_divide_by_two(number_as_vector);
  }
  return (n == "0") ? sort_nat::c0() : static_cast< data_expression const& >(sort_nat::cnat(sort_pos::pos(n)));
#endif
}

/// \brief Determines whether n is a natural constant
/// \param n A data expression
inline bool is_natural_constant(const data_expression& n)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  return (sort_nat::is_most_significant_digit_nat_application(n) &&
          is_machine_number(sort_nat::arg(n))) 
         ||
         (sort_nat::is_concat_digit_application(n) &&
          sort_nat::is_natural_constant(sort_nat::arg1(n))
         );
#else
  return sort_nat::is_c0_function_symbol(n) ||
         (sort_nat::is_cnat_application(n) &&
          sort_pos::is_positive_constant(sort_nat::arg(n))
         );
#endif
}

/// \brief Return the string representation of a natural number
/// \param n A data expression
/// \pre is_natural_constant(n)
/// \return String representation of n
inline std::string natural_constant_as_string(const data_expression& n_in)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  std::vector<std::size_t> number_vector;
  data_expression n=n_in;
  
  while (sort_nat::is_concat_digit_application(n))
  { 
    number_vector.push_back(atermpp::down_cast<machine_number>(sort_nat::arg2(n)).value());
    n = sort_nat::arg1(n);
  }

  assert(sort_nat::is_most_significant_digit_nat_application(n));
  number_vector.push_back(atermpp::down_cast<machine_number>(sort_nat::arg(n)).value());
    
  if (detail::is_zero_number_vector(number_vector))
  { 
    return "0";
  }
  
  return detail::number_vector_as_string(number_vector);
#else
  assert(is_natural_constant(n_in));
  if (sort_nat::is_c0_function_symbol(n_in))
  {
    return "0";
  }
  else
  {
    return sort_pos::positive_constant_as_string(sort_nat::arg(n_in));
  }
#endif
}

/// \brief Return the NUMERIC_VALUE representation of a natural number
/// \param n A data expression
/// \pre is_natural_constant(n)
/// \return NUMERIC_VALUE representation of n
template <class NUMERIC_TYPE>
inline NUMERIC_TYPE natural_constant_to_value(const data_expression& n)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  if constexpr (std::is_integral_v<NUMERIC_TYPE>)
  {
    if (is_concat_digit_application(n))
    {
      throw mcrl2::runtime_error("Number " + pp(n) + " is too large to transform to a machine number.");
    }
    assert(is_most_significant_digit_nat_application(n));
    return atermpp::down_cast<machine_number>(sort_nat::arg(n)).value();
  }
  else
  {
    if (is_concat_digit_application(n))
    {
      return natural_constant_to_value<NUMERIC_TYPE>(sort_nat::arg1(n))*pow(2.0,64) +
             atermpp::down_cast<machine_number>(sort_nat::arg2(n)).value();
    }
    assert(is_most_significant_digit_nat_application(n));
    return atermpp::down_cast<machine_number>(sort_nat::arg(n)).value();
  }
#else
  assert(is_natural_constant(n));
  if (sort_nat::is_c0_function_symbol(n))
  {
    return static_cast<NUMERIC_TYPE>(0);
  }
  else
  {
    return sort_pos::positive_constant_to_value<NUMERIC_TYPE>(sort_nat::arg(n));
  }
#endif
}

#ifdef MCRL2_ENABLE_MACHINENUMBERS
/// \brief Translate a positive constant to a natural constant;
/// \returns A natural constant with the same value as the positive constant;
inline data_expression transform_positive_constant_to_nat(const data_expression& n)
{ 
  assert(n.sort()==sort_pos::pos());
  if (sort_pos::is_concat_digit_application(n))
  {
    return sort_nat::concat_digit(transform_positive_constant_to_nat(sort_pos::arg1(n)),sort_nat::arg2(n));
  }

  assert(sort_pos::is_most_significant_digit_application(n));
  return sort_nat::most_significant_digit_nat(sort_pos::arg(n));
}

/// \brief Translate a positive constant to a natural constant;
/// \returns A natural constant with the same value as the positive constant;
inline data_expression transform_positive_number_to_nat(const data_expression& n)
{ 
  if (sort_pos::is_positive_constant(n))
  {
    return transform_positive_constant_to_nat(n);
  }
  if (is_function_symbol(n))
  {
    const function_symbol& f=atermpp::down_cast<function_symbol>(n);
    const std::string& name(f.name());
    if (name.find_first_not_of("-/0123456789") == std::string::npos) // crude but efficient
    {
      return sort_nat::nat(name);
    }
  }
  return application(sort_nat::pos2nat(),n);
}
#endif

} // namespace sort_nat

namespace sort_int
{

/// \brief Constructs expression of type pos from an integral type
template <typename T>
inline data_expression int_(T t)
  requires(std::is_integral_v<T> && std::is_unsigned_v<T>) { return sort_int::cint(sort_nat::nat(t)); }

/// \brief Constructs expression of type pos from an integral type.
template <typename T>
inline data_expression int_(T t)
  requires(std::is_integral_v<T> && std::is_signed_v<T>) {
    if (t < 0)
    {
      return sort_int::cneg(sort_pos::pos(std::make_unsigned_t<T>(-t)));
    }
    return sort_int::cint(sort_nat::nat(std::make_unsigned_t<T>(t)));
  }

/// \brief Constructs expression of type Int from a string.
/// \param n A string.
/// \pre n is of the form ([-]?[0...9][0...9]+)([0...9]+).
inline data_expression int_(const std::string& n)
{
  if (n[0] == '-')
  {
    return sort_int::cneg(sort_pos::pos(n.substr(1)));
  }
  return sort_int::cint(sort_nat::nat(n));
}

/// \brief Determines whether n is an integer constant.
/// \param n A data expression.
inline bool is_integer_constant(const data_expression& n)
{
  return (sort_int::is_cint_application(n) &&
          sort_nat::is_natural_constant(sort_int::arg(n))) ||
         (sort_int::is_cneg_application(n) &&
          sort_pos::is_positive_constant(sort_int::arg(n))
         );
}

/// \brief Return the string representation of an integer number.
/// \param n A data expression.
/// \pre is_integer_constant(n).
/// \return String representation of n.
inline std::string integer_constant_as_string(const data_expression& n)
{
  assert(is_integer_constant(n));
  if (sort_int::is_cint_application(n))
  {
    return sort_nat::natural_constant_as_string(sort_int::arg(n));
  }
  else
  {
    return "-" + sort_pos::positive_constant_as_string(sort_int::arg(n));
  }
}

/// \brief Return the NUMERIC_VALUE representation of an integer number.
/// \param n A data expression.
/// \pre is_integer_constant(n).
/// \return NUMERIC_VALUE representation of n.
template <class NUMERIC_VALUE>
inline NUMERIC_VALUE integer_constant_to_value(const data_expression& n)
{
  assert(is_integer_constant(n));
  if (sort_int::is_cint_application(n))
  {
    return sort_nat::natural_constant_to_value<NUMERIC_VALUE>(sort_int::arg(n));
  }
  else
  {
    return -sort_pos::positive_constant_to_value<NUMERIC_VALUE>(sort_int::arg(n));
  }
}
} // namespace sort_int

namespace sort_real
{
/// \brief Constructs expression of type Real from an integral type
/// \param t An expression of type T.
template <typename T>
inline data_expression real_(T t)
  requires std::is_integral_v<T>
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  return sort_real::creal(sort_int::int_(t), sort_pos::pos(1));
#else
  return sort_real::creal(sort_int::int_(t), sort_pos::c1());
#endif
}

/// \brief Constructs expression of type Real from an integral type
/// \param numerator numerator.
/// \param denominator denominator.
template <typename T>
inline data_expression real_(T numerator, T denominator)
  requires std::is_integral_v<T>
{
  return sort_real::creal(sort_int::int_(numerator), sort_pos::pos(denominator));
}

/// \brief Constructs expression of type Real from two number strings.
/// \param numerator numerator.
/// \param denominator denominator.
/// \return numerator / denominator.
inline data_expression real_(const std::string& numerator, const std::string& denominator)
{
  return sort_real::creal(sort_int::int_(numerator), sort_pos::pos(denominator));
}

/// \brief Constructs expression of type Real from a string
/// \param n A string
/// \pre n is of the form (-[1...9][0...9]+)([0...9]+)
inline data_expression real_(const std::string& n)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  return sort_real::creal(sort_int::int_(n), sort_pos::pos(1));
#else
  return sort_real::creal(sort_int::int_(n), sort_pos::c1());
#endif
}

/// \brief Yields the real value of a data expression.
/// \param r A data expression of sort real in normal form.
template <class NUMERIC_TYPE>
inline NUMERIC_TYPE value(const data_expression& r,
    std::enable_if_t<std::is_floating_point_v<NUMERIC_TYPE>>* = nullptr)
{
  if (is_creal_application(r))
  {
    const application& a = atermpp::down_cast<application>(r);
    return sort_int::integer_constant_to_value<NUMERIC_TYPE>(a[0]) / 
           sort_pos::positive_constant_to_value<NUMERIC_TYPE>(a[1]);
  }
  throw runtime_error("Expected a closed term of type real " + pp(r) + ".");
}

} // namespace sort_real

/// \brief Construct numeric expression from a string representing a number in decimal notation
/// \param s A sort expression
/// \param n A string
/// \pre n is of the form [1]?[0...9]+
inline data_expression number(const sort_expression& s, const std::string& n)
{
  if (s == sort_pos::pos())
  {
    return sort_pos::pos(n);
  }
  else if (s == sort_nat::nat())
  {
    return sort_nat::nat(n);
  }
  else if (s == sort_int::int_())
  {
    return sort_int::int_(n);
  }

  return sort_real::real_(n);
}

/// \brief Returns true if and only if s1 == s2, or if s1 is a less specific numeric type than s2
///
/// \param[in] s1 a sort expression
/// \param[in] s2 a sort expression
inline bool is_convertible(const sort_expression& s1, const sort_expression& s2)
{
  if (s1 != s2)
  {
    if (s2 == sort_real::real_())
    {
      return s1 == sort_int::int_() || s1 == sort_nat::nat() || s1 == sort_pos::pos();
    }
    else if (s2 == sort_int::int_())
    {
      return s1 == sort_nat::nat() || s1 == sort_pos::pos();
    }
    else if (s2 == sort_nat::nat())
    {
      return s1 == sort_pos::pos();
    }

    return false;
  }

  return true;
}

} // namespace mcrl2::data

#endif

