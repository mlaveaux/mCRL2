// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

/** \file
 *
 * \brief This file contains a class big_natural_number that stores big positive numbers of arbitrary size.
 *        It has all common operations that one can expect on big numbers.
 * \author Jan Friso Groote
 */


#ifndef MCRL2_UTILITIES_BIG_NUMBERS_H
#define MCRL2_UTILITIES_BIG_NUMBERS_H

#include "mcrl2/utilities/hash_utility.h"

#include <algorithm>
#include <limits>
#include <ranges>
#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

// Prototype.

namespace mcrl2::utilities
{
class big_natural_number;

inline std::string pp(const big_natural_number& l);

} // namespace mcrl2::utilities

namespace mcrl2::utilities
{
namespace detail
{

  // Calculate <carry,result>:=n1+n2+carry. The carry can be either 0 or 1, both
  // at the input and the output.
  inline std::size_t add_single_number(const std::size_t n1, const std::size_t n2, std::size_t& carry)
  {
    assert(carry<=1);
    std::size_t result=n1+n2+carry;
    if (carry==0)
    {
      if (result<n1)
      {
        carry=1;
      }
    }
    else // carry==1
    {
      if (result>n1)
      {
        carry=0;
      }
    }
    return result;
  }

  // Calculate <carry,result>:=n1-n2-carry. The carry can be either 0 or 1, both
  // at the input and the output. If the carry is 1, this indicated that 1 must be subtracted.
  inline std::size_t subtract_single_number(const std::size_t n1, const std::size_t n2, std::size_t& carry)
  {
    assert(carry<=1);
    std::size_t result=n1-n2-carry;
    if (carry==0)
    {
      if (result>n1)
      {
        carry=1;
      }
    }
    else // carry==1
    {
      if (result<n1)
      {
        carry=0;
      }
    }
    return result;
  }
  
  // Calculate <carry,result>:=n1*n2+carry, where the lower bits of the calculation
  // are stored in the result, and the higher bits are stored in carry.
  inline std::size_t multiply_single_number(const std::size_t n1, const std::size_t n2, std::size_t& multiplication_carry)
  {
    // TODO: It is more concise and efficient to use 128bit machine calculation when available.
    const int no_of_bits_per_digit=std::numeric_limits<std::size_t>::digits;

    // split input numbers into no_of_bits_per_digit/2 digits
    std::size_t n1ls = n1 & ((1LL<<(no_of_bits_per_digit/2))-1);
    std::size_t n1ms = n1 >> (no_of_bits_per_digit/2);
    std::size_t n2ls = n2 & ((1LL<<(no_of_bits_per_digit/2))-1);
    std::size_t n2ms = n2 >> (no_of_bits_per_digit/2);
    
    // First calculate the result of the least significant no_of_bits_per_digit.
    std::size_t local_carry=0;
    std::size_t result = add_single_number(n1ls*n2ls,multiplication_carry,local_carry);
    std::size_t cumulative_carry=local_carry;
    local_carry=0;
    result=add_single_number(result,((n1ms*n2ls)<<(no_of_bits_per_digit/2)),local_carry);
    cumulative_carry=cumulative_carry+local_carry;
    local_carry=0;
    result=add_single_number(result,((n1ls*n2ms)<<(no_of_bits_per_digit/2)),local_carry);
    cumulative_carry=cumulative_carry+local_carry;

    // Now calculate the result of the most significant no_of_bits_per_digit.
    multiplication_carry=cumulative_carry;
    multiplication_carry=multiplication_carry+((n1ms*n2ls)>>(no_of_bits_per_digit/2));
    multiplication_carry=multiplication_carry+((n1ls*n2ms)>>(no_of_bits_per_digit/2));
    multiplication_carry=multiplication_carry+n1ms*n2ms;

    return result;
  }
  
  // Calculate <result,remainder>:=(remainder * 2^64 + p) / q assuming the result
  // fits in 64 bits. More concretely, q>remainder. 
  inline std::size_t divide_single_number(const std::size_t p, const std::size_t q, std::size_t& remainder)
  {
    assert(q>remainder);
    const int no_of_bits_per_digit=std::numeric_limits<std::size_t>::digits;

    // Split input numbers into no_of_bits_per_digit/2 digits.
    // First get the least significant part.
    std::size_t pms = (p >> (no_of_bits_per_digit/2)) + (remainder << (no_of_bits_per_digit/2));
    
    // First divide the most significant part by q.
    std::size_t resultms = pms/q;
    assert((resultms >> (no_of_bits_per_digit/2)) == 0);
    std::size_t remainderms = pms%q;

    // Now obtain the least significant part.
    std::size_t pls = (p & ((1LL<<(no_of_bits_per_digit/2))-1)) + (remainderms << (no_of_bits_per_digit/2));
    
    // Second divide the least significant part by q.
    std::size_t resultls = pls/q;
    remainder = pls%q;
    assert((resultls >> (no_of_bits_per_digit/2)) == 0);

    return resultls + (resultms << (no_of_bits_per_digit/2));
  }

  // Calculate the common greates divisor of p and q. 
  inline std::size_t greatest_common_divisor(std::size_t p, std::size_t q)
  {
    if (p==0) { return q; }
    if (q==0) { return p; }
    if (p<q)
    {
      std::swap(p,q);
    }
    // p is the largest.
    std::size_t remainder=p % q; 
    while (remainder>0)
    {
      p=q; 
      q=remainder;
      remainder=p % q;
    }
    return q; 
  }

  inline void remove_common_divisor(std::size_t& p, std::size_t& q)
  {
    const std::size_t gcd=greatest_common_divisor(p,q);
    p = p/gcd;
    q = q/gcd;
    assert(greatest_common_divisor(p,q)==1); 
  }

} // namespace detail

// class big_natural_number;
inline std::string pp(const big_natural_number& l);

class big_natural_number
{
    friend std::hash<big_natural_number>;
    friend inline void swap(big_natural_number& x, big_natural_number& y) noexcept;

  protected:
    // Numbers are stored as std::size_t words, with the most significant number last. 
    // Note that the number representation is not unique. Numbers have no trailing
    // zero's, i.e., this->back()!=0 (if this->size()>0). Therefore their representation is unique.
    std::vector<std::size_t> m_number;

    /* Multiply the current number by n and add the carry */
    void multiply_by(std::size_t n, std::size_t carry)
    {
      for(std::size_t& i: m_number)
      {
        i=detail::multiply_single_number(i,n,carry);
      }
      if (carry)
      {
        /* Add an extra digit with the carry */
        m_number.push_back(carry);
      }
      remove_significant_digits_that_are_zero();
      is_well_defined();
    }

    // Remove zero's at the end of m_number vector that are
    // leading and irrelevant zero's in the number representation.
    void remove_significant_digits_that_are_zero()
    {
      for( ; m_number.size()>0 && m_number.back()==0 ; )
      {
        m_number.pop_back();
      }
    }

    // Check that the number is well defined, in the sense
    // that there are no most significant digits that are zero.
    void is_well_defined() const
    {
      assert(m_number.size()==0 || m_number.back()!=0);
    }

    // \brief This functions prints a number in internal represenation. This function is useful and only meant for debugging.
    void print_number(const std::string& s) const
    {
      std::cerr << s << "  " << m_number.size() << "\n";
      for(std::size_t i: m_number)
      {
        std::cerr << i << " ";
      }
      std::cerr << "\n---------------------\n";
    }

  public:
    // Default constructor. The value is 0 by default.
    explicit big_natural_number() = default;

    // Constructor based on a std::size_t. The value of the number will be n.
    explicit big_natural_number(const std::size_t n)
    {
      if (n>0)
      {
        m_number.push_back(n);
      }
      is_well_defined();
    }

    // Constructor based on a string. The string is a digital number.
    big_natural_number(const std::string& s)
    {
      for(char c: s)
      {
        if (!isdigit(c))
        {
          throw mcrl2::runtime_error("Number " + s + " contains symbol '" + c + "' which is not a digit.");
        }
        if (c>='0')
        {
          multiply_by(10,std::size_t(c-'0'));
        }
      }
      is_well_defined();
    }

    // Add a digit the this natural number. 
    void push_back(const std::size_t n)
    {
      m_number.push_back(n);
    }

    /** \brief Returns whether this number equals zero.
        \details This is more efficient than checking x==big_natural_number(0).
    */
    bool is_zero() const
    {
      is_well_defined();
      return m_number.size()==0;
    }

    /** \brief Returns whether this number equals a number of std::size_t.
        \details This is more efficient than checking x==big_natural_number(1).
    */
    bool is_number(std::size_t n) const
    {
      is_well_defined();
      if (n==0) 
      {
        return m_number.size()==0;
      }
      return m_number.size()==1 && m_number.front()==n;
    }

    /** \brief Sets the number to zero.
        \details This is more efficient than using an assignment x=0.
    */
    void clear()
    {
      is_well_defined();
      m_number.clear();
    }

    /** \brief Transforms this number to a std::size_t, provided it is sufficiently small.
               If not an mcrl2::runtime_error is thrown. 
    */
    explicit operator std::size_t() const
    {
      is_well_defined();
      if (m_number.size()>1)
      {
        throw mcrl2::runtime_error("It is not possible to transform a big natural number into a machine size number if it does not fit.");
      }
      if (m_number.size()==0)
      {
        return 0;
      }
      return m_number.front();
    }

    /** \brief Give the number of digits in this big number.
    */
    std::size_t size() const
    {
      return m_number.size();
    }

    /** \brief Give the n-th digit where the most significant digit is positioned last. 
    */
    std::size_t operator[](const std::size_t n) const
    {
      assert(n<m_number.size());
      return m_number.at(n);
    }

    /* \brief Standard comparison operator.
    */
    bool operator==(const big_natural_number& other) const
    {
      is_well_defined();
      other.is_well_defined(); 
      return m_number==other.m_number;
    }

    /* \brief Standard comparison operator.
    */
    bool operator!=(const big_natural_number& other) const
    {
      return !this->operator==(other);
    } 

    /* \brief Standard comparison operator.
    */
    bool operator<(const big_natural_number& other) const
    {
      is_well_defined();
      other.is_well_defined(); 
      if (m_number.size()<other.m_number.size())
      {
        return true;
      }
      if (m_number.size()>other.m_number.size())
      {
        return false;
      }
      assert(m_number.size()==other.m_number.size());
      std::vector<std::size_t>::const_reverse_iterator j=other.m_number.rbegin();
      for(std::vector<std::size_t>::const_reverse_iterator i=m_number.rbegin(); i!=m_number.rend(); ++i, ++j)
      {
        if (*i < *j)
        {
          return true;
        }
        if (*i > *j)
        {
          return false;
        }
      }
      // The numbers are equal.
      assert(m_number==other.m_number);
      return false;
    }

    /* \brief Standard comparison operator.
    */
    bool operator<=(const big_natural_number& other) const
    {
      return !this->operator>(other);
    }

    /* \brief Standard comparison operator.
    */
    bool operator>(const big_natural_number& other) const
    {
      return other.operator<(*this);
    }

    /* \brief Standard comparison operator.
    */
    bool operator>=(const big_natural_number& other) const
    {
      return other.operator<=(*this);
    }

    /* Divide the current number by n. If there is a remainder return it. */
    std::size_t divide_by(std::size_t n)
    {
      std::size_t remainder=0;
      for(std::vector<std::size_t>::reverse_iterator i=m_number.rbegin(); i!=m_number.rend(); ++i)
      {
        *i=detail::divide_single_number(*i,n,remainder);
      }
      remove_significant_digits_that_are_zero();
      is_well_defined();
      return remainder;
    }

    // Add the argument to this big natural number 
    void add(const big_natural_number& other)
    {
      is_well_defined();
      other.is_well_defined();

      // big_natural_number result;
      // result.m_number.reserve((std::max)(m_number.size(),other.m_number.size()));
      std::size_t carry=0;
      for(std::size_t i=0; i < (std::max)(m_number.size(),other.m_number.size()); ++i)
      {
        if (i>=m_number.size())
        {
          m_number.push_back(detail::add_single_number(0,other.m_number[i],carry));
        }
        else if (i>=other.m_number.size())
        {
          m_number[i]=detail::add_single_number(m_number[i],0,carry);
        }
        else
        {
          m_number[i]=detail::add_single_number(m_number[i],other.m_number[i],carry);
        }
      }
      if (carry>0)
      {
        m_number.push_back(carry);
      }
      is_well_defined();
    }


    /* \brief Standard addition operator.
     */
    big_natural_number operator+(const big_natural_number& other) const
    {
      is_well_defined();
      other.is_well_defined();

      big_natural_number result=*this;
      result.add(other);
      return result;
    }

    /* \brief Standard subtraction. 
       \detail Subtract other from this number. Throws an exception if 
               the result is negative and cannot be represented. 
    */
    void subtract(const big_natural_number& other)
    {
      is_well_defined();
      other.is_well_defined();

      assert(m_number.size()==0 || m_number.back()!=0);
      assert(other.m_number.size()==0 || other.m_number.back()!=0);
      if (m_number.size()<other.m_number.size())
      {
        throw mcrl2::runtime_error("Subtracting numbers " + pp(*this) + " and " + pp(other) + " yields a non representable negative result (1).");
      }
      std::size_t carry=0;
      for(std::size_t i=0; i<m_number.size(); ++i)
      {
        if (i>=other.m_number.size())
        {
          m_number[i]=detail::subtract_single_number(m_number[i],0,carry);
        }
        else
        {
          m_number[i]=detail::subtract_single_number(m_number[i],other.m_number[i],carry);
        }
      }
      if (carry>0)
      {
        throw mcrl2::runtime_error("Subtracting numbers " + pp(*this) + " and " + pp(other) + " yields a non representable negative result (2).");
      }
      remove_significant_digits_that_are_zero();
      is_well_defined();
    } 

    /* \brief Standard subtraction operator. Throws an exception if the result
     *        is negative and cannot be represented.
     */
    big_natural_number operator-(const big_natural_number& other) const
    {
      is_well_defined();
      other.is_well_defined();
      big_natural_number result=*this;
      result.subtract(other);
      return result;
    }

    

    /* \brief Efficient multiplication operator that does not declare auxiliary vectors.
       \detail Initially result must be zero. At the end: result equals (*this)*other+result.
               The calculation_buffer does not need to be initialised. 
     */
    void multiply(const big_natural_number& other,
                  big_natural_number& result,
                  big_natural_number& calculation_buffer_for_multiplicand) const
    {
      is_well_defined();
      other.is_well_defined();
      std::size_t offset=0; 
      for(std::size_t digit: other.m_number)
      {
        // Move n2 to the multiplicand and multiply it with base^offset.
        calculation_buffer_for_multiplicand.clear();
        for(std::size_t i=0; i< offset; ++i) { calculation_buffer_for_multiplicand.m_number.push_back(0); }
        for(std::size_t n: m_number)
        { 
          calculation_buffer_for_multiplicand.m_number.push_back(n);
        }

        // Multiply the multiplicand with digit.
        calculation_buffer_for_multiplicand.multiply_by(digit,0);
        // Add the multiplicand to the result.
        result.add(calculation_buffer_for_multiplicand);
        offset++;
      }
      result.is_well_defined();
    }

    /* \brief Standard multiplication operator. This is not very efficient as it declares two temporary vectors.
     */
    big_natural_number operator*(const big_natural_number& other) const
    {
      big_natural_number result;
      big_natural_number buffer;
      multiply(other,result,buffer);
      return result;
    } 

    // This is an auxiliary function getting the n-th digit,
    // where digit 0 is the least significant one.
    // It is not necessary that n is in range.
    /* std::size_t getdigit(const std::size_t n) const
    {
      if (n>=m_number.size()) 
      {
        return 0;
      }
      return m_number[n];
    } */

    /* \brief Efficient divide operator that does not declare auxiliary vectors.
       \detail Initially result must be zero. At the end: result equals (*this)*other+result.
               The calculation_buffer does not need to be initialised. 
               The algorithm uses standard "primary school" division, except that
               the digits in this case are 64 bits numbers. The calculation is tricky
               as the most significant digit of this may be a remaining digit from
               the previous calculation.
     */
    void div_mod(const big_natural_number& other,
                 big_natural_number& result,
                 big_natural_number& remainder,
                 big_natural_number& calculation_buffer_divisor) const
    {
      is_well_defined();
      other.is_well_defined();
      assert(!other.is_zero());

      if (m_number.size()==1 && other.m_number.size()==1)
      {
        std::size_t n=m_number.front()/other.m_number.front();      // Calculate div.
        if (n==0)
        {
          result.clear();
        }
        else 
        {
          result.m_number.resize(1);
          result.m_number[0]=n;
        }
        n=m_number.front() % other.m_number.front(); // Calculate mod. 
        if (n==0)
        {
          remainder.clear();
        }
        else
        {
          remainder.m_number.resize(1);
          remainder.m_number[0]=n;
        }
        result.is_well_defined();
        remainder.is_well_defined();
        return;
      }

      // TODO: The procedure below works bitwise, as no efficient division algorithm has yet
      // been implemented. A natural candidate is the algorithm in "Per Brinch Hansen, Multiple
      // length division revisited: A tour of the minefield. Software practice and experience 24,
      // 579-601, 1994.
      result.clear();
      remainder=*this;

      if (m_number.size()<other.m_number.size())
      {
        result.is_well_defined();
        remainder.is_well_defined();
        return; 
      }
      const int no_of_bits_per_digit=std::numeric_limits<std::size_t>::digits;

      big_natural_number divisor;
      calculation_buffer_divisor.clear();
      for(std::size_t i=0; i< (1+m_number.size())-other.m_number.size(); ++i) { calculation_buffer_divisor.m_number.push_back(0); }

      // Place 0 digits at least significant position of the calculation_buffer_divisor to make it of comparable length as the remainder.
      for(std::size_t i: other.m_number)
      {
        calculation_buffer_divisor.m_number.push_back(i);
      }
      calculation_buffer_divisor.remove_significant_digits_that_are_zero();
      
      for(std::size_t i=0; i<=no_of_bits_per_digit*(m_number.size()-other.m_number.size()+1); ++i)
      {
        if (remainder<calculation_buffer_divisor)
        {
          // We cannot subtract the calculation_buffer_divisor from the remainder.
          result.multiply_by(2,0); // result=result*2
        }
        else
        {
          // We subtract the calculation_buffer_divisor from the remainder.
          result.multiply_by(2,1); // result=result*2 + 1
          remainder.subtract(calculation_buffer_divisor);
        }
        calculation_buffer_divisor.divide_by(2); // Shift the calculation_buffer_divisor one bit to the left.
      }
      
      result.remove_significant_digits_that_are_zero();
      result.is_well_defined();
      remainder.is_well_defined();
    }

    /* void div_mod(const big_natural_number& other,
                 big_natural_number& result,
                 big_natural_number& remainder,
                 big_natural_number& calculation_buffer_subtractor) const
    {
      is_well_defined();
      other.is_well_defined();
      assert(!other.is_zero());

      result.clear();
      remainder=*this;

      if (m_number.size()<other.m_number.size())
      {
        return;
      }
this->print_number("div_mod: this: ");
other.print_number("div_mod: other: ");
      std::size_t remaining_digit=0;
      for(std::size_t n=remainder.m_number.size()-1; n>=other.m_number.size(); n--)
      {
std::cerr << "Wat is n " << n << "\n";
        std::size_t r=remainder.getdigit(n+1);
        std::size_t divisor=detail::divide_single_number(
                                            remainder.getdigit(n),
                                            other.m_number.back(),
                                            r);

result.print_number("div_mod: result: ");
remainder.print_number("div_mod: remainder: ");
std::cerr << "dmwhile: divisor: " << divisor <<"\n";
calculation_buffer_subtractor.is_well_defined();
        calculation_buffer_subtractor.m_number=std::vector<std::size_t>(n-other.m_number.size(),0); Inefficient, want constructie
        for(std::size_t i: other.m_number)
        {
          calculation_buffer_subtractor.m_number.push_back(i);
        }
        calculation_buffer_subtractor.multiply_by(divisor,0);
calculation_buffer_subtractor.print_number("div_mod: subtractor: ");
        if (remainder<calculation_buffer_subtractor)
        {
          divisor=divisor-1;
std::cerr << "dmwhile: divisor adapted: " << divisor <<"\n";
          calculation_buffer_subtractor.m_number=std::vector<std::size_t>(n-other.m_number.size(),0);
          for(std::size_t i: other.m_number)
          {
            calculation_buffer_subtractor.m_number.push_back(i);
          }
          calculation_buffer_subtractor.multiply_by(divisor,0);
        }
        result.m_number.push_back(divisor);
remainder.print_number("remainder before assert");         
calculation_buffer_subtractor.print_number("subtractor before assert");         
        assert(remainder>=calculation_buffer_subtractor);
        std::size_t old_remainder_size=remainder.m_number.size();
        remainder.subtract(calculation_buffer_subtractor);
        remaining_digit=(remaining_digit?
                            remainder.m_number.size()+1==old_remainder_size:
                            remainder.m_number.size()==old_remainder_size);
      }
result.print_number("div_mod: result voor swap: ");
      if (result.m_number.size()==0)
      {
        return;
      }
      // Result must be reverted.
      std::size_t begin = 0; 
      std::size_t end = result.m_number.size()-1;
      while (begin<end)
      {
        std::swap(result.m_number[begin],result.m_number[end]);
        begin++; 
        end--;
      }
      result.remove_significant_digits_that_are_zero();
result.print_number("div_mod: result na swap: ");
      result.is_well_defined();
      remainder.is_well_defined();
    } */

    /* \brief Standard division operator. This is currently implemented by a bit wise subtraction. Can be optimized by a 64 bit calculation.
       \detail. This routine is not particularly efficient as it declares three temporary vectors.
     */
    big_natural_number operator/(const big_natural_number& other) const
    {
      // Division by zero is not allowed.
      if (other.is_zero())
      {
        throw mcrl2::runtime_error("Division by zero.");
      }
      // Zero divided by something is zero.
      if (is_zero())
      {
        return *this;
      }
      
      // Often numbers only consist of one digit. Deal with this using machine division.
      if (m_number.size()==1 && other.m_number.size()==1)
      {
        return big_natural_number(m_number.front()/other.m_number.front());
      }
      
      // Otherwise do a multiple digit division.
      big_natural_number result;
      big_natural_number remainder;
      big_natural_number buffer;
      div_mod(other,result,remainder,buffer);
      return result;
    } 

    /* \brief Standard modulo operator. This is currently implemented by a bit wise subtraction. Can be optimized by a 64 bit calculation.
       \detail. This routine is not particularly efficient as it declares three temporary vectors.
     */
    big_natural_number operator%(const big_natural_number& other) const
    {
      // Modulo zero is yields the value itself. 
      // Zero modulo  something is zero.
      if (other.is_zero() || is_zero())
      {
        return *this;
      }
      
      // Often numbers only consist of one digit. Deal with this using machine division.
      if (m_number.size()==1 && other.m_number.size()==1)
      {
        return big_natural_number(m_number.front()%other.m_number.front());
      }

      big_natural_number result;
      big_natural_number remainder;
      big_natural_number buffer;
      div_mod(other,result,remainder,buffer);
      return remainder;

    } 
};

inline std::ostream& operator<<(std::ostream& ss, const big_natural_number& l)
{
  static big_natural_number n; // This code is not re-entrant, but it avoids declaring a vector continuously. 
  n=l;  
  std::string s; // This string contains the number in reverse ordering.
  for( ; !n.is_zero() ; )
  {
    std::size_t remainder = n.divide_by(10); /* This divides n by 10. */
    s.push_back(static_cast<char>('0'+remainder));
  }
  if (s.empty())
  {
    ss << "0";
  }
  else
  {
    for (char i : std::ranges::reverse_view(s))
    {
      ss << i;
    }
  }
  return ss;
}


/** \brief Standard overload of swap.
 **/
inline void swap(big_natural_number& x, big_natural_number& y) noexcept
{
  x.m_number.swap(y.m_number);
}

} // namespace mcrl2::utilities

namespace std
{

template <>
struct hash< mcrl2::utilities::big_natural_number >
{
  std::size_t operator()(const mcrl2::utilities::big_natural_number& n) const
  {
    hash<std::vector<std::size_t> > hasher;
    return hasher(n.m_number);
  }
};

  
} // namespace std

namespace mcrl2::utilities
{
/* \brief A pretty print operator on action labels, returning it as a string.
*/
inline std::string pp(const big_natural_number& l)
{
  std::stringstream s;
  s << l;
  return s.str();
}
} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_BIG_NUMBERS_H


