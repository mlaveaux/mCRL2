// Author(s): Jeroen Keiren, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/application.h
/// \brief The class application.

#ifndef MCRL2_DATA_VARIABLE_H
#include "mcrl2/data/variable.h"
#endif

#ifndef MCRL2_DATA_APPLICATION_H
#define MCRL2_DATA_APPLICATION_H

// #include "mcrl2/atermpp/aterm_list.h"

namespace mcrl2
{

namespace data
{

namespace detail
{

// The class and the specialisation below are intended to
// make the type of the term_appl_prepend_iterator a forward iterator,
// unless it is based on an input iterator, in which case it should be an
// input iterator.
template <class IteratorTag>
class prepend_iterator_tag_convertor
{
  public:
    typedef std::forward_iterator_tag iterator_category;
};

template <>
class prepend_iterator_tag_convertor<std::input_iterator_tag>
{
  public:
    typedef std::input_iterator_tag iterator_category;
};

// Iterator for term_appl which prepends a data_expression to a list convertible to data_expressions.
template <class Iterator >
class term_appl_prepend_iterator
{
  private:
    typedef std::iterator_traits<int*> traits;

  public:
    // The value_type.
    typedef data_expression value_type;
    // The reference type.
    typedef const data_expression& reference;
    // The pointer type.
    typedef const data_expression* pointer;
    // Difference type
    typedef ptrdiff_t difference_type;
    // The iterator category.
    // The iterator category is a forward_iterator, unless Iterator is an input iterator, in which case
    // it is an input iterator.
    typedef typename prepend_iterator_tag_convertor<traits::iterator_category>::iterator_category iterator_category;

  protected:
    Iterator m_it;
    pointer m_prepend;

  private:
    // Prevent the use of the following operators in this class, including the
    // postfix increment.
    reference operator[](difference_type n) const;
    term_appl_prepend_iterator operator++(int);
    term_appl_prepend_iterator& operator--();
    term_appl_prepend_iterator& operator--(int);
    term_appl_prepend_iterator& operator+=(difference_type n);
    term_appl_prepend_iterator& operator-=(difference_type n);
    term_appl_prepend_iterator operator+(difference_type n) const;
    term_appl_prepend_iterator operator-(difference_type n) const;
    difference_type operator-(const term_appl_prepend_iterator& other) const;
    difference_type distance_to(const term_appl_prepend_iterator& other) const;

  public:
    // Constructor.
    term_appl_prepend_iterator(Iterator it,
                               pointer prepend=nullptr)
      : m_it(it), m_prepend(prepend)
    {}

    // The copy constructor.
    term_appl_prepend_iterator(const term_appl_prepend_iterator& other)
      : m_it(other.m_it),
        m_prepend(other.m_prepend)
    {
    }

    // The assignment operator.
    // other The term to be assigned.
    // Returns a reference to the assigned iterator.
    term_appl_prepend_iterator& operator=(const term_appl_prepend_iterator& other)
    {
      m_it=other.m_it;
      m_prepend=other.m_prepend;
      return *this;
    }

    // The dereference operator.
    // Return the dereferenced term.
    reference operator*() const
    {
      if (m_prepend)
      {
        return *m_prepend;
      }
      return *(this->m_it);
    }

    // Dereference the current iterator.
    pointer operator->() const
    {
      if (m_prepend)
      {
        return m_prepend;
      }
      return &*(this->m_it);
    }

    // Prefix increment.
    // Returns the iterator after it is incremented.
    term_appl_prepend_iterator& operator++()
    {
      if (m_prepend)
      {
        m_prepend = nullptr;
      }
      else
      {
        ++(this->m_it);
      }
      return *this;
    }

    // Equality of iterators.
    bool operator ==(const term_appl_prepend_iterator& other) const
    {
      return m_prepend==other.m_prepend && this->m_it==other.m_it;
    }

    // Inequality of iterators.
    bool operator !=(const term_appl_prepend_iterator& other) const
    {
      return !(*this==other);
    }

    // Comparison of iterators.
    bool operator <(const term_appl_prepend_iterator& other) const
    {
      return m_prepend < other.m_prepend || (m_prepend==other.m_prepend && this->m_it<other.m_it);
    }

    // Comparison of iterators.
    bool operator <=(const term_appl_prepend_iterator& other) const
    {
      return m_prepend < other.m_prepend || (m_prepend==other.m_prepend && this->m_it<=other.m_it);
    }

    // Comparison of iterators.
    bool operator >(const term_appl_prepend_iterator& other) const
    {
      return other<*this;
    }

    // Comparison of iterators.
    bool operator >=(const term_appl_prepend_iterator& other) const
    {
      return other<=*this;
    }
};

// The class below transforms a function that is to be applied to
// the arguments of an application into a function that is not applied
// to the head, and only applied to the arguments. 
template <class ArgumentConverter>
class skip_function_application_to_head
{
  protected:
    ArgumentConverter& m_f;
    std::size_t m_current_index;
    const bool m_skip_head;
  public:
    skip_function_application_to_head(ArgumentConverter&f, const bool skip_head )
     : m_f(f),
       m_current_index(0),
       m_skip_head(skip_head)
    {}

    data_expression operator()(const data_expression& d)
    {
      if (m_skip_head && m_current_index++==0)
      {
        return d;
      }
      return m_f(d);
    }
};

inline bool contains_untyped_sort(const sort_expression& s)
{
  if (is_untyped_sort(s))
  {
    return true;
  }
  if (is_function_sort(s))
  {
    const function_sort& fs=atermpp::down_cast<function_sort>(s);
    if (contains_untyped_sort(fs.codomain()))
    {
      return true;
    }
    for(const sort_expression& sd: fs.domain())
    {
      if (contains_untyped_sort(sd))
      {
        return true;
      }
    }
  }
  if (is_container_sort(s))
  {
    if (contains_untyped_sort(atermpp::down_cast<container_sort>(s).element_sort()))
    {
      return true;
    }
  }
  assert(is_structured_sort(s) || is_basic_sort(s));
  return true;
}

template <class CONTAINER>
inline bool check_whether_sorts_match(const data_expression& head, const CONTAINER& l)
{
  if (contains_untyped_sort(head.sort()))
  {
    // Most likely head is a just parsed, untyped object.
    return true;
  }
  function_sort fs(head.sort());
  if ((fs.domain().size()==1 && contains_untyped_sort(fs.domain().front())) ||
      (l.size()==1 && contains_untyped_sort(l.front().sort())))
  {
    // This is most likely an application from or to an Rewritten@@term, used in the jitty rewriter
    // to indicate that a term is in normal form.
    return true;
  }
  // Check that the sorts of the function domain of the head matches those of the arguments.
  if (fs.domain().size()!=l.size())
  {
    return false;
  }
  typename CONTAINER::const_iterator i=l.begin();
  for(const sort_expression& s: fs.domain())
  {
    if (!contains_untyped_sort(i->sort()) && s!=i->sort())
    {
      return false;
    }
    ++i;
  }
  return true;
}

} // namespace detail

/// \brief An application of a data expression to a number of arguments
class application: public data_expression
{
  public:
    /// \brief Default constructor.
    application()
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(0)))
    {}

    /// \brief Constructor.
    application(const data_expression& head,
                const data_expression& arg1)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(2),head,arg1))
    {
      assert(detail::check_whether_sorts_match<data_expression_list>(head, {arg1}));
    }

    /// \brief Constructor.
    application(const data_expression& head,
                const data_expression& arg1,
                const data_expression& arg2)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(3),head,arg1,arg2))
    {
      assert(detail::check_whether_sorts_match<data_expression_list>(head, {arg1, arg2}));
    }

    /// \brief Constructor.
    application(const data_expression& head,
                const data_expression& arg1,
                const data_expression& arg2,
                const data_expression& arg3)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(4),head,arg1,arg2,arg3))
    {
      assert(detail::check_whether_sorts_match<data_expression_list>(head, {arg1, arg2, arg3}));
    }

    /// \brief Constructor.
    application(const data_expression& head,
                const data_expression& arg1,
                const data_expression& arg2,
                const data_expression& arg3,
                const data_expression& arg4)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(5),head,arg1,arg2,arg3,arg4))
    {
      assert(detail::check_whether_sorts_match<data_expression_list>(head, {arg1, arg2, arg3, arg4}));
    }

    /// \brief Constructor
    application(const data_expression& head,
                const data_expression& arg1,
                const data_expression& arg2,
                const data_expression& arg3,
                const data_expression& arg4,
                const data_expression& arg5)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(6),head,arg1,arg2,arg3,arg4,arg5))
    {
      assert(detail::check_whether_sorts_match<data_expression_list>(head, {arg1, arg2, arg3, arg4, arg5}));
    }

    /// \brief Constructor
    application(const data_expression& head,
                const data_expression& arg1,
                const data_expression& arg2,
                const data_expression& arg3,
                const data_expression& arg4,
                const data_expression& arg5,
                const data_expression& arg6)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(7),head,arg1,arg2,arg3,arg4,arg5,arg6))
    {
      assert(detail::check_whether_sorts_match<data_expression_list>(head, {arg1, arg2, arg3, arg4, arg5, arg6}));
    }

    /// \brief Constructor.
    /// \param term A term
    explicit application(const atermpp::aterm& term)
      : data_expression(term)
    {
      assert(core::detail::check_term_DataAppl(*this));
    }

    /// \brief Constructor.
    template <typename Container>
    application(const data_expression& head,
                const Container& arguments,
                typename atermpp::enable_if_container<Container, data_expression>::type* = nullptr)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(arguments.size() + 1),
                                         detail::term_appl_prepend_iterator<typename Container::const_iterator>(arguments.begin(), &head),
                                         detail::term_appl_prepend_iterator<typename Container::const_iterator>(arguments.end())))
    {
      assert(arguments.size()>0);
      assert(detail::check_whether_sorts_match(head,arguments));
    }

  private:
    // forbid the use of iterator, which is silently inherited from
    // aterm_appl. Modifying the arguments of an application through the iterator
    // is not allowed!
    typedef data_expression::iterator iterator;

  public:

    /// \brief An iterator to traverse the arguments of an application.
    /// \details There is a subtle difference with the arguments of an iterator on
    ///          the arguments of an aterm_appl from which an application is derived.
    ///          As an application has a head as its first argument, the iterator
    ///          of the aterm_appl starts at this head, where the iterator of the
    ///          application starts at the first argument. This also means that
    ///          t[n] for t an application is equal to t[n+1] if t is interpreted as an
    ///          aterm_appl.
    typedef atermpp::term_appl_iterator<data_expression> const_iterator;

    /// \brief Constructor.
    template <typename FwdIter>
    application(const data_expression& head,
                FwdIter first,
                FwdIter last,
                typename std::enable_if< !std::is_base_of<data_expression, FwdIter>::value>::type* = nullptr)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(std::distance(first, last) + 1),
                                         detail::term_appl_prepend_iterator<FwdIter>(first, &head),
                                         detail::term_appl_prepend_iterator<FwdIter>(last)))
    {
      assert(first!=last);
      assert(detail::check_whether_sorts_match(head,data_expression_list(begin(), end())));
    }

    /// \brief Constructor.
    template <typename FwdIter>
    application(const std::size_t arity,
                const data_expression& head,
                FwdIter first,
                FwdIter last,
                typename std::enable_if< !std::is_base_of<data_expression, FwdIter>::value>::type* = 0)
      : data_expression(atermpp::term_appl<aterm>(core::detail::function_symbol_DataAppl(arity + 1),
                                         detail::term_appl_prepend_iterator<FwdIter>(first, &head),
                                         detail::term_appl_prepend_iterator<FwdIter>(last)))
    {
      assert(arity>0);
      assert(std::distance(first, last)==arity);
      assert(detail::check_whether_sorts_match(head,data_expression_list(begin(), end())));
    }


    /// \brief Constructor.
    /// \details Construct at term head(arg_first,...,arg_last) where convert_arguments
    ///          has been applied to the head and all the arguments. 
    /// \parameter head This is the new head for the application.
    /// \parameter first This is a forward iterator yielding the first argument.
    /// \parameter last  This is an iterator beyond the last argument.
    /// \parameter convert_arguments This is a function applied to optionally the head and the arguments.
    /// \parameter skip_first_argument A boolean which is true if the function must not be applied to the head.
    template <typename FwdIter, class ArgumentConverter>
    application(const data_expression& head,
                FwdIter first,
                FwdIter last,
                ArgumentConverter convert_arguments,
                const bool skip_first_argument=false,
                typename std::enable_if< !std::is_base_of<data_expression, FwdIter>::value>::type* = nullptr,
                typename std::enable_if< !std::is_base_of<data_expression, ArgumentConverter>::value>::type* = nullptr)
      : data_expression(atermpp::term_appl<aterm>(
                                 core::detail::function_symbol_DataAppl(std::distance(first, last) + 1),
                                 detail::term_appl_prepend_iterator<FwdIter>(first, &head),
                                 detail::term_appl_prepend_iterator<FwdIter>(last),
                                 detail::skip_function_application_to_head(convert_arguments,skip_first_argument)))
    {
      assert(first!=last);
      assert(detail::check_whether_sorts_match(head,data_expression_list(begin(), end())));
    } 

    /// Move semantics
    application(const application&) noexcept = default;
    application(application&&) noexcept = default;
    application& operator=(const application&) noexcept = default;
    application& operator=(application&&) noexcept = default;

    /// \brief Get the function at the head of this expression.
    const data_expression& head() const
    {
      return atermpp::down_cast<const data_expression>(atermpp::aterm_appl::operator[](0));
    }

    /// \brief Get the i-th argument of this expression.
    const data_expression& operator[](std::size_t index) const
    {
      assert(index<size());
      return atermpp::down_cast<const data_expression>(atermpp::aterm_appl::operator[](index+1));
    }

    /// \brief Returns an iterator pointing to the first argument of the
    ///        application.
    const_iterator begin() const
    {
      return atermpp::detail::aterm_appl_iterator_cast<data_expression>(atermpp::aterm_appl::begin()+1);
    }

    /// \brief Returns an iterator pointing past the last argument of the
    ///        application.
    const_iterator end() const
    {
      return atermpp::detail::aterm_appl_iterator_cast<data_expression>(atermpp::aterm_appl::end());
    }

    /// \return The number of arguments of this application.
    std::size_t size() const
    {
      return atermpp::aterm_appl::size() - 1;
    }
};

/// \brief swap overload
inline void swap(application& t1, application& t2)
{
  t1.swap(t2);
}

//--- start generated class application ---//
// prototype declaration
std::string pp(const application& x);

/// \brief Outputs the object to a stream
/// \param out An output stream
/// \param x Object x
/// \return The output stream
inline
std::ostream& operator<<(std::ostream& out, const application& x)
{
  return out << data::pp(x);
}
//--- end generated class application ---//

inline
const data_expression& unary_operand(const application& x)
{
  return x[0];
}

inline
const data_expression& binary_left(const application& x)
{
  return x[0];
}

inline
const data_expression& binary_right(const application& x)
{
  return x[1];
}

inline
const data_expression& unary_operand1(const data_expression& x)
{
  const application& y = atermpp::down_cast<application>(x);
  return y[0];
}

inline
const data_expression& binary_left1(const data_expression& x)
{
  const application& y = atermpp::down_cast<application>(x);
  return y[0];
}

inline
const data_expression& binary_right1(const data_expression & x)
{
  const application& y = atermpp::down_cast<application>(x);
  return y[1];
}

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_APPLICATION_H
