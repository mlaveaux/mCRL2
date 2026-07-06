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


namespace mcrl2::data
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
    using iterator_category = std::forward_iterator_tag;
};

template <>
class prepend_iterator_tag_convertor<std::input_iterator_tag>
{
  public:
    using iterator_category = std::input_iterator_tag;
};

// Iterator for term_appl which prepends a data_expression to a list convertible to data_expressions.
template <class Iterator >
class term_appl_prepend_iterator
{
  private:
    using traits = std::iterator_traits<int*>;

  public:
    // The value_type.
    using value_type = data_expression;
    // The reference type.
    using reference = const data_expression&;
    // The pointer type.
    using pointer = const data_expression*;
    // Difference type
    using difference_type = ptrdiff_t;
    // The iterator category.
    // The iterator category is a forward_iterator, unless Iterator is an input iterator, in which case
    // it is an input iterator.
    using iterator_category = typename prepend_iterator_tag_convertor<traits::iterator_category>::iterator_category;

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
    term_appl_prepend_iterator& operator=(const term_appl_prepend_iterator& other) = default;

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
    std::size_t m_current_index = 0;
    const bool m_skip_head;
  public:
    skip_function_application_to_head(ArgumentConverter& f, const bool skip_head)
        : m_f(f),

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

// The class below transforms a function that is to be applied to
// the arguments of an application into a function that is not applied
// to the head, and only applied to the arguments.
template <class ArgumentConverter>
class skip_function_application_to_head_assignment
{
  protected:
    ArgumentConverter& m_f;
    std::size_t m_current_index = 0;
    const bool m_skip_head;
  public:
    skip_function_application_to_head_assignment(ArgumentConverter& f, const bool skip_head)
        : m_f(f),

          m_skip_head(skip_head)
    {}

    void operator()(data_expression& result, const data_expression& d)
    {
      if (m_skip_head && m_current_index++==0)
      {
        result=d;
        return;
      }
      m_f(result,d);
      return;
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

inline const data_expression& evaluate_lambda_data_expression(const data_expression& t)
{
  return t;
}

template <class TERM>
  requires(std::is_invocable_r_v<void, TERM, data_expression&>)
inline data_expression evaluate_lambda_data_expression(const TERM& t)
{
  data_expression result;
  t(result);
  return result;
}

template <class TERM>
  requires(std::is_invocable_r_v<const data_expression, TERM, void>)
inline data_expression evaluate_lambda_data_expression(const TERM& t)
{
  return t();
}

template <class HEAD, class CONTAINER>
inline bool check_whether_sorts_match(const HEAD& head_lambda, const CONTAINER& l)
{
  data_expression head = evaluate_lambda_data_expression(head_lambda);
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

inline data_expression_list get_arguments()
{
  return data_expression_list();
}

template <class HEAD, class... ARGUMENTS>
inline
data_expression_list get_arguments(const HEAD& h, const ARGUMENTS&... args)
{
  data_expression_list result=get_arguments(args...);
  result.push_front(evaluate_lambda_data_expression(h));
  return result;
}

} // namespace detail

/// \brief An application of a data expression to a number of arguments
class application: public data_expression
{
  public:
    /// \brief Default constructor.
    application()
      : data_expression(atermpp::aterm(core::detail::function_symbol_DataAppl(0)))
    {}

    /// \brief Constructor.
    template <typename... Terms>
    application(const data_expression& head, const data_expression& arg1, const Terms&... other_arguments)
      requires(std::conjunction_v<std::is_convertible<Terms, data_expression>...>)
        : data_expression(atermpp::aterm(core::detail::function_symbol_DataAppl(sizeof...(Terms) + 2),
              head,
              arg1,
              other_arguments...))
    {
      assert(detail::check_whether_sorts_match(head, detail::get_arguments(arg1, other_arguments...)));
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
      requires(std::ranges::forward_range<Container> &&
               std::is_convertible_v<std::ranges::range_value_t<Container>, data_expression>)
    application(const data_expression& head,
                const Container& arguments)
      : data_expression(atermpp::aterm(core::detail::function_symbol_DataAppl(arguments.size() + 1),
                                         detail::term_appl_prepend_iterator<typename Container::const_iterator>(arguments.begin(), &head),
                                         detail::term_appl_prepend_iterator<typename Container::const_iterator>(arguments.end())))
    {
      assert(arguments.size()>0);
      assert(detail::check_whether_sorts_match(head,arguments));
    }

  private:
    // forbid the use of iterator, which is silently inherited from
    // aterm. Modifying the arguments of an application through the iterator
    // is not allowed!
    using iterator = data_expression::iterator;

  public:

    /// \brief An iterator to traverse the arguments of an application.
    /// \details There is a subtle difference with the arguments of an iterator on
    ///          the arguments of an aterm from which an application is derived.
    ///          As an application has a head as its first argument, the iterator
    ///          of the aterm starts at this head, where the iterator of the
    ///          application starts at the first argument. This also means that
    ///          t[n] for t an application is equal to t[n+1] if t is interpreted as an
    ///          aterm.
    using const_iterator = atermpp::term_appl_iterator<data_expression>;

    /// \brief Constructor.
    template <typename FwdIter>
      requires(!std::is_base_of_v<data_expression, FwdIter>)
    application(const data_expression& head,
        FwdIter first,
        FwdIter last)
        : data_expression(atermpp::aterm(core::detail::function_symbol_DataAppl(std::distance(first, last) + 1),
              detail::term_appl_prepend_iterator<FwdIter>(first, &head),
              detail::term_appl_prepend_iterator<FwdIter>(last)))
    {
      assert(first!=last);
      assert(detail::check_whether_sorts_match(head,data_expression_list(begin(), end())));
    }

    /// \brief Constructor.
    template <typename FwdIter>
      requires(!std::is_base_of_v<data_expression, FwdIter>)
    application(const std::size_t arity,
        const data_expression& head,
        FwdIter first,
        FwdIter last)
        : data_expression(atermpp::aterm(core::detail::function_symbol_DataAppl(arity + 1),
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
    /// \param head This is the new head for the application.
    /// \param first This is a forward iterator yielding the first argument.
    /// \param last  This is an iterator beyond the last argument.
    /// \param convert_arguments This is a function applied to optionally the head and the arguments.
    /// \param skip_first_argument A boolean which is true if the function must not be applied to the head.
    template <typename FwdIter, class ArgumentConverter>
      requires(!std::is_base_of_v<data_expression, FwdIter> && !std::is_base_of_v<data_expression, ArgumentConverter>)
    application(const data_expression& head,
        FwdIter first,
        FwdIter last,
        ArgumentConverter convert_arguments,
        const bool skip_first_argument = false)
        : data_expression(atermpp::aterm(core::detail::function_symbol_DataAppl(std::distance(first, last) + 1),
              detail::term_appl_prepend_iterator<FwdIter>(first, &head),
              detail::term_appl_prepend_iterator<FwdIter>(last),
              detail::skip_function_application_to_head(convert_arguments, skip_first_argument)))
    {
      assert(first!=last);
      assert(detail::check_whether_sorts_match(head,data_expression_list(begin(), end())));
    }

    /// \brief Constructor.
    /// \details Construct at term head(arg_first,...,arg_last) where convert_arguments
    ///          has been applied to the head and all the arguments.
    /// \param head This is the new head for the application.
    /// \param first This is a forward iterator yielding the first argument.
    /// \param last  This is an iterator beyond the last argument.
    /// \param convert_arguments This is a function applied to optionally the head and the arguments.
    /// \param skip_first_argument A boolean which is true if the function must not be applied to the head.
    template <typename FwdIter, class ArgumentConverter>
      requires(!std::is_base_of_v<data_expression, FwdIter> &&
               !std::is_base_of_v<data_expression, ArgumentConverter> &&
               std::is_same_v<std::invoke_result_t<ArgumentConverter, data_expression&, typename FwdIter::value_type>,
                   void>)
    application(const data_expression& head,
        FwdIter first,
        FwdIter last,
        ArgumentConverter convert_arguments,
        const bool skip_first_argument = false)

        : data_expression(atermpp::aterm(core::detail::function_symbol_DataAppl(std::distance(first, last) + 1),
              detail::term_appl_prepend_iterator<FwdIter>(first, &head),
              detail::term_appl_prepend_iterator<FwdIter>(last),
              detail::skip_function_application_to_head_assignment(convert_arguments, skip_first_argument)))
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
      return atermpp::down_cast<data_expression>(atermpp::aterm::operator[](0));
    }

    /// \brief Get the i-th argument of this expression.
    const data_expression& operator[](std::size_t index) const
    {
      assert(index<size());
      return atermpp::down_cast<data_expression>(atermpp::aterm::operator[](index+1));
    }

    /// \brief Returns an iterator pointing to the first argument of the
    ///        application.
    const_iterator begin() const
    {
      return atermpp::detail::aterm_appl_iterator_cast<data_expression>(atermpp::aterm::begin()+1);
    }

    /// \brief Returns an iterator pointing past the last argument of the
    ///        application.
    const_iterator end() const
    {
      return atermpp::detail::aterm_appl_iterator_cast<data_expression>(atermpp::aterm::end());
    }

    /// \return The number of arguments of this application.
    std::size_t size() const
    {
      return atermpp::aterm::size() - 1;
    }
};

/// \brief swap overload
inline void swap(application& t1, application& t2) noexcept
{
  t1.swap(t2);
}

/// \brief Make function for an application.
/// \param result variable into which the application is constructed.
inline void make_application(atermpp::aterm& result)
{
  atermpp::make_term_appl(result,core::detail::function_symbol_DataAppl(1));
}

/// \brief Constructor.
/// \param result variable into which the application is constructed.
template <typename HEAD,
    typename TERM,
    typename... Terms>
  requires(std::disjunction_v<typename std::is_convertible<HEAD, data_expression>,
        typename std::is_invocable_r<void, HEAD, data_expression&>,
        typename std::is_invocable_r<const data_expression, HEAD, void>> &&
      std::conjunction_v<std::disjunction<typename std::is_convertible<Terms, data_expression>,
          typename std::is_invocable_r<void, Terms, data_expression&>,
          typename std::is_invocable_r<const data_expression, Terms, void>>...>)
inline void
make_application(atermpp::aterm& result, const HEAD& head, const TERM& arg1, const Terms&... other_arguments)
{
  assert(detail::check_whether_sorts_match(head, detail::get_arguments(arg1, other_arguments...)));
  atermpp::make_term_appl(result,
                          core::detail::function_symbol_DataAppl(sizeof...(Terms)+2),
                          head,
                          arg1,
                          other_arguments...);
}

/// \brief Constructor.
/// \param result variable into which the application is constructed.
template <typename Container>
  requires(std::ranges::forward_range<Container> &&
           std::is_convertible_v<std::ranges::range_value_t<Container>, data_expression>)
inline void make_application(data_expression& result,
                      const data_expression& head,
                      const Container& arguments)
{
  assert(arguments.size()>0);
  assert(detail::check_whether_sorts_match(head,arguments));
  atermpp::make_term_appl(result,
                          core::detail::function_symbol_DataAppl(arguments.size() + 1),
                          detail::term_appl_prepend_iterator<typename Container::const_iterator>(arguments.begin(), &head),
                          detail::term_appl_prepend_iterator<typename Container::const_iterator>(arguments.end()));
}

/// \brief Constructor.
/// \param result variable into which the application is constructed.
template <typename FwdIter>
  requires(!std::is_base_of_v<data_expression, FwdIter>)
inline void make_application(atermpp::aterm& result,
    const data_expression& head,
    FwdIter first,
    FwdIter last)
{
  assert(first!=last);
  atermpp::make_term_appl(result,
                          core::detail::function_symbol_DataAppl(std::distance(first, last) + 1),
                          detail::term_appl_prepend_iterator<FwdIter>(first, &head),
                          detail::term_appl_prepend_iterator<FwdIter>(last));
  assert(detail::check_whether_sorts_match(head,data_expression_list(
                                                   atermpp::down_cast<application>(result).begin(),
                                                   atermpp::down_cast<application>(result).end())));
}

/// \brief Constructor.
/// \param result variable into which the application is constructed.
template <typename FwdIter>
  requires(!std::is_base_of_v<data_expression, FwdIter>)
inline void make_application(atermpp::aterm& result,
    const std::size_t arity,
    const data_expression& head,
    FwdIter first,
    FwdIter last)
{
  assert(arity>0);
  assert(std::distance(first, last)==arity);
  atermpp::make_term_appl(result,
                          core::detail::function_symbol_DataAppl(arity + 1),
                          detail::term_appl_prepend_iterator<FwdIter>(first, &head),
                          detail::term_appl_prepend_iterator<FwdIter>(last));
  assert(detail::check_whether_sorts_match(head,data_expression_list(
                                                   atermpp::down_cast<application>(result).begin(),
                                                   atermpp::down_cast<application>(result).end())));
}


/// \brief Constructor.
/// \details Construct at term head(arg_first,...,arg_last) where convert_arguments
///          has been applied to the head and all the arguments.
/// \param result variable into which the application is constructed.
/// \param head This is the new head for the application.
/// \param first This is a forward iterator yielding the first argument.
/// \param last  This is an iterator beyond the last argument.
/// \param convert_arguments This is a function applied to optionally the head and the arguments.
/// \param skip_first_argument A boolean which is true if the function must not be applied to the head.
template <typename FwdIter, class ArgumentConverter>
  requires(!std::is_base_of_v<data_expression, FwdIter> &&
           !std::is_base_of_v<data_expression, ArgumentConverter> &&
           std::is_same_v<std::invoke_result_t<ArgumentConverter, typename FwdIter::value_type>,
               data_expression>)
inline void make_application(atermpp::aterm& result,
    const data_expression& head,
    FwdIter first,
    FwdIter last,
    ArgumentConverter convert_arguments,
    const bool skip_first_argument = false)
{
  assert(first!=last);
  atermpp::make_term_appl(result,
                          core::detail::function_symbol_DataAppl(std::distance(first, last) + 1),
                          detail::term_appl_prepend_iterator<FwdIter>(first, &head),
                          detail::term_appl_prepend_iterator<FwdIter>(last),
                          detail::skip_function_application_to_head(convert_arguments,skip_first_argument));
  assert(detail::check_whether_sorts_match(head,data_expression_list(
                                                   atermpp::down_cast<application>(result).begin(),
                                                   atermpp::down_cast<application>(result).end())));
}

/// \brief Constructor.
/// \details Construct at term head(arg_first,...,arg_last) where convert_arguments
///          has been applied to the head and all the arguments.
/// \param result variable into which the application is constructed.
/// \param head This is the new head for the application.
/// \param first This is a forward iterator yielding the first argument.
/// \param last  This is an iterator beyond the last argument.
/// \param convert_arguments This is a function applied to optionally the head and the arguments.
/// \param skip_first_argument A boolean which is true if the function must not be applied to the head.
template <typename FwdIter, class ArgumentConverter>
  requires(!std::is_base_of_v<data_expression, FwdIter> &&
           !std::is_base_of_v<data_expression, ArgumentConverter> &&
           std::is_same_v<std::invoke_result_t<ArgumentConverter, data_expression&, typename FwdIter::value_type>,
               void>)
static inline void make_application(atermpp::aterm& result,
    const data_expression& head,
    FwdIter first,
    FwdIter last,
    ArgumentConverter convert_arguments,
    const bool skip_first_argument = false)
{
  atermpp::make_term_appl(
                          result,
                          core::detail::function_symbol_DataAppl(std::distance(first, last) + 1),
                          detail::term_appl_prepend_iterator<FwdIter>(first, &head),
                          detail::term_appl_prepend_iterator<FwdIter>(last),
                          detail::skip_function_application_to_head_assignment(convert_arguments,skip_first_argument));
}


//--- start generated class application ---//
// prototype declaration
std::string pp(const application& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
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

} // namespace mcrl2::data



#endif // MCRL2_DATA_APPLICATION_H
