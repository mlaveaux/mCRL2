// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/builder.h
/// \brief add your file description here.

#ifndef MCRL2_CORE_BUILDER_H
#define MCRL2_CORE_BUILDER_H

#include <set>
#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm_list.cxx"
  #include "mcrl2/atermpp/concepts.cxx"
#else
  import atermpp;
#endif

namespace mcrl2::core
{

#ifdef MCRL2_DEBUG_EXPRESSION_BUILDER
inline void msg(const std::string& s)
{
  std::cout << "--- " << s << " ---" << std::endl;
}
#else
inline void msg(const std::string&)
{}
#endif

/**
 * \brief expression builder that visits all sub expressions
 *
 * Types:
 *  \arg Derived the type of a derived class, as per CRTP
 *
 **/
template <typename Derived>
struct builder
{
  // Enter object
  template <typename T>
  void enter(const T&)
  {}

  // Leave object
  template <typename T>
  void leave(const T&)
  {}

  template <typename T>
    requires (!std::ranges::range<T>)
  void update(T& x)
  {
    msg("non-container visit");
    T result;
    static_cast<Derived*>(this)->apply(result, x);
    x = result;
  }

  // container visit
  template <std::ranges::range T>
  void update(T& x)
  {
    msg("container visit");
    for (auto& v: x)
    {
      static_cast<Derived*>(this)->update(v);
    }
  }

  // aterm set visit
  template <typename T>
  void update(std::set<T>& x)
  {
    msg("set visit");
    std::set<T> result;
    for (T v: x)
    {
      static_cast<Derived*>(this)->update(v);
      result.insert(v);
    }
    result.swap(x);
  }

  // term_list visit copy
  template <atermpp::IsATerm T, atermpp::IsATerm U>
  void apply(atermpp::term_list<T>& result, const atermpp::term_list<U>& x)
  {
    msg("term_list traversal");
    atermpp::make_term_list<T>(static_cast<atermpp::term_list<T>&>(result), 
                               x.begin(), 
                               x.end(), 
                               [&](T& result, const U& v) { static_cast<Derived*>(this)->apply(result, v); } );
  }
};


// apply a builder without additional template arguments
template <template <class> class Builder>
class apply_builder: public Builder<apply_builder<Builder> >
{
  using super = Builder<apply_builder<Builder>>;

public:

  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;
};

template <template <class> class Builder>
apply_builder<Builder>
make_apply_builder()
{
  return apply_builder<Builder>();
}

// apply a builder with one additional template argument
template <template <class> class Builder, class Arg1>
class apply_builder_arg1: public Builder<apply_builder_arg1<Builder, Arg1> >
{
  using super = Builder<apply_builder_arg1<Builder, Arg1>>;

public:
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  apply_builder_arg1(const Arg1& arg1)
      : super(arg1)
  {}
};

template <template <class> class Builder, class Arg1>
apply_builder_arg1<Builder, Arg1>
make_apply_builder_arg1(const Arg1& arg1)
{
  return apply_builder_arg1<Builder, Arg1>(arg1);
}

// apply a builder with two additional template arguments
template <template <class> class Builder, class Arg1, class Arg2>
class apply_builder_arg2: public Builder<apply_builder_arg2<Builder, Arg1, Arg2> >
{
  using super = Builder<apply_builder_arg2<Builder, Arg1, Arg2>>;

public:
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  apply_builder_arg2(const Arg1& arg1, const Arg2& arg2)
      : super(arg1, arg2)
  {}
};

template <template <class> class Builder, class Arg1, class Arg2>
apply_builder_arg2<Builder, Arg1, Arg2>
make_apply_builder_arg2(const Arg1& arg1, const Arg2& arg2)
{
  return apply_builder_arg2<Builder, Arg1, Arg2>(arg1, arg2);
}


// apply a builder without additional template arguments
template <template <class> class Builder, class Function>
struct update_apply_builder: public Builder<update_apply_builder<Builder, Function> >
{
  using super = Builder<update_apply_builder<Builder, Function>>;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  using argument_type = typename Function::argument_type; ///< The argument_type is required to restrict the overloads of apply(x) to that type.
  using result_type = typename Function::result_type; ///< The argument_type is required to restrict the overloads of apply(x) to that type.
  const Function& f_;
  
  template <class T>
  void apply(T& result, const argument_type& x) // -> decltype(f_(x))
  {
    // static_cast<result_type&>(result) = f_(x);
    result = static_cast<T>(f_(x));
  }

  update_apply_builder(const Function& f)
    : f_(f)
  {}
};

template <template <class> class Builder, class Function>
update_apply_builder<Builder, Function>
make_update_apply_builder(const Function& f)
{
  return update_apply_builder<Builder, Function>(f);
}

// apply a builder with one additional template argument
template <template <class> class Builder, class Function, class Arg1>
class update_apply_builder_arg1: public Builder<update_apply_builder_arg1<Builder, Function, Arg1> >
{
  using super = Builder<update_apply_builder_arg1<Builder, Function, Arg1>>;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;
  
  using argument_type = typename Function::argument_type; ///< The argument_type is required to restrict the overloads of apply(x) to that type.
  using result_type = typename Function::result_type;
  const Function& f_;
  
  template <class T>
  void apply(T& result, const argument_type& x) // -> decltype(f_(x))
  {
    result = f_(x);
  }

  update_apply_builder_arg1(const Function& f, const Arg1& arg1):
    super(arg1),
    f_(f)
  {}
};

template <template <class> class Builder, class Function, class Arg1>
update_apply_builder_arg1<Builder, Function, Arg1>
make_update_apply_builder_arg1(const Function& f)
{
  return update_apply_builder_arg1<Builder, Function, Arg1>(f);
}

} // namespace mcrl2::core



#endif // MCRL2_CORE_BUILDER_H
