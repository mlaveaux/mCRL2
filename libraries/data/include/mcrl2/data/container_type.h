// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/container.h
/// \brief The class container_type.

#ifndef MCRL2_DATA_CONTAINER_TYPE_H
#define MCRL2_DATA_CONTAINER_TYPE_H

#include "mcrl2/core/detail/default_values.h"
#include "mcrl2/core/detail/soundness_checks.h"

namespace mcrl2::data
{

//--- start generated classes ---//
/// \\brief Container type
class container_type: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    container_type()
      : atermpp::aterm(core::detail::default_values::SortConsType)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit container_type(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_SortConsType(*this));
    }

    /// Move semantics
    container_type(const container_type&) noexcept = default;
    container_type(container_type&&) noexcept = default;
    container_type& operator=(const container_type&) noexcept = default;
    container_type& operator=(container_type&&) noexcept = default;
};

/// \\brief list of container_types
using container_type_list = atermpp::term_list<container_type>;

/// \\brief vector of container_types
using container_type_vector = std::vector<container_type>;

// prototype declaration
std::string pp(const container_type& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const container_type& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(container_type& t1, container_type& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Container type for lists
class list_container: public container_type
{
  public:
    /// \\brief Default constructor X3.
    list_container()
      : container_type(core::detail::default_values::SortList)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit list_container(const atermpp::aterm& term)
      : container_type(term)
    {
      assert(core::detail::check_term_SortList(*this));
    }

    /// Move semantics
    list_container(const list_container&) noexcept = default;
    list_container(list_container&&) noexcept = default;
    list_container& operator=(const list_container&) noexcept = default;
    list_container& operator=(list_container&&) noexcept = default;
};

/// \\brief Test for a list_container expression
/// \\param x A term
/// \\return True if \\a x is a list_container expression
inline
bool is_list_container(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortList;
}

// prototype declaration
std::string pp(const list_container& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const list_container& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(list_container& t1, list_container& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Container type for sets
class set_container: public container_type
{
  public:
    /// \\brief Default constructor X3.
    set_container()
      : container_type(core::detail::default_values::SortSet)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit set_container(const atermpp::aterm& term)
      : container_type(term)
    {
      assert(core::detail::check_term_SortSet(*this));
    }

    /// Move semantics
    set_container(const set_container&) noexcept = default;
    set_container(set_container&&) noexcept = default;
    set_container& operator=(const set_container&) noexcept = default;
    set_container& operator=(set_container&&) noexcept = default;
};

/// \\brief Test for a set_container expression
/// \\param x A term
/// \\return True if \\a x is a set_container expression
inline
bool is_set_container(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortSet;
}

// prototype declaration
std::string pp(const set_container& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const set_container& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(set_container& t1, set_container& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Container type for bags
class bag_container: public container_type
{
  public:
    /// \\brief Default constructor X3.
    bag_container()
      : container_type(core::detail::default_values::SortBag)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit bag_container(const atermpp::aterm& term)
      : container_type(term)
    {
      assert(core::detail::check_term_SortBag(*this));
    }

    /// Move semantics
    bag_container(const bag_container&) noexcept = default;
    bag_container(bag_container&&) noexcept = default;
    bag_container& operator=(const bag_container&) noexcept = default;
    bag_container& operator=(bag_container&&) noexcept = default;
};

/// \\brief Test for a bag_container expression
/// \\param x A term
/// \\return True if \\a x is a bag_container expression
inline
bool is_bag_container(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortBag;
}

// prototype declaration
std::string pp(const bag_container& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const bag_container& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(bag_container& t1, bag_container& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Container type for finite sets
class fset_container: public container_type
{
  public:
    /// \\brief Default constructor X3.
    fset_container()
      : container_type(core::detail::default_values::SortFSet)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit fset_container(const atermpp::aterm& term)
      : container_type(term)
    {
      assert(core::detail::check_term_SortFSet(*this));
    }

    /// Move semantics
    fset_container(const fset_container&) noexcept = default;
    fset_container(fset_container&&) noexcept = default;
    fset_container& operator=(const fset_container&) noexcept = default;
    fset_container& operator=(fset_container&&) noexcept = default;
};

/// \\brief Test for a fset_container expression
/// \\param x A term
/// \\return True if \\a x is a fset_container expression
inline
bool is_fset_container(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortFSet;
}

// prototype declaration
std::string pp(const fset_container& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const fset_container& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(fset_container& t1, fset_container& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Container type for finite bags
class fbag_container: public container_type
{
  public:
    /// \\brief Default constructor X3.
    fbag_container()
      : container_type(core::detail::default_values::SortFBag)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit fbag_container(const atermpp::aterm& term)
      : container_type(term)
    {
      assert(core::detail::check_term_SortFBag(*this));
    }

    /// Move semantics
    fbag_container(const fbag_container&) noexcept = default;
    fbag_container(fbag_container&&) noexcept = default;
    fbag_container& operator=(const fbag_container&) noexcept = default;
    fbag_container& operator=(fbag_container&&) noexcept = default;
};

/// \\brief Test for a fbag_container expression
/// \\param x A term
/// \\return True if \\a x is a fbag_container expression
inline
bool is_fbag_container(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortFBag;
}

// prototype declaration
std::string pp(const fbag_container& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const fbag_container& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(fbag_container& t1, fbag_container& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated classes ---//

} // namespace mcrl2::data



namespace std
{
  
template<>
struct hash<mcrl2::data::container_type>
{
    std::size_t operator()(const mcrl2::data::container_type& v) const
    {
      const hash<atermpp::aterm> hasher;
      return hasher(v);
    }
};

} // namespace std
  


#endif // MCRL2_DATA_CONTAINER_TYPE_H
