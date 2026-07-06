// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/sort_specification.h
/// \brief This file describes the a sort_specification,
//         which contains declared sorts and sort aliases.

#ifndef MCRL2_DATA_SORT_SPECIFICATION_H
#define MCRL2_DATA_SORT_SPECIFICATION_H

#include "mcrl2/data/alias.h"
#include "mcrl2/data/structured_sort.h"

namespace mcrl2::data
{

/// \cond INTERNAL_DOC
namespace detail
{

  template < typename Container, typename T >
  inline
  void remove(Container& container, const T& t)
  {
    typename Container::iterator i = std::find(container.begin(), container.end(), t);
    if(i != container.end())
    {
      container.erase(i);
    }
  }

} // end namespace detail
/// \endcond

class sort_specification;
sort_expression normalize_sorts(const sort_expression& x, const data::sort_specification& sortspec);

class sort_specification
{
  protected:

    /// \brief This boolean indicates whether the variables
    /// m_normalised_constructors, m_mappings, m_equations, m_normalised_sorts,
    /// m_normalised_aliases are up to date with respect to changes of this sort specification.
    mutable bool m_normalised_sorts_are_up_to_date;

    /// \brief The variable below indicates whether a surrounding
    ///        data specification is up to data with respect to
    ///        sort normalisation and available sorts. This is set to false if
    ///        an alias or a new sort is added.
    mutable bool m_normalised_data_is_up_to_date;

    /// \brief The basic sorts and structured sorts in the specification.
    basic_sort_vector m_user_defined_sorts;

    /// \brief Set containing all the sorts, including the system defined ones.
    mutable std::set<sort_expression> m_normalised_sorts;

    /// \brief The sorts that occur are needed in this sort specification but are
    ///        not explicitly defined as user defined sorts. An example is the sort Nat
    ///        when declaring the use of a sort List(Nat).
    ///        The normalised sorts, constructors, mappings and equations are complete
    ///        with respect to these sorts.
    std::set< sort_expression > m_sorts_in_context;

    /// \brief The basic sorts and structured sorts in the specification.
    alias_vector m_user_defined_aliases;

    /// \brief Table containing how sorts should be mapped to normalised sorts.
    mutable std::map< sort_expression, sort_expression > m_normalised_aliases;

  public:

    /// \brief Default constructor
    sort_specification()
     : m_normalised_sorts_are_up_to_date(false),
       m_normalised_data_is_up_to_date(false)
    {
      add_predefined_basic_sorts();
    }

    sort_specification(const basic_sort_vector& sorts, const alias_vector& aliases)
     : m_normalised_sorts_are_up_to_date(false),
       m_normalised_data_is_up_to_date(false)
    {
      add_predefined_basic_sorts();

      for(const basic_sort& sort: sorts)
      {
        add_sort(sort);
      }
      for(const alias& a: aliases)
      {
        add_alias(a);
      }
    }


    /// \brief Adds a sort to this specification
    /// \param[in] s A sort expression.
    void add_sort(const basic_sort& s)
    {
      if (std::find(m_user_defined_sorts.begin(),m_user_defined_sorts.end(),s)==m_user_defined_sorts.end())
      {
        m_user_defined_sorts.push_back(s);
        import_system_defined_sort(s);
        sorts_are_not_necessarily_normalised_anymore();
      }
    }

    /// \brief Adds a sort to this specification, and marks it as system
    ///        defined
    ///
    /// \param[in] s A sort expression.
    /// \pre s does not yet occur in this specification.
    /// \post is_system_defined(s) = true
    /// \note this operation does not invalidate iterators of sorts_const_range
    void add_system_defined_sort(const sort_expression& s)
    {
      import_system_defined_sort(s);
    }

    ///\brief Adds the sort s to the context sorts
    /// \param[in] s a sort expression. It is
    /// added to m_sorts_in_context. For this sort standard functions are generated
    /// automatically (if, <,<=,==,!=,>=,>) and if the sort is a standard sort,
    /// the necessary constructors, mappings and equations are added to the data type.
    void add_context_sort(const sort_expression& s)
    {
      import_system_defined_sort(s);
      sorts_are_not_necessarily_normalised_anymore();
    }

    ///\brief Adds the sorts in c to the context sorts
    /// \param[in] c a container of sort expressions. These are
    /// added to m_sorts_in_context. For these sorts standard functions are generated
    /// automatically (if, <,<=,==,!=,>=,>) and if the sorts are standard sorts,
    /// the necessary constructors, mappings and equations are added to the data type.
    template <typename Container>
    void add_context_sorts(const Container& c, typename atermpp::enable_if_container<Container>::type* = nullptr)
    {
      for(typename Container::const_iterator i=c.begin(); i!=c.end(); ++i)
      {
        add_context_sort(*i);
      }
    }

    /// \brief Removes sort from the user defined sorts in the specification.
    /// Note that this does not remove aliases for the sort, and it does not remove
    /// constructors, mappings and equations, and also keeps the sort as
    /// defined in the context.
    /// \param[in] s A sort expression.
    /// \post s does not occur in this specification.
    void remove_sort(const sort_expression& s)
    {
      detail::remove(m_user_defined_sorts,s);
      sorts_are_not_necessarily_normalised_anymore();
    }

    /// \brief Gets the normalised sort declarations including those that are system defined.
    ///        This is the set with all sorts that occur in a data_specification, or that
    ///        have been registered as sorts used in the context.
    ///
    /// \details The time complexity of this operation is constant, except when
    ///      the data specification has been changed, in which case it can be that
    ///      the sorts must be renormalised. This operation is linear in the size of
    ///      the specification.
    /// \return The sort normalised using the aliases occurring in this specification,
    ///         including the built in sorts such as Bool and Nat, and complex
    ///         sorts that are used in the user defined aliases and sorts, but
    ///         also that are registered as sorts used in the context of the specification.
    inline
    const std::set<sort_expression>& sorts() const
    {
      normalise_sort_specification_if_required();
      return m_normalised_sorts;
    }

    /// \brief Return the user defined context sorts of the current specification.
    /// \details Time complexity is constant.
    /// \return The set with sorts used in the context.
    const std::set<sort_expression>& context_sorts() const
    {
      return m_sorts_in_context;
    }

    /// \brief Gets all sorts defined by a user (excluding the system defined sorts).
    ///
    /// \details The time complexity of this operation is constant.
    /// \return The user defined sort declaration.
    inline
    const basic_sort_vector& user_defined_sorts() const
    {
      return m_user_defined_sorts;
    }

    /// \brief Adds an alias (new name for a sort) to this specification
    /// \param[in] a an alias
    void add_alias(const alias& a)
    {
      if (std::find(m_user_defined_aliases.begin(),m_user_defined_aliases.end(),a)==m_user_defined_aliases.end())
      {
        m_user_defined_aliases.push_back(a);
        import_system_defined_sort(a.name());
        import_system_defined_sort(a.reference());
        sorts_are_not_necessarily_normalised_anymore();
      }
    }


    /// \brief Removes a user defined //alias from specification.
    /// \details This also removes the defined sort of this alias from the set of user defined sorts.
    ///          This routine does not check whether the alias, or name was in the user defined sets.
    void remove_alias(const alias& a)
    {
      detail::remove(m_user_defined_sorts, a.name());
      detail::remove(m_user_defined_aliases, a);
      sorts_are_not_necessarily_normalised_anymore();
    }


    /// \brief Gets the user defined aliases.
    /// \details The time complexity is constant.
    inline
    const alias_vector& user_defined_aliases() const
    {
      return m_user_defined_aliases;
    }


    /// \brief Gets a normalisation mapping that maps each sort to its unique normalised sort
    /// \details Sorts that are mapped to itself are not include in the mapping.
    ///    This map is required in functions with the name normalize_sorts.
    ///    When in a specification sort aliases are used, like sort A=B or
    ///    sort Tree=struct leaf | node(Tree,Tree) then there are different representations
    ///    for each sort. The normalisation mapping maps each sort to a unique representant.
    ///    Moreover, it is this unique sort that it provides in internal mappings.
    const std::map< sort_expression, sort_expression >& sort_alias_map() const
    {
      normalise_sort_specification_if_required();
      return m_normalised_aliases;
    }

    bool operator==(const sort_specification& other) const
    {
      return m_user_defined_sorts==other.m_user_defined_sorts &&
             m_sorts_in_context == other.m_sorts_in_context &&
             m_user_defined_aliases==other.m_user_defined_aliases;
    }

  // Below are auxiliary functions for this class.
  protected:
    void sorts_are_not_necessarily_normalised_anymore() const
    {
      m_normalised_sorts_are_up_to_date=false;
      data_is_not_necessarily_normalised_anymore();
    }

    void data_is_not_necessarily_normalised_anymore() const
    {
      m_normalised_data_is_up_to_date=false;
    }

    void normalise_sort_specification_if_required() const
    {
      if (!m_normalised_sorts_are_up_to_date)
      {
        m_normalised_sorts_are_up_to_date=true;
        m_normalised_sorts.clear();
        reconstruct_m_normalised_aliases();
        for (const sort_expression& s: m_sorts_in_context)
        {
          m_normalised_sorts.insert(normalize_sorts(s,*this));
        }
        for (const sort_expression& s: m_user_defined_sorts)
        {
          m_normalised_sorts.insert(normalize_sorts(s,*this));
        }
        data_is_not_necessarily_normalised_anymore();
      }
    }

    void add_predefined_basic_sorts();

    template <class CONTAINER>
    void import_system_defined_sorts(const CONTAINER& sorts)
    {
      for(const sort_expression& sort: sorts)
      {
        import_system_defined_sort(sort);
      }
    }

    /// \brief Adds the system defined sorts in a sequence.
    ///       The second argument is used to check which sorts are added, to prevent
    ///       useless repetitions of additions of sorts.
    /// \details The function normalise_sorts imports for the given sort_expression sort all sorts, constructors,
    ///       mappings and equations that belong to this sort to the @c normalised sets in this
    ///       data type. E.g. for the sort Nat of natural numbers, it is required that Pos
    ///       (positive numbers) are defined.
    void import_system_defined_sort(const sort_expression& sort);

    // The function below recalculates m_normalised_aliases, such that
    // it forms a confluent terminating rewriting system using which
    // sorts can be normalised.
    void reconstruct_m_normalised_aliases() const;

    // The function below checks whether there is an alias loop, e.g. aliases
    // of the form A=B; B=A; or more complex A=B->C; B=Set(D); D=List(A); Loops
    // through structured sorts are allowed. If a loop is detected, an exception
    // is thrown.
    void check_for_alias_loop(const sort_expression& s,
      std::set<sort_expression> sorts_already_seen,
      bool toplevel = true) const;
}
;

} // namespace mcrl2::data



#endif // MCRL2_DATA_SORT_SPECIFICATION_H
