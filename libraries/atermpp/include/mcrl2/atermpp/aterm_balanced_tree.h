// Author(s): Jeroen van der Wulp
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_BALANCED_TREE_H
#define MCRL2_ATERMPP_ATERM_BALANCED_TREE_H

#include "mcrl2/atermpp/aterm_appl.h"

#include <boost/iterator/iterator_facade.hpp>

namespace atermpp
{

namespace
{
  global_function_symbol g_empty("@empty@", 0);
  global_function_symbol g_tree_node("@node@", 2);
  aterm_appl g_empty_tree(g_empty);
}

/// \brief Read-only balanced binary tree of terms.
template <typename Term>
class term_balanced_tree : public aterm_appl
{
  protected:

    static const atermpp::function_symbol& tree_empty_function() { return g_empty; }
    static const function_symbol& tree_node_function() { return g_tree_node; }
    static const aterm_appl& empty_tree() { return g_empty_tree; }

    static void make_node(aterm& term, const term_balanced_tree& left_tree, const term_balanced_tree& right_tree)
    {
      detail::g_thread_term_pool().create_appl(term, tree_node_function(), left_tree, right_tree);
    }

    template < typename ForwardTraversalIterator, class Transformer >
    void make_tree(aterm& term, ForwardTraversalIterator& p, const std::size_t size, Transformer transformer )
    {
      if (size>1)
      {
        std::size_t left_size = (size + 1) >> 1; // size/2 rounded up.
        term_balanced_tree left_tree;
        make_tree(left_tree, p, left_size, transformer);
        std::size_t right_size = size >> 1; // size/2 rounded down.
        term_balanced_tree right_tree;
        make_tree(right_tree, p, right_size, transformer);

        make_node(term, left_tree, right_tree);
        return;
      }

      if (size==1)
      {
        term = down_cast<aterm_appl>(static_cast<const aterm&>(transformer(*(p++))));
        return;
      }

      assert(size==0);
      term = empty_tree();
    }

    explicit term_balanced_tree(detail::_term_appl* t)
         : term_appl(static_cast<detail::_term_appl*>(t))
    {}

  public:

    /// The type of object, T stored in the term_balanced_tree.
    typedef Term value_type;
    
    /// Pointer to T.
    typedef Term* pointer;

    /// Reference to T.
    typedef Term& reference;

    /// Const reference to T.
    typedef const Term const_reference;

    /// An unsigned integral type.
    typedef std::size_t size_type;

    /// A signed integral type.
    typedef ptrdiff_t difference_type;

    /// \brief Default constructor. Creates an empty tree.
    term_balanced_tree()
      : aterm_appl(empty_tree())
    {}

    /// \brief Copy constructor.
    term_balanced_tree(const term_balanced_tree&) noexcept = default;

    /// \brief Move constructor.
    term_balanced_tree(term_balanced_tree&&) noexcept = default;

    /// \brief Assignment operator.
    term_balanced_tree& operator=(const term_balanced_tree&) noexcept = default;

    /// \brief Move assign operator.
    term_balanced_tree& operator=(term_balanced_tree&&) noexcept = default;

    /// \brief Construction from aterm.
    explicit term_balanced_tree(const aterm& tree) 
       : aterm_appl(down_cast<aterm_appl>(tree))
    {
    }

    /// \brief Creates an term_balanced_tree with a copy of a range.
    /// \param first The start of a range of elements.
    /// \param size The size of the range of elements.
    template<typename ForwardTraversalIterator>
    term_balanced_tree(ForwardTraversalIterator first, const std::size_t size)
    {
        make_tree(*this, first, size, [](const Term& t) { return t; });
    }

    /// \brief Creates an term_balanced_tree with a copy of a range, where a transformer is applied to each term
    ///        before adding it to the tree..
    /// \param[in] first The start of a range of elements.
    /// \param[in] size The size of the range of elements.
    /// \param[in] transformer A class with an operator() that is applied to each term before adding it to the tree.
    template<typename ForwardTraversalIterator, typename Transformer>
    term_balanced_tree(ForwardTraversalIterator first, const std::size_t size, Transformer transformer)
    {
        make_tree(*this, first, size, transformer);
    }

    /// \brief Get the left branch of the tree
    /// \details It is assumed that the tree is a node with a left branch.
    /// \return A reference t the left subtree of the current tree
    const term_balanced_tree<Term>& left_branch() const
    {
      assert(is_node());
      return down_cast<const term_balanced_tree<Term>>(aterm_appl::operator[](0));
    }

    /// \brief Get the left branch of the tree
    /// \details It is assumed that the tree is a node with a left branch.
    /// \return A reference t the left subtree of the current tree
    const term_balanced_tree<Term>& right_branch() const
    {
      assert(is_node());
      return down_cast<const term_balanced_tree<Term>>(aterm_appl::operator[](1));
    }

    /// \brief Element indexing operator.
    /// \param position Index in the tree.
    /// \details This operation behaves linearly with respect to container size,
    ///          because it must calculate the size of the container. The operator
    ///          element_at behaves logarithmically.
    const Term& operator[](std::size_t position) const
    {
      return element_at(position, size());
    } 

    /// \brief Get an element at the indicated position.
    /// \param position The required position
    /// \param size The number of elements in the tree.
    ///                    This is required to make the complexity logarithmic.
    /// \details By providing the size this operation is logarithmic. If a wrong
    ///         size is provided the outcome is not determined. See also operator [].
    /// \return The element at the indicated position.
    const Term& element_at(std::size_t position, std::size_t size) const
    {
      assert(size == this->size());
      assert(position < size);

      if (size>1)
      {
        std::size_t left_size = (size + 1) >> 1;

        return (position < left_size) ?
               left_branch().element_at(position, left_size) :
               right_branch().element_at(position-left_size, size - left_size);
      }

      return vertical_cast<Term>(static_cast<const aterm&>(*this));
    }

    /// \brief Returns the size of the term_balanced_tree.
    /// \details This operator is linear in the size of the balanced tree.
    /// \return The size of the tree.
    size_type size() const
    {
      if (is_node())
      {
        return left_branch().size() + right_branch().size();
      }
      return (empty()) ? 0 : 1;
    }

    /// \brief Returns true if tree is empty.
    /// \return True iff the tree is empty.
    bool empty() const
    {
      return m_term->function() == tree_empty_function();
    }

    /// \brief Returns true iff the tree is a node with a left and right subtree.
    /// \return True iff the tree is a node with a left and right subtree.
    bool is_node() const
    {
      return function() == tree_node_function();
    }

    class iterator: public boost::iterator_facade<
      iterator,                            // Derived
      const Term,                          // Value
      boost::forward_traversal_tag,        // CategoryOrTraversal
      const Term&                          // Reference
      >
    {
      private:
        using Tree = term_balanced_tree<Term>;
    
        friend class boost::iterator_core_access;
    
        static constexpr std::size_t maximal_size_of_stack = 20;      // We assume here that a tree never has more than 2^20 leaves, o
                                                           // equivalently that states consist of not more than 2^20 data_expressions.
        unprotected_aterm m_stack[maximal_size_of_stack];
        std::size_t m_top_of_stack;                             // First element in the stack that is empty.
    
        /// \brief Dereference operator
        /// \return The value that the iterator references
        const Term& dereference() const
        {
          assert(m_top_of_stack > 0);
          return static_cast<const Term&>(m_stack[m_top_of_stack-1]);
        }
    
        /// \brief Equality operator
        bool equal(const iterator& other) const
        {
          if (m_top_of_stack != other.m_top_of_stack)
          {
            return false; 
          }
          
          for(std::size_t i = 0; i < m_top_of_stack; ++i)
          { 
            if (m_stack[i] != other.m_stack[i])
            {
              return false;
            }
          }
          return true;
        }
    
        /// \brief Increments the iterator
        void increment()
        {
          --m_top_of_stack;
          if (m_top_of_stack>0)
          {
            unprotected_aterm current = m_stack[m_top_of_stack-1];
            if (current.function() != Tree::tree_node_function())
            {
              // This subtree is empty.
              return;
            }

            --m_top_of_stack;
            do
            {
              m_stack[m_top_of_stack++] = static_cast<Tree&>(current).right_branch();
              current = static_cast<Tree&>(current).left_branch();
            }
            while (current.function() == Tree::tree_node_function());
    
            m_stack[m_top_of_stack++] = current;
          }
        }
    
        void initialise(const term_balanced_tree<Term>& tree)
        {
          if (tree.empty())
          {
            return;
          }

          unprotected_aterm current = tree;
    
          while (current.function() == Tree::tree_node_function())
          {
            assert(m_top_of_stack + 1 < maximal_size_of_stack);
            m_stack[m_top_of_stack++] = static_cast<Tree&>(current).right_branch();
            current = static_cast<Tree&>(current).left_branch();
          }

          assert(m_top_of_stack + 1 < maximal_size_of_stack);
          m_stack[m_top_of_stack++] = current;
        }
    
      public:
    
        iterator()
          : m_top_of_stack(0)
        { }
    
        iterator(const term_balanced_tree<Term>& tree)
          : m_top_of_stack(0)
        {
          initialise(tree);
        } 
    
        iterator(const iterator& other) 
           : m_top_of_stack(other.m_top_of_stack)
        { 
          for(std::size_t i = 0; i < m_top_of_stack; ++i)
          {
            m_stack[i] = other.m_stack[i];
          }
        }
    
    };

    /// \brief Returns an iterator pointing to the beginning of the term_balanced_tree.
    /// \return The beginning of the list.
    iterator begin() const
    {
      return iterator(*this);
    }

    /// \brief Returns an iterator pointing to the end of the term_balanced_tree.
    /// \return The end of the list.
    iterator end() const
    {
      return iterator();
    }

};

/// \brief A term_balanced_tree with elements of type aterm.
typedef term_balanced_tree<aterm> aterm_balanced_tree;

template <class Term>
std::string pp(const term_balanced_tree<Term> t)
{
  std::stringstream ss;
  for(typename term_balanced_tree<Term>::iterator i = t.begin(); i != t.end(); ++i)
  {
    if (i!=t.begin()) 
    {
      ss << ", ";
    }
    ss << pp(*i);
  }
  return ss.str();
}
} // namespace atermpp

namespace std
{
/// \brief Swaps two balanced trees.
/// \details This operation is more efficient than exchanging terms by an assignment,
///          as swapping does not require to change the protection of terms.
/// \param t1 The first term
/// \param t2 The second term
template <class T>
inline void swap(atermpp::term_balanced_tree<T>& t1, atermpp::term_balanced_tree<T>& t2)
{
  t1.swap(t2);
}

/// \brief Standard hash function.
template<class T>
struct hash<atermpp::term_balanced_tree<T> >
{
  std::size_t operator()(const atermpp::term_balanced_tree<T>& t) const
  {
    return std::hash<atermpp::aterm>()(t);
  }
};
} // namespace std


#endif // MCRL2_ATERMPP_ATERM_BALANCED_TREE_H
