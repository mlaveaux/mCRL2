Common functionality
====================
The mCRL2 Library contains some algorithms that are applicable to several libraries.
The most prominent example is a traversal framework and a family of find functions
that has been implemented on top of that. This chapter discusses some of these
algorithms.

Traversal functions
-------------------
The mCRL2 Library contains a range of generic traversal functions. There are
different types of travelsals possible:

.. table:: Traversal types

   ================  ========================================================
   Type              Meaning
   ================  ========================================================
   normal traversal  Traverse an object completely
   sort traversal    Traverses only the parts of an object that contain sorts
   ================  ========================================================


Each library defines classes that support these traversal types, namely
``traverser`` and ``builder``. They are the most commonly occurring traversals.
When writing an algorithm that needs to traverse an object, one can usually
implement it using one of the predefined traversal classes.

When deriving from a traverser class, the ``enter`` and ``leave`` functions must be made visible using

.. code-block:: c++

   using super::enter;
   using super::leave;

.. note::

   It is possible to change the way an object is traversed by overriding member
   functions by providing specialisations of ``operator()``, ``enter()``, or
   ``leave()`` in the derived class.

Implementation details
^^^^^^^^^^^^^^^^^^^^^^

The traverser classes have the following structure:

.. code-block:: c++

   template < typename Derived >
   class traverser
   {
     public:
       void operator()(function_symbol const& e);
       void operator()(data_expression const& e);
       ...
   };

   template < typename Derived >
   class sort_traverser : public traverser< Derived >;

The default implementation for ``operator()`` calls the function ``enter`` upon
entering and ``leave`` upon leaving an object. By default these functions have an
empty body, so they can be optimized away by the compiler. The user can override
these functions in a base class to add behavior. If an object has any sub objects,
they can be visited by applying the ``operator()`` function to them. In the example
below this is done with the ``name()`` attribute of a ``variable``.

.. code-block:: c++

   void operator()(variable const& e)
   {
     ...
     static_cast< Derived& >(*this).enter(e);
     static_cast< Derived& >(*this)(e.name());
     static_cast< Derived& >(*this).leave(e);
     ...
   }

Note that the `Curiously recurring template pattern <http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern>`_
needs to be applied in the implementation of ``traverser``, since it needs access
to protected interfaces of derived classes.

The traversal framework can be reused in other libraries. Each library that
needs traversal support defines its own traverser classes in the same way,
extending the base ``traverser`` template with overloads for its specific
types.

Examples
^^^^^^^^
The following are three examples of applications of traversers and builders.

The first example shows a traverser that computes the sum of all integers.
Applied to a vector of integers, this will compute the sum of all integers in
the vector.

.. literalinclude:: examples/traverser1.cpp

The example is extended with different traversers, and again applied to a
vector of integers in the following example.

.. literalinclude:: examples/traverser2.cpp

Traversers are read only, i.e. they do not modify a term. The following example
shows a builder, that increments each integer.

.. literalinclude:: examples/builder1.cpp

Applications
^^^^^^^^^^^^
The Traversal Framework is a collection of classes that support traversal over hierarchical
mCRL2 data structures. Many algorithms are built on top of this library. The most important
examples are the find and replace functions. The goal of the traversal framework is to provide
algorithms with a uniform interface, and to prevent repetitive code for traversing mCRL2 data
structures.

Find functions
""""""""""""""

The mCRL2 Library contains a wide range of algorithms built on top of the traversers and builders.
An overview of the most commonly used traverser based functions is given below.

.. note::

   The find, replace, and other traverser-based functions listed below are **not** part of the
   ``mcrl2::core`` namespace. They are generic patterns, each implemented as a set of overloads
   in the library that owns the relevant type. For example, ``find_free_variables`` on an LPS
   object lives in ``mcrl2::lps``, not ``mcrl2::core``. Always include the correct library header
   and use the matching namespace overload to avoid compilation errors or silent runtime failures.

The mCRL2 Library contains a range of generic find functions. These are:

.. table:: Generic find functions

   =========================  ==============================================
   Function                   Meaning
   =========================  ==============================================
   ``find_basic_sorts``       Finds basic sorts that occur in an object
   ``find_data_expressions``  Finds data expressions that occur in an object
   ``find_free_variables``    Finds free variables that occur in an object
   ``find_identifiers``       Finds identifiers that occur in an object
   ``find_sort_expressions``  Finds sort expressions that occur in an object
   ``find_variables``         Finds variables that occur in an object
   =========================  ==============================================

These functions can be applied to a wide range of objects in the mCRL2 Library,
ranging from ``data_expression`` to ``action_summand``.

.. table:: Find functions for data specifications

   ====================  ===================================================================
   Function              Meaning
   ====================  ===================================================================
   ``find_constructor``  Finds a constructor
   ``find_equations``    Finds equations with a given expression as head on one of its sides
   ``find_mapping``      Finds a mapping
   ``find_sort``         Finds a sort
   ====================  ===================================================================

.. table:: Generic replace functions

   ===================================  ================================================================================================================================
   Function                             Meaning
   ===================================  ================================================================================================================================
   ``replace_sort_expressions``         Applies a substitution to sort expressions that occur in an object. Optionally it can be applied to nested expressions.
   ``replace_variables``                Applies a substitution to variables that occur in an object within the context of a data expression.
   ``replace_free_variables``           Applies a substitution to free variables that occur in an object within the context of a data expression.
   ``replace_data_expressions``         Applies a substitution to data expressions that occur in an object. Optionally it can be applied to nested expressions.
   ``replace_propositional_variables``  Applies a substitution to propositional variable instantiations that occur in an object within the context of a pbes expression.
   ===================================  ================================================================================================================================

.. table:: Other generic traverser based functions

   ===========================  ==========================================================================================
   Function                     Meaning
   ===========================  ==========================================================================================
   ``translate_user_notation``  Applies a transformation to numbers, sets and bags that is needed after parsing.
   ``normalize_sorts``          Brings embedded sorts into normal form.
   ``rewrite``                  Applies a rewriter to embedded expressions, optionally in combination with a substitution.
   ===========================  ==========================================================================================


These functions can be applied to a wide range of objects in the mCRL2 Library, ranging from ``data::data_expression`` to
``lps::action_summand``. Each library has it's own overloads of these functions. For example to find data variables in an
object of the LPS library, the function ``lps::find_variables`` must be used.

.. warning::

  It is important to choose the right overload of generic traverser functions. Failing to do so may
  result in compilation errors, or runtime exceptions.



Generic programming techniques
------------------------------

The file ``mcrl2/data/detail/container_utility.h`` contains two utility functions
``enable_if_container`` and ``disable_if_container`` for determining if a type is a container type.
It can be used to create an interface that does different things for containers and for
other types:

.. code-block:: c++

   template <typename Container>
   void operator()(const Container& t, typename data::detail::enable_if_container<Container>::type* = 0)
   {
     traverse_container(t);
   }

   template <typename Term>
   void operator()(const Term& t, typename data::detail::disable_if_container<Term>::type* = 0)
   {
     traverse(t);
   }

Static polymorphism
-------------------

For the traversal functions in mCRL2 a technique called `static
polymorphism <http://en.wikipedia.org/wiki/Template_metaprogramming#Static_polymorphism>`_ is used.

Below an example is given how static polymorphism can in principle be
used to implement find functions in the mCRL2 Library.
In this case, A models data::detail::traverser, B models lps::detail::traverser
and C models data::detail::find_helper.

.. code-block:: c++

   #include <iostream>

   // base class
   template <typename Derived>
   struct A
   {
     void a()
     {
       std::cout << "A::a()" << std::endl;
     }
   };

   // extended class
   template <typename Derived>
   struct B: public A<Derived>
   {
     void b()
     {
       std::cout << "B::b()" << std::endl;
       static_cast<Derived&>(*this).a();
     }
   };

   // override the a() function
   template <template <class> class T>
   struct override: public T<override<T> >
   {
     void a()
     {
       std::cout << "override::a()" << std::endl;
     }
   };

   int main()
   {
     override<B> f;
     f.a();
     f.b();

     return 0;
   }

Concepts
--------

In the tables below, X is an substitution type, a and b are constant objects of
type X, s are mutable objects of type X, V is X::variable_type, and v is a
constant object of type V. E is X::expression_type, and e is a constant object
of type E.

A class or built-in type X models the Substitution concept if it models Unary
Function and the following expressions are valid and respect the stated
semantics.

.. table:: Substitution Requirements

   ==========  ===========
   Expression  Return Type
   ==========  ===========
   ``s(v)``    E
   ==========  ===========

A class or built-in type X models the Mutable Substitution concept if it models
the Substitution concept and if X::variable_type can be assigned to
X::expression_type. The following expressions are valid and respect the stated
semantics.

.. table:: Mutable Substitution Requirements (in addition to Substitution)

   ============  ===========  ==============================
   Expression    Return Type  Assertion/ Pre-/Post-condition
   ============  ===========  ==============================
   ``s[v] = e``               post: s(v) == e
   ``s[v] = v``               post: s(v) == v
   ============  ===========  ==============================

The ``mutable_substition`` template class implements the Mutable Substitution
concept. In addition, it has a number additional member functions in common
with a ``std::map``.
