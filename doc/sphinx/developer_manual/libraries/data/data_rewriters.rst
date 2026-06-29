Data rewriters
##############

A rewriter is a function that rewrites terms using a number of rewrite rules.
In the mCRL2 toolset a class ``data::rewriter`` is available that operates on data expressions,
and that is initialized using a data specification. The equations of the data specification
are interpreted as rewrite rules from left to right. An example is given below.

.. code-block:: c++

   // rewrite two data expressions, and check if they are the same
   rewriter r;
   data_expression d1 = parse_data_expression("2+7");
   data_expression d2 = parse_data_expression("4+5");
   assert(d1 != d2);
   assert(r(d1) == r(d2));

.. note::

   Rewriters can be used to determine equivalence between data expressions. In general this
   problem is undecidable. Only if :math:`r(d_1) = r(d_2)` can one conclude that the
   expressions ``d1`` and ``d2`` are equivalent; otherwise the answer is unknown.

For efficiency reasons a rewriter can be invoked with an optional substitution function
:math:`\sigma` as a second argument, where :math:`\sigma` maps data variables to data
expressions. The function :math:`\sigma` must satisfy the property that for all data variables
``v``:

.. math::

   \sigma(v) = r(\sigma(v)).

Under this condition the following property holds:

.. math::

   r(d, \sigma) = r(\sigma(d)).

Note that in general the computation of :math:`r(d, \sigma)` can be done more efficiently than
the computation of :math:`r(\sigma(d))`. In the mCRL2 toolset substitution functions are used
that take constant time.

An example of rewriting with a substitution function is given below.

.. code-block:: c++

   rewriter r;

   // Create a substitution sequence sigma with two substitutions: [m:=3, n:=4]
   std::string var_decl = "m, n: Pos;\n";
   mutable_map_substitution sigma;
   sigma[parse_data_expression("m", var_decl)] = r(parse_data_expression("3"));
   sigma[parse_data_expression("n", var_decl)] = r(parse_data_expression("4"));

   data::data_expression d1 = parse_data_expression("m+n", var_decl);
   data::data_expression d2 = parse_data_expression("7");
   assert(r(d1, sigma) == r(d2));

.. caution::

   The current implementation of rewriting with substitutions to data variables is
   inefficient. The interface of the underlying ``Rewriter`` class needs to be adapted
   to get rid of these inefficiencies.

Rewriter Concept
----------------

In the mCRL2 toolset a ``Rewriter`` is a concept with the following requirements:

.. table:: Associated types

   =======================  ====================================================
   Expression               Meaning
   =======================  ====================================================
   Rewriter::term_type      the type of the terms on which the rewriter operates
   Rewriter::variable_type  the type of the variables
   =======================  ====================================================

A ``SubstitutionFunction`` is a function that maps variables to terms. Let ``sigma`` be a
substitution function, and let ``v`` be an expression of type ``Rewriter::variable_type``.

.. table:: Expression semantics for SubstitutionFunction

   ==========  ===================================================
   Expression  Result
   ==========  ===================================================
   sigma(v)    Returns an expression of type ``Rewriter::term_type``
   ==========  ===================================================

Let ``r`` be a ``Rewriter``, let ``d`` be an expression of type ``Rewriter::term_type``
and let ``sigma`` be a ``SubstitutionFunction``.

.. list-table:: Expression semantics for Rewriter
   :widths: 15 85
   :header-rows: 1

   * - Expression
     - Result
   * - ``r(d)``
     - Returns an expression of type ``Rewriter::term_type`` that is the result
       of applying the rewriter ``r`` to term ``d``.
   * - ``r(d, sigma)``
     - Returns an expression of type ``Rewriter::term_type`` that is the result
       of applying the rewriter ``r`` to term ``d``, while on the fly applying
       the substitution function ``sigma`` to all data variables in ``d``.

The mathematical foundations underlying these operations — formal definitions of
substitutions, rewrite rules, and the abstract rewrite algorithm — are described in
:doc:`Rewriter implementation notes <data_rewriter_implementation>`.

Algorithms using a rewriter
----------------------------

Algorithms that use a rewriter are typically implemented with a template parameter
for the rewriter. An example of this is the constelm algorithm of the LPS library:

.. code-block:: c++

    template <typename Rewriter>
    specification constelm(const specification& spec, Rewriter r, bool verbose = false)
    {
      ...
    }

The algorithm may only assume that the requirements of the Rewriter Concept hold,
with proper choices for the nested variable and term types.

.. seealso::

   :doc:`Rewriter implementation notes <data_rewriter_implementation>` — mathematical
   foundations for the Rewriter concept (HRS theory, match trees, rewriting strategies).

   :doc:`Capture-avoiding substitutions <data_substitutions>` — formal definition of the
   substitution operation used during rewriting.

   :doc:`Data enumerator <data_enumerator>` — an algorithm built on top of the rewriter
   that enumerates values satisfying a predicate.
