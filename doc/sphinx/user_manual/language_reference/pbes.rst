.. _language-pbes:

Parameterised Boolean Equation Systems
======================================

Parameterised Boolean equation systems (PBESs) can be used to encode model
checking problems (such as verifying first-order modal :math:`\mu`-calculus
formulae on LPSs, implemented in the tool :ref:`tool-lps2pbes`; the translation
of the first-order modal :math:`\mu`-calculus model checking problem on LPSs is
documented in the `PBES library manual
<../../developer_manual/libraries/pbes/pbes_library.html>`_. Furthermore, PBESs
can be used to code equivalence and preorder relations on processes with data,
see [CPPW07a]_.

PBES expression
---------------

The right hand sides of equations in a PBES are predicate formulae, adhering to
the following syntax.

.. dparser:: PropVarInst PbesExpr

The ``val`` operator stands for the value of a boolean data expression,
``true`` and ``false`` are the booleans true and false, ``!`` is negation,
``&&`` stands for and, ``||`` for or and ``=>`` for implication. The rules
starting with ``forall`` and ``exists`` stand for univeral and existential
quantification.

.. todo:: priorities

The following restrictions apply to propositional variables:

* monotonicity: every occurrence of a propositional variable should be in a
  scope such that the number of ``!`` operators plus the number of left hand
  sides of ``=>`` is even.
* no overloading: it is not allowed to declare two propositional variables with
  the same name but with a different type.

In mathematical notation, a predicate formula is defined as

.. math::

  \varphi ::= b \mid X(e) \mid \varphi \land \varphi \mid \varphi \lor \varphi \mid \forall d \colon D . \varphi \mid \exists d \colon D . \varphi

where :math:`b` is a Boolean expression, :math:`d` is a sorted data variable and
:math:`e` is a data expression of the sort of variable :math:`X`.

PBES equation
-------------

Parameterised Boolean equations are fixed point equations with a propositional
variable declaration as left hand side and a predicate formula as right hand
side. A propositional variable declaration is a sorted predicate variable,
with a finite number of sorted data variables.

.. dparser:: PropVarDecl PbesEqnDecl PbesEqnSpec

In mathematical notation, we write
:math:`(\mu X(d:D) = \varphi)` or :math:`(\nu X(d:D) = \varphi)`
for least and greatest fixpoint equations, where :math:`\varphi` is a
predicate formula.

PBES specification
------------------

A PBES specification contains a sequence of parameterised Boolean equations,
preceded by the ``pbes`` keyword. Furthermore, an initial propositional variable
instantiation must be specified after the ``init`` keyword.

.. dparser:: PbesInit PbesSpec

Files containing a PBES specification can be parsed using :ref:`tool-txt2pbes`.

Transforming PBESs
------------------
Several operations on PBESs can be done without these operations influencing the
solution to the equation system. Such operations include "migration",
"substitution", which form the basis for the so-called Gauß elimination
strategy for solving PBESs. Let :math:`\mathcal{E}`, :math:`\mathcal{F}` and
:math:`\mathcal{G}` denote arbitrary PBESs.
Substitution, for instance is based on the following transformation:

.. math::

  \mathcal{E} (\sigma X(d:D) = \varphi) \mathcal{F} (\sigma' Y(e: E) = \psi) \mathcal{G}

to

.. math::

  \mathcal{E} (\sigma X(d:D) = \varphi[Y := \lambda e: E . \psi]) \mathcal{F} (\sigma' Y(e: E) = \psi) \mathcal{G}

A note of warning: substitution in the other direction (i.e. substituting
:math:`\varphi` for :math:`X` in the equation for :math:`Y`) is *not allowed*
as it affects the solution to the PBES. The PBES library provides the basic
facilities for performing a substitution such as
:math:`\varphi[Y := \lambda e:E. \psi]`, in which every occurrence of :math:`Y`
in :math:`\varphi` is replaced by the predicate :math:`\psi`.

Migration, which is a transformation defined by the following correspondence:

.. math::

  \mathcal{E} (\sigma X(d \colon D) = \varphi) \mathcal{F} \mathcal{G}

to

.. math::

  \mathcal{E} \mathcal{F} (\sigma X(d \colon D) = \varphi) \mathcal{G}

is only allowed when :math:`\varphi` contains *no* predicate variables. Such a
predicate formula is called *simple*, and an equation for which its right-hand
side expression is a simple predicate formula is called *solved*. The PBES
library offers methods to check whether an equation is solved and whether a
predicate formula is simple.

Solving PBESs
-------------
The PBES library provides the means to construct PBESs and modify these. As may
be clear, one is most-often interested in the solution of a PBES, as it provides
the answer to some verification task. There are two main approaces to solving
PBESs:

* Symbolic approximation, combined with Gauß elimination
* Enumerative

Currently, the following strategies have been implemented for solving PBESs:

* Enumerative, by translation to BES, implemented in :ref:`tool-pbessolve`.
* Enumerative, by translation to parity games, implemented in :ref:`tool-pbespgsolve`.
* Symbolic approximation, implemented in :ref:`tool-pbesfixpointsolve`.

Symbolic approximation + Gauß elimination
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As a running example, consider the following PBES:

.. math::

  \mu X(b \colon Bool) & = b \lor X(\neg b) \lor Y(b)\\
  \nu Y(b \colon Bool) & = X(b) \land Y(b)

Gauß Elimination basically employs the migration and substitution
transformations to solve the global PBES, whereas symbolic approximation tries
to solve a single equation by means of an approximation procedure, in which the
approximants are represented by predicate formulae. For instance, the following
sequence of approximations is needed for computing the solution to :math:`Y`:

.. math::

  Y_0 & = true \\
  Y_1 & = (X(b) \land Y(b))[Y := \lambda b \colon Bool. true] \\
      & = X(b) \\
  Y_2 & = (X(b) \land Y(b))[Y := \lambda b \colon Bool. X(b)] \\
      & = X(b) \land X(b) \\
      & = X(b)

Since the approximation process stabilises at the second approximant, the
solution to :math:`Y` is the predicate formula :math:`X(b)`. A solution that is
found by means of approximation can be plugged into the original PBES without
changing the solution to the PBES; in this case, this results in the following
PBES:

.. math::

  \mu X(b \colon Bool) & = b \lor X(\neg b) \lor Y(b) \\
  \nu Y(b \colon Bool) & = X(b)

Substitution then gives the following equivalent PBES:

  \mu X(b \colon Bool) & = b \lor X(\neg b) \lor X(b) \\
  \nu Y(b \colon Bool) & = X(b)

Observe that the equation for :math:`X` is *closed*, meaning that it does not
refer to predicate variables, other than :math:`X`. Solving the equation for
:math:`X` using symbolic approximation, we get:

.. math::

  X_0 & = false \\
  X_1 & = (b \lor X(\neg b) \lor X(b))[ X := \lambda b \colon Bool . false]\\
      & = b\\
  X_2 & = (b \lor X(\neg b) \lor X(b))[ X := \lambda b \colon Bool . b]\\
      & = b \lor \neg b \lor b\\
      & = true

Since there is no predicate formula weaker than :math:`true`, the solution to
:math:`X` is also :math:`true`. Replacing the solution :math:`true` for the
predicate formula in the equation for :math:`X` results in the following
equivalent equation system:

.. math::

  \mu X(b \colon Bool) & = true \\
  \nu Y(b \colon Bool) & = X(b)

Using migration, and, subsequently a substitution, the following **solved** PBES
is obtained:

.. math::

  \nu Y(b \colon Bool) & = true \\
  \mu X(b \colon Bool) & = true

Suppose we would be interested in knowning whether :math:`X(false)` would be
:math:`true` or :math:`false` then requires looking at the solved PBES and
results in the answer :math:`true` for :math:`X(false)`.

Enumerative
^^^^^^^^^^^
Again, as a running example, consider the following PBES:

.. math::

  \mu X(b \colon Bool) & = b \lor X(\neg b) \lor Y(b)\\
  \nu Y(b \colon Bool) & = X(b) \land Y(b)

The enumerative approach explores the equations of a PBES on demand. Suppose we
are interested in knowning whether :math:`X(false)` would be :math:`true` or
:math:`false`. This question can be answered by looking at the equations that
are needed for :math:`X(false)`. This can be found out by the following
procedure:

* replace the data variable :math:`b` with :math:`false` in the predicate
  formula for :math:`X`
* simplify the resulting expression,
* introduce an equation for :math:`X_{false}`, encoding :math:`X(false)`, which
  has the resulting expression as its right-hand side,
* recursively compute all equations for the predicate variables instances that
  occur in the resulting expression.
* as a final step: order every equation according to the ordering of the original PBES.

For the example, this yields the following strategy:

.. math::

    & (b \lor X(\neg b) \lor Y(b) )[b := false]\\
  = & X(true) \lor Y(false)

Introduce an equation :math:`(\mu X_{false} = X_{true} \lor Y_{false} )` and
continue with the computation for the equations for :math:`X(true)` and
:math:`Y(false)`. This yields two more equations:
:math:`(\mu X_{true} = true )` and
:math:`(\nu Y_{false} = X_{false} \land Y_{false})`. The resulting equations are
ordered with respect to the ordering of the original PBES, leading to the
following PBES:

.. math::

  \mu X_{false} & = X_{true} \lor Y_{false} \\ 
  \mu X_{true} & = true  \\
  \nu Y_{false} & = X_{false} \land Y_{false}

The resulting PBES is a BES, for which several well-documented algorithms exist
for computing the solution. The solution to :math:`X(false)` is effectively
encoded by the variable :math:`X_{false}`.


References
----------
.. [CPPW07a] T. Chen, B. Ploeger, J. van de Pol and T.A.C. Willemse. Equivalence
            Checking for Infinite Systems using Parameterized Boolean Equation
            Systems. In L. Caires, V.T. Vasconcelos (eds.), Concurrency Theory,
            18th International Conference, CONCUR 2007, Lisbon, Portugal,
            Lecture Notes in Computer Science 4703, Springer-Verlag, pp.
            120--135, 2007.
            `(DOI) <http://dx.doi.org/10.1007/978-3-540-74407-8_9>`_
