.. index:: lps2pbes

.. _tool-lps2pbes:

lps2pbes
========

The tool ``lps2pbes`` reads a modal formula as well as a linear process
and generates a parameterised boolean equation system (PBES) of which the solution
of the initial variable indicates whether the
formula is valid in the initial state of the transition system. The generated PBES can be solved
using tools such as :ref:`tool-pbessolve`, :ref:`tool-pbessolve` or :ref:`tool-pbespgsolve`.

When using particular formulas, for instance such as::

   [true*]<a.b.c.d.e.f.g.h>true

then the standard translation to PBESs can yield a very PBES which is very elaborate to generate and can
become large. This is due to the fact
that the subformula ``<a.b.c.d.e.f.g.h>`` is translated into one PBES equation with a huge right hand side.
This right hand side essentially reflects for any state of the linear process whether a trace ``a.b.c.d.e.f.g.h`` is possible.
When using the flag ``--preprocess-modal-operators`` the formula is first transformed into the equivalent
formula::

   [true*]mu X1.<a>mu X2.<b>mu X3.<c>mu X4.<d>mu X5.<d>mu X6.<e>mu X7.<f>mu X8.<g>mu X9.<h>true

This formula replaces the single very large equation by 9 ones, where the right hand sides only contain
the information whether a single action can be done. This is generally faster and yields a substantially
smaller PBES. Note that elaborate generations of PBESs can already occur when using subformulas of the
shape ``[a.b]...`` or ``<a.b>...``.

The tool :ref:`tool-pbessolve` is capable of generating a counter example in the form of labelled transition systems,
provided the PBES is generated
using lps2pbes (or :ref:`tool-lts2pbes`) using the ``--counter-example`` flag. The generated PBES is
more complicated
and may be harder to solve. Yet, these counter examples are very helpful in determining whether formulas do
not hold. If formulas are valid, this flag can also be used to determine witnesses, i.e., evidence when formulas
are valid. The tool :ref:`tool-pbessolve` can generate counter examples without the use of this flag. These
counter-examples are solely based on the provided PBES and they must be manually be related to the original
transition system.

Note that there is also an option ``--structured`` which can be used to generate boolean equation systems that
do not contain both the conjunction and the disjunction operators among PBES variables in the right hand side.
This flag can lead to a substantially larger number of equation (but linear in the original formula).

Note that the option ``--check-only`` can be used to check whether the passed
formula is valid w.r.t. parsing and type checking, without actually generating
the PBES.

.. mcrl2_manual:: lps2pbes
