.. index:: ltscompare

.. _tool-ltscompare:

ltscompare
==========

The tool ``ltscompare`` compares two transition systems that both must have
either the extension .aut, .lts or .fsm. The transition systems should not
be probabilistic.

The transition systems can either be compared using an equivalence relation
(option ``--equivalence``) or a preorder (option ``--preorder``). The list of
available equivalences is provided below.

There are two useful options. One allows to generate counter examples (option
``--counter-example``). The counter example consist of a Hennessy-Milner formula
which is true in the first input LTS and false in the second. Counter-examples
are implemented and tested for the `bisim`, `branching-bisim` and `trace`
equivalence options. 

The second useful option is to hide some actions while doing the comparisons
(option ``--tau=`` followed by a comma separated list of actions). Counter examples
will only be distinguishing for the input transition systems with the hiding operation
applied.
