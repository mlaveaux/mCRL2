.. index:: pbesstategraph

.. _tool-pbesstategraph:

pbesstategraph
==============

This tool performs a dead variable analysis on parameterised Boolean equation systems. Dead variables are reset to a constant.
The analysis performed in this tool is generally more powerful than the one performed by :ref:`tool-pbesparelm`, but it is also more resource intensive.

Example:

.. math::

   \begin{array}{l}
   \nu X(i,j,k,l{:}\mathbb{N}) = (i \neq 1 \lor j \neq 1 \lor X(2,j,k,l+1)) \land \forall m{:}\mathbb{N} . Z(i,2,m+k,k)\\
   \mu Y(i,j,k,l{:}\mathbb{N}) = k = 1 \lor (i = 2 \land X(1,j,k,l))\\
   \nu Z(i,j,k,l{:}\mathbb{N}) = (k < 10 \lor j = 2) \land (j \neq 2 \lor Y(1,1,l,1)) \land Y(2,2,1,l)\\
   ~\\
   \mathbf{init}\ X(1,1,1,1)
   \end{array}

Instantiation of this PBES using tools such as :ref:`tool-pbessolve` or :ref:`tool-pbessolve` does not terminate.
Also, :ref:`tool-pbesparelm` is only able to remove parameter :math:`l` from equations :math:`X` and :math:`Y`, and parameter :math:`i` from :math:`Z`. The presence of the predicate variable instantiations :math:`X(2,j,k,l+1)` and :math:`Z(i,2,m+k,k)` in :math:`X`'s equation means the solution to :math:`X(1,1,1,1)` depends on the solutions to :math:`X(2,1,1,2)` and :math:`X(1,2,v+1,1)` for all values :math:`v`, so instantiation will not terminate (not even after applying :ref:`tool-pbesparelm`, since parameter :math:`k` in the equation for :math:`Z` was not removed).

The dead variable analysis will (among others) observe that the third parameter of :math:`Z` does not affect the solution when :math:`j = 2`. The method will observe this, and replace :math:`Z(i,2,m+k,k)` by, for example, :math:`Z(i,2,1,k)`.

The resulting PBES will be:

.. math::

   \begin{array}{l}
   \nu X(i,j,k,l{:}\mathbb{N}) = (i \neq 1 \lor j \neq 1 \lor X(2,j,k,1)) \land \forall m{:}\mathbb{N} . Z(0,2,0,k)\\
   \mu Y(i,j,k,l{:}\mathbb{N}) = k = 1 \lor (i = 2 \land X(1,j,k,1))\\
   \nu Z(i,j,k,l{:}\mathbb{N}) = (k < 10 \lor j = 2) \land (j \neq 2 \lor Y(1,1,l,0)) \land Y(2,2,1,0)\\
   ~\\
   \mathbf{init}\ X(1,1,1,1)
   \end{array}

To perform the dead variable analysis, the algorithm computes a set of control flow variables. By default, local control flow graphs are constructed in order to avoid a combinatorial explosion of the control flow graph. The option `-g/--use-global-variant` uses a global control flow graph instead. The latter option may be able to detect more dead variables, but may use (much) more memory.

The algorithm underlying the tool, including the detection of control flow variables and the construction of local and global control flow graphs, is described in detail in [KWW14]_.

.. mcrl2_manual:: pbesstategraph
