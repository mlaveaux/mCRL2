Process Library Implementation Notes
======================================

*Author: Wieger Wesselink*

Processes
----------

Process expressions in mCRL2 are expressions built according to the following syntax:

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Expression
     - C++ equivalent
     - ATerm grammar
   * - :math:`a(e)`
     - ``action(a, e)``
     - ``Action``
   * - :math:`P(e)`
     - ``process(P, e)``
     - ``Process``
   * - :math:`P(d:=e)`
     - ``process_assignment(P, d:=e)``
     - ``ProcessAssignment``
   * - :math:`\delta`
     - ``delta()``
     - ``Delta``
   * - :math:`\tau`
     - ``tau()``
     - ``Tau``
   * - :math:`\sum_{d} x`
     - ``sum(d, x)``
     - ``Sum``
   * - :math:`\partial_B(x)`
     - ``block(B, x)``
     - ``Block``
   * - :math:`\tau_B(x)`
     - ``hide(B, x)``
     - ``Hide``
   * - :math:`\rho_R(x)`
     - ``rename(R, x)``
     - ``Rename``
   * - :math:`\Gamma_C(x)`
     - ``comm(C, x)``
     - ``Comm``
   * - :math:`\nabla_V(x)`
     - ``allow(V, x)``
     - ``Allow``
   * - :math:`x \mid y`
     - ``sync(x, y)``
     - ``Sync``
   * - :math:`x \mathbin{@} t`
     - ``at_time(x, t)``
     - ``AtTime``
   * - :math:`x \cdot y`
     - ``seq(x, y)``
     - ``Seq``
   * - :math:`c \rightarrow x`
     - ``if_then(c, x)``
     - ``IfThen``
   * - :math:`c \rightarrow x \diamond y`
     - ``if_then_else(c, x, y)``
     - ``IfThenElse``
   * - :math:`x \ll y`
     - ``binit(x, y)``
     - ``BInit``
   * - :math:`x \parallel y`
     - ``merge(x, y)``
     - ``Merge``
   * - :math:`x \lmerge y`
     - ``lmerge(x, y)``
     - ``LMerge``
   * - :math:`x + y`
     - ``choice(x, y)``
     - ``Choice``

where the types of the symbols are as follows:

.. list-table::
   :header-rows: 0
   :widths: auto

   * - :math:`a, b`
     - strings (action names)
   * - :math:`P`
     - a process identifier
   * - :math:`e`
     - a sequence of data expressions
   * - :math:`d`
     - a sequence of data variables
   * - :math:`B`
     - a set of strings (action names)
   * - :math:`R`
     - a sequence of rename expressions
   * - :math:`C`
     - a sequence of communication expressions
   * - :math:`V`
     - a sequence of multi actions
   * - :math:`t`
     - a data expression of type real
   * - :math:`x, y`
     - process expressions
   * - :math:`c`
     - a data expression of type bool

A rename expression is of the form :math:`a \rightarrow b`, with :math:`a` and :math:`b`
action names. A multi action is of the form :math:`a_1 \mid \cdots \mid a_n`, with
:math:`a_i` actions. A communication expression is of the form
:math:`b_1 \mid \cdots \mid b_n \rightarrow b`, with :math:`b` and :math:`b_i` action names.

Restrictions
~~~~~~~~~~~~~

A multi action is a multi set of actions. The left hand sides of the communication
expressions in :math:`C` must be unique. Also the left hand sides of the rename expressions
in :math:`R` must be unique.

Linear process expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Linear process expressions are a subset of process expressions satisfying the following
grammar:

.. code-block:: none

   <linear process expression> ::= choice(<linear process expression>, <linear process expression>)
                                 | <summand>

   <summand>                   ::= sum(<variables>, <alternative>)
                                 | <conditional action prefix>
                                 | <conditional deadlock>

   <conditional action prefix> ::= if_then(<condition>, <action prefix>)
                                 | <action prefix>

   <action prefix>             ::= seq(<timed multiaction>, <process reference>)
                                 | <timed multiaction>

   <timed multiaction>         ::= at_time(<multiaction>, <time stamp>)
                                 | <multiaction>

   <multiaction>               ::= tau()
                                 | <action>
                                 | sync(<multiaction>, <multiaction>)

   <conditional deadlock>      ::= if_then(<condition>, <timed deadlock>)
                                 | <timed deadlock>

   <timed deadlock>            ::= delta()
                                 | at_time(delta(), <time stamp>)

   <process reference>         ::= process(<process identifier>, <data expressions>)
                                 | process_assignment(<process identifier>, <data assignments>)

Guarded process expressions
-----------------------------

We define the predicate :math:`\mathit{is\_guarded}` for process expressions as follows:
:math:`\mathit{is\_guarded}(p) = \mathit{is\_guarded}(p, \emptyset)`

.. math::

   \begin{array}{lll}
   \mathit{is\_guarded}(a(e),\, W) & = & \mathit{true} \\
   \mathit{is\_guarded}(\delta,\, W) & = & \mathit{true} \\
   \mathit{is\_guarded}(\tau,\, W) & = & \mathit{true} \\
   \mathit{is\_guarded}(P(e),\, W) & = &
     \begin{array}{l}
       \left\{
       \begin{array}{ll}
         \mathit{false} & \text{if } P \in W \\
         \mathit{is\_guarded}(p,\, W \cup \{P\}) & \text{if } P \notin W
       \end{array}
       \right. \\
       \text{where } P(d) = p \text{ is the equation corresponding to } P(e)
     \end{array} \\[6pt]
   \mathit{is\_guarded}(p+q,\, W) & = & \mathit{is\_guarded}(p,\, W) \wedge \mathit{is\_guarded}(q,\, W) \\
   \mathit{is\_guarded}(p \cdot q,\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(c \rightarrow p,\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(c \rightarrow p \diamond q,\, W) & = & \mathit{is\_guarded}(p,\, W) \wedge \mathit{is\_guarded}(q,\, W) \\
   \mathit{is\_guarded}(\Sigma_{d:D}\, p,\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(p \mathbin{@} t,\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(p \ll q,\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(p \parallel q,\, W) & = & \mathit{is\_guarded}(p,\, W) \wedge \mathit{is\_guarded}(q,\, W) \\
   \mathit{is\_guarded}(p \lmerge q,\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(p \mid q,\, W) & = & \mathit{is\_guarded}(p,\, W) \wedge \mathit{is\_guarded}(q,\, W) \\
   \mathit{is\_guarded}(\rho_R(p),\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(\partial_B(p),\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(\tau_I(p),\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(\Gamma_C(p),\, W) & = & \mathit{is\_guarded}(p,\, W) \\
   \mathit{is\_guarded}(\nabla_V(p),\, W) & = & \mathit{is\_guarded}(p,\, W)
   \end{array}

N.B. This specification assumes that process names are unique. In mCRL2 process names can be
overloaded, therefore in the implementation :math:`W` contains *process identifiers* (i.e. both
the process name and the sorts of the arguments) instead of process names.

Alphabet reduction
-------------------

Alphabet reduction is a preprocessing step for linearization. It is a transformation on process
expressions that preserves branching bisimulation. The transformation distributes renaming
operators (``allow``, ``hide``, ``rename``, ``block``, and ``comm``) over the other process
operators, with the goal of reducing the number of multi-actions that arise during intermediate
steps of the linearization procedure.

Notations
~~~~~~~~~~

In this text action names are represented using :math:`a, b, \ldots` and multi action names
using :math:`\alpha, \beta, \ldots`. So in general we have
:math:`\alpha = a_1 \mid \ldots \mid a_n`. In alphabet reduction data parameters play a minor
role, therefore we choose a notation in which data parameters are omitted. We use the
abbreviation :math:`\overline{a} = a(e_1, \ldots, e_n)` to denote an action, and
:math:`\overline{\alpha} = \overline{a_1} \mid \ldots \mid \overline{a_n}` to denote a multi
action, where :math:`e_1, \ldots, e_n` are data expressions. Note that a multi action is a
multiset (or bag) of actions and a multi action name is a multiset of names. We write
:math:`\alpha\beta` as shorthand for :math:`\alpha \cup \beta` and :math:`a\beta` for
:math:`\{a\} \cup \beta`. Sets of multi action names are represented using
:math:`A, A_1, A_2, \ldots`. A communication :math:`C` maps multi action names to action
names, and is denoted as :math:`\{\alpha_1 \rightarrow a_1, \ldots, \alpha_n \rightarrow a_n\}`.
A renaming :math:`R` is a substitution on action names, and is denoted as
:math:`R = \{a_1 \rightarrow b_1, \ldots, a_n \rightarrow b_n\}`. A block set :math:`B` is a
set of action names. A hide set :math:`I` is a set of action names.

Definitions
~~~~~~~~~~~~

We define multi actions :math:`\overline{\alpha}` using the following grammar:

.. math::

   \overline{\alpha} := \overline{a} \shortmid \overline{\alpha} \mid \overline{a},

where :math:`\overline{a}` is an action, and where :math:`\shortmid` is used to distinguish
alternatives.

We define pCRL terms :math:`p` using the following grammar:

.. math::

   p ::= \overline{a} \shortmid P \shortmid \delta \shortmid \tau \shortmid p+p \shortmid
         p \cdot p \shortmid c \rightarrow p \shortmid c \rightarrow p \diamond p \shortmid
         \Sigma_{d:D} p \shortmid p \mathbin{@} t \shortmid p \ll p,

and parallel mCRL terms :math:`q` using the following grammar:

.. math::

   q ::= p \shortmid q \parallel q \shortmid q \lmerge q \shortmid q \mid q \shortmid
         \rho_R(q) \shortmid \partial_B(q) \shortmid \tau_I(q) \shortmid
         \Gamma_C(q) \shortmid \nabla_V(q).

**Remark 1**

Note that there is an unfortunate overload of the :math:`|`-operator in both multi actions and
process expressions. This has consequences for the implementation, since there is no clean
distinction between parallel and non-parallel operators.

**Remark 2**

The mCRL2 language also has a construct
:math:`P(d_{i_1} = e_{i_1}, \ldots, d_{i_k} = e_{i_k})`, but this is just a shorthand
notation. Therefore we will ignore it in this text.

Alphabet operations
~~~~~~~~~~~~~~~~~~~~

Let :math:`A`, :math:`A_1` and :math:`A_2` be sets of multi action names. Then we define

.. math::

   \begin{array}{lll}
   A^{\subseteq} & = & \{\alpha \mid \exists \beta.\, \alpha\beta \in A\} \\
   A_1 A_2 & = & \{\alpha\beta \mid \alpha \in A_1 \text{ and } \beta \in A_2\} \\
   A_1 \leftarrowtail A_2 & = & \{\alpha \mid \exists \beta.\, \alpha\beta \in A_1 \text{ and } \beta \in A_2\}
   \end{array}

Note that :math:`\beta` can take the value :math:`\tau` in the definition of
:math:`A_1 \leftarrowtail A_2`, which implies :math:`A_1 \subset A_1 \leftarrowtail A_2`.
The set :math:`A^{\subseteq}` has an exponential size, so whenever possible it should not be
computed explicitly.

Let :math:`C` be a communication set, then we define

.. math::

   \begin{array}{lll}
   C(A) & = & \bigcup_{\alpha \in A} \operatorname{Comm}(C,\, \alpha) \\
   C^{-1}(A) & = & \bigcup_{\alpha \in A} \operatorname{CommInverse}(C,\, \alpha) \\
   \mathit{filter}_{\nabla}(C,\, A) & = &
     \{\gamma \rightarrow c \in C \mid \exists_{\alpha \in A}.\, \gamma \subset \alpha\}
   \end{array}

where :math:`\operatorname{Comm}` and :math:`\operatorname{CommInverse}` are defined using
pseudo code as follows:

.. math::

   \begin{array}{l}
   \operatorname{Comm}(C,\, \alpha) \\
   R := \{\alpha\} \\
   \mathbf{for}\ \gamma \rightarrow c \in C\ \mathbf{do} \\
   \qquad \mathbf{if}\ \exists \beta.\, \alpha = \beta\gamma
     \ \mathbf{then}\ R := R \cup \operatorname{Comm}(C,\, \beta c) \\
   \mathbf{return}\ R
   \end{array}

.. math::

   \begin{array}{l}
   \operatorname{CommInverse}(C,\, \alpha_1,\, \alpha_2) \\
   R := \{\alpha_1 \alpha_2\} \\
   \mathbf{for}\ \gamma \rightarrow c \in C\ \mathbf{do} \\
   \qquad \mathbf{if}\ \exists \beta.\, \alpha_1 = \beta c
     \ \mathbf{then}\ R := R \cup \operatorname{CommInverse}(C,\, \beta,\, \alpha_2\gamma) \\
   \mathbf{return}\ R
   \end{array}

Note that :math:`C^{-1}(\alpha) = \operatorname{CommInverse}(C,\, \alpha,\, \tau)`.

Let :math:`R` be a rename set, then we define

.. math::

   \begin{array}{lll}
   R(\alpha)      & = & \{R(\alpha_i) \mid \alpha_i \in \alpha\} \\
   R^{-1}(\alpha) & = & \{\beta \mid R(\beta) = \alpha\} \\
   R(A)           & = & \{R(\alpha) \mid \alpha \in A\} \\
   R^{-1}(A)      & = & \{R^{-1}(\alpha) \mid \alpha \in A\}
   \end{array}

Let :math:`I` be a hide set, then we define

.. math::

   \begin{array}{lll}
   \tau_I(A)      & = & \{\beta \mid \exists_{\alpha \in A,\, \gamma \in I^*}.\,
                          \alpha = \beta\gamma \wedge \beta \cap I = \emptyset\} \\
   \tau_I^{-1}(A) & = & \partial_I(A)\, I^*
   \end{array}

Let :math:`B` be a block set, then we define

.. math::

   \partial_B(A) = \{\alpha \in A \mid \alpha \cap B = \emptyset\}

We define a mapping :math:`act` that extracts the individual action names of a set of multi
action names:

.. math::

   \begin{array}{lll}
   act(a_1 \mid \cdots \mid a_n) & = & \{a_1 \mid \cdots \mid a_n\} \\
   act(A)                        & = & \bigcup_{\alpha \in A} act(\alpha)
   \end{array}

The mapping :math:`\alpha`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We define the mapping :math:`\alpha` as follows. The value :math:`\alpha(p, \emptyset)` is an
over-approximation of the alphabet of process expression :math:`p`.

.. math::

   \begin{array}{lll}
   \alpha(\overline{a},\, W) & = & \{a\} \\[4pt]
   \alpha(P,\, W) & = &
     \begin{array}{l}
       \left\{
       \begin{array}{ll}
         \emptyset & \text{if } P \in W \\
         \alpha(p,\, W \cup \{P\}) & \text{if } P \notin W,
       \end{array}
       \right. \\
       \text{where } P = p \text{ is the equation of } P
     \end{array} \\[10pt]
   \alpha(\delta,\, W) & = & \emptyset \\
   \alpha(\tau,\, W)   & = & \{\tau\} \\
   \alpha(p+q,\, W)    & = & \alpha(p,\, W) \cup \alpha(q,\, W) \\
   \alpha(p \cdot q,\, W) & = & \alpha(p,\, W) \cup \alpha(q,\, W) \\
   \alpha(c \rightarrow p,\, W) & = & \alpha(p,\, W) \\
   \alpha(c \rightarrow p \diamond q,\, W) & = & \alpha(p,\, W) \cup \alpha(q,\, W) \\
   \alpha(\Sigma_{d:D}\, p,\, W) & = & \alpha(p,\, W) \\
   \alpha(p \mathbin{@} t,\, W) & = & \alpha(p,\, W) \\
   \alpha(p \ll q,\, W) & = & \alpha(p,\, W) \cup \alpha(q,\, W) \\
   \alpha(p \parallel q,\, W) & = &
     \alpha(p,\, W) \cup \alpha(q,\, W) \cup \alpha(p,\, W)\alpha(q,\, W) \\
   \alpha(p \lmerge q,\, W) & = &
     \alpha(p,\, W) \cup \alpha(q,\, W) \cup \alpha(p,\, W)\alpha(q,\, W) \\
   \alpha(p \mid q,\, W) & = & \alpha(p,\, W)\alpha(q,\, W) \\
   \alpha(\rho_R(p),\, W) & = & R(\alpha(p,\, W)) \\
   \alpha(\partial_B(p),\, W) & = & \partial_B(\alpha(p,\, W)) \\
   \alpha(\tau_I(p),\, W) & = & \tau_I(\alpha(p,\, W)) \\
   \alpha(\Gamma_C(p),\, W) & = & C(\alpha(p,\, W)) \\
   \alpha(\nabla_V(p),\, W) & = & \alpha(p,\, W) \cap (V \cup \{\tau\})
   \end{array}

*Example 1*

If :math:`C = \{a \mid b \rightarrow c\}`, then
:math:`\alpha(\Gamma_C(a(1) \mid b(2))) = \{a, b, c, a \mid b\}`. Note that the action
:math:`c` does not occur in the transition system of this process expression.

*Example 2*

In the computation of
:math:`\{a_1, a_2, \ldots, a_{20}\} \cap \alpha(a_1 \parallel a_2 \parallel \cdots \parallel a_{20})`
the above mentioned optimization is really needed.

Computation of the alphabet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When computing :math:`A \cap \alpha(p, W)` for some multi action name set :math:`A`, it may be
beneficial to apply an optimization. This is done to keep intermediate expressions small. We
introduce :math:`\alpha(p, W, A) = A \cap \alpha(p, W)`, and define it as follows:

.. math::

   \begin{array}{lll}
   \alpha(\overline{a},\, W,\, A) & = &
     \left\{
     \begin{array}{ll}
       \{a\} & \text{if } a \in A \\
       \emptyset & \text{if } a \notin A
     \end{array}
     \right. \\[6pt]
   \alpha(P,\, W,\, A) & = &
     \begin{array}{l}
       \left\{
       \begin{array}{ll}
         \emptyset & \text{if } P \in W \\
         \alpha(p,\, W \cup \{P\},\, A) & \text{if } P \notin W,
       \end{array}
       \right. \\
       \text{where } P = p \text{ is the equation of } P
     \end{array} \\[10pt]
   \alpha(p+q,\, W,\, A) & = & \alpha(p,\, W,\, A) \cup \alpha(q,\, W,\, A) \\
   \alpha(p \cdot q,\, W,\, A) & = & \alpha(p,\, W,\, A) \cup \alpha(q,\, W,\, A) \\
   \alpha(c \rightarrow p,\, W,\, A) & = & \alpha(p,\, W,\, A) \\
   \alpha(c \rightarrow p \diamond q,\, W,\, A) & = &
     \alpha(p,\, W,\, A) \cup \alpha(q,\, W,\, A) \\
   \alpha(\Sigma_{d:D}\, p,\, W,\, A) & = & \alpha(p,\, W,\, A) \\
   \alpha(p \mathbin{@} t,\, W,\, A) & = & \alpha(p,\, W,\, A) \\
   \alpha(p \ll q,\, W,\, A) & = & \alpha(p,\, W,\, A) \cup \alpha(q,\, W,\, A) \\
   \alpha(p \parallel q,\, W,\, A) & = &
     \alpha(p,\, W,\, A) \cup \alpha(q,\, W,\, A) \cup
     \alpha(p,\, W,\, A^{\subseteq})\alpha(q,\, W,\, A^{\subseteq}) \\
   \alpha(p \lmerge q,\, W,\, A) & = &
     \alpha(p,\, W,\, A) \cup \alpha(q,\, W,\, A) \cup
     \alpha(p,\, W,\, A^{\subseteq})\alpha(q,\, W,\, A^{\subseteq}) \\
   \alpha(p \mid q,\, W,\, A) & = &
     \alpha(p,\, W,\, A^{\subseteq})\alpha(q,\, W,\, A^{\subseteq})
   \end{array}

More efficient computation of the alphabet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The computation of :math:`\alpha(p, W, A)` can be done more efficiently. We define the
function :math:`\mathit{proc}(p, W)` as follows:

.. math::

   \begin{array}{lll}
   \mathit{proc}(\overline{a},\, W) & = & \emptyset \\[4pt]
   \mathit{proc}(P,\, W) & = &
     \left\{
     \begin{array}{ll}
       \emptyset & \text{if } P \in W \\
       \{P\} \cup \mathit{proc}(p,\, W) & \text{if } P \notin W
     \end{array}
     \right. \\[6pt]
   \mathit{proc}(p+q,\, W) & = & \mathit{proc}(p,\, W) \cup \mathit{proc}(q,\, W) \\
   \mathit{proc}(p \cdot q,\, W) & = & \mathit{proc}(p,\, W) \cup \mathit{proc}(q,\, W) \\
   \mathit{proc}(c \rightarrow p,\, W) & = & \mathit{proc}(p,\, W) \\
   \mathit{proc}(c \rightarrow p \diamond q,\, W) & = &
     \mathit{proc}(p,\, W) \cup \mathit{proc}(q,\, W) \\
   \mathit{proc}(\Sigma_{d:D}\, p,\, W) & = & \mathit{proc}(p,\, W) \\
   \mathit{proc}(p \mathbin{@} t,\, W) & = & \mathit{proc}(p,\, W)
   \end{array}

Using this function we can change the computation of :math:`\alpha(p, W, A)` at three places:

.. math::

   \begin{array}{lll}
   \alpha(p+q,\, W,\, A) & = &
     \alpha(p,\, W,\, A) \cup \alpha(q,\, W \cup \mathit{proc}(p, W),\, A) \\
   \alpha(p \cdot q,\, W,\, A) & = &
     \alpha(p,\, W,\, A) \cup \alpha(q,\, W \cup \mathit{proc}(p, W),\, A) \\
   \alpha(c \rightarrow p \diamond q,\, W,\, A) & = &
     \alpha(p,\, W,\, A) \cup \alpha(q,\, W \cup \mathit{proc}(p, W),\, A)
   \end{array}

Note that the value :math:`\mathit{proc}(p, W)` can be computed on the fly during the
computation of :math:`\alpha(p, W, A)`.

Bounded alphabet
~~~~~~~~~~~~~~~~~

In practice one often wants to compute :math:`\alpha(p, A) = \alpha(\nabla_A(p))`. This can
be computed more efficiently as follows:

.. math::

   \begin{array}{lll}
   \alpha(\overline{a},\, A) & = &
     \left\{
     \begin{array}{ll}
       \{a\} & \text{if } a \in A \\
       \emptyset & \text{if } a \notin A
     \end{array}
     \right. \\[6pt]
   \alpha(P,\, A) & = & \alpha(p,\, A),\quad \text{where } P = p \text{ is the equation of } P \\
   \alpha(p+q,\, A) & = & \alpha(p,\, A) \cup \alpha(q,\, A) \\
   \alpha(p \cdot q,\, A) & = & \alpha(p,\, A) \cup \alpha(q,\, A) \\
   \alpha(c \rightarrow p,\, A) & = & \alpha(p,\, A) \\
   \alpha(c \rightarrow p \diamond q,\, A) & = & \alpha(p,\, A) \cup \alpha(q,\, A) \\
   \alpha(\Sigma_{d:D}\, p,\, A) & = & \alpha(p,\, A) \\
   \alpha(p \mathbin{@} t,\, A) & = & \alpha(p,\, A) \\
   \alpha(p \ll q,\, A) & = & \alpha(p,\, A) \cup \alpha(q,\, A) \\
   \alpha(p \parallel q,\, A) & = &
     \alpha(p,\, A) \cup \alpha(q,\, A) \cup
     \alpha(p,\, A^{\subseteq})\alpha(q,\, A \leftarrowtail \alpha(p,\, A^{\subseteq})) \\
   \alpha(p \lmerge q,\, A) & = &
     \alpha(p,\, A) \cup \alpha(q,\, A) \cup
     \alpha(p,\, A^{\subseteq})\alpha(q,\, A \leftarrowtail \alpha(p,\, A^{\subseteq})) \\
   \alpha(p \mid q,\, A) & = &
     \alpha(p,\, A^{\subseteq})\alpha(q,\, A \leftarrowtail \alpha(p,\, A^{\subseteq})) \\
   \alpha(\rho_R(p),\, A) & = & R(\alpha(p,\, R^{-1}(A))) \\
   \alpha(\partial_B(p),\, A) & = & \alpha(p,\, \partial_B(A)) \\
   \alpha(\tau_I(p),\, A) & = & \tau_I(\alpha(p,\, \tau_I^{-1}(A))) \\
   \alpha(\Gamma_C(p),\, A) & = & C(\alpha(p,\, C^{-1}(A))) \\
   \alpha(\nabla_V(p),\, A) & = & \alpha(p,\, A \cap V)
   \end{array}

The mappings :math:`push`, :math:`push_{\nabla}` and :math:`push_{\partial}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We define mappings :math:`push`, :math:`push_{\nabla}` and :math:`push_{\partial}` such that
:math:`push(p)` is bisimulation equivalent to :math:`p`,
:math:`push_{\nabla}(A, p)` is bisimulation equivalent to :math:`\nabla_A(p)`, and
:math:`push_{\partial}(B, p)` is bisimulation equivalent to :math:`\partial_B(p)`. The goal
of these mappings is to push allow and block expressions deeply inside process expressions.
It is important to know that an allow set :math:`A` in the expression :math:`\nabla_A(p)`
implicitly contains the empty multi action :math:`\tau`. Let
:math:`\mathcal{E} = \{P_1(d) = p_1, \ldots, P_n(d) = p_n\}` be a sequence of process
equations.

.. math::

   \begin{array}{lll}
   push(p)              & = & p \quad \text{if } p \text{ is a pCRL expression} \\
   push(p \parallel q)  & = & push(p) \parallel push(q) \\
   push(p \lmerge q)    & = & push(p) \lmerge push(q) \\
   push(p \mid q)       & = & push(p) \mid push(q) \\
   push(\rho_R(p))      & = & \rho_R(push(p)) \\
   push(\partial_B(p))  & = & push_{\partial}(B,\, p) \\
   push(\tau_I(p))      & = & \tau_I(push(p)) \\
   push(\Gamma_C(p))    & = & \Gamma_C(push(p)) \\
   push(\nabla_V(p))    & = & push_{\nabla}(V,\, p)
   \end{array}

We assume that :math:`P_A^{\nabla}` is a unique name for every
:math:`P \in \{P_1, \ldots, P_n\}`, multi action name set :math:`A` and sequence of data
expressions :math:`e`.

.. math::

   \begin{array}{lll}
   push_{\nabla}(A,\, \overline{a}) & = &
     \left\{
     \begin{array}{ll}
       \overline{a} & \text{if } N(\overline{a}) \in A \\
       \delta & \text{otherwise}
     \end{array}
     \right. \\[6pt]
   push_{\nabla}(A,\, P(e)) & = &
     \begin{array}{l}
       P_A^{\nabla}(e),\quad \text{where } P(d) = p \text{ is the equation of } P, \\
       \text{and where } P_A^{\nabla}(d) = push_{\nabla}(A,\, p) \text{ is a new equation}
     \end{array} \\[8pt]
   push_{\nabla}(A,\, \delta) & = & \delta \\
   push_{\nabla}(A,\, \tau)   & = & \tau \\
   push_{\nabla}(A,\, p+q)    & = & \nabla_A(p+q) \\
   push_{\nabla}(A,\, p \cdot q) & = & \nabla_A(p \cdot q) \\
   push_{\nabla}(A,\, c \rightarrow p) & = & \nabla_A(c \rightarrow p) \\
   push_{\nabla}(A,\, c \rightarrow p \diamond q) & = & \nabla_A(c \rightarrow p \diamond q) \\
   push_{\nabla}(A,\, \Sigma_{d:D} p) & = & \nabla_A(\Sigma_{d:D} p) \\
   push_{\nabla}(A,\, p \mathbin{@} t) & = & \nabla_A(p \mathbin{@} t) \\
   push_{\nabla}(A,\, p \ll q) & = & \nabla_A(p \ll q) \\
   push_{\nabla}(A,\, p \parallel q) & = &
     \nabla_A(p' \parallel q')\ \text{where }
     \left\{
     \begin{array}{lll}
       p' & = & push_{\nabla}(A^{\subseteq},\, p) \\
       q' & = & push_{\nabla}(A \leftarrowtail \alpha(p'),\, q)
     \end{array}
     \right. \\[6pt]
   push_{\nabla}(A,\, p \lmerge q) & = &
     \nabla_A(p' \lmerge q')\ \text{where }
     \left\{
     \begin{array}{lll}
       p' & = & push_{\nabla}(A^{\subseteq},\, p) \\
       q' & = & push_{\nabla}(A \leftarrowtail \alpha(p'),\, q)
     \end{array}
     \right. \\[6pt]
   push_{\nabla}(A,\, p \mid q) & = &
     \nabla_A(p' \mid q')\ \text{where }
     \left\{
     \begin{array}{lll}
       p' & = & push_{\nabla}(A^{\subseteq},\, p) \\
       q' & = & push_{\nabla}(A \leftarrowtail \alpha(p'),\, q)
     \end{array}
     \right. \\[6pt]
   push_{\nabla}(A,\, \rho_R(p)) & = &
     \rho_R(p')\ \text{where } p' = push_{\nabla}(R^{-1}(A),\, p) \\
   push_{\nabla}(A,\, \partial_B(p)) & = & push_{\nabla}(\partial_B(A),\, p) \\
   push_{\nabla}(A,\, \tau_I(p)) & = &
     \tau_I(p')\ \text{where } p' = push_{\nabla}(\tau_I^{-1}(A),\, p) \\
   push_{\nabla}(A,\, \Gamma_C(p)) & = &
     \mathsf{allow}(A,\, \Gamma_C(p'))\ \text{where } p' = push_{\nabla}(C^{-1}(A),\, p) \\
   push_{\nabla}(A,\, \nabla_V(p)) & = & push_{\nabla}(A \cap V,\, p)
   \end{array}

**Optimizations**

During the computation of :math:`push_{\nabla}` the following optimizations are applied in
the right hand side of each equation:

.. math::

   \begin{array}{lll}
   \nabla_A(p) & = &
     \left\{
     \begin{array}{ll}
       p & \text{if } (A \cup \{\tau\}) \cap \alpha(p) = \alpha(p) \\
       \nabla_{A \cap \alpha(p)}(p) & \text{otherwise}
     \end{array}
     \right. \\[6pt]
   \nabla_{\emptyset}(p) & = &
     \left\{
     \begin{array}{ll}
       \tau   & \text{if } p = \tau \\
       \delta & \text{otherwise}
     \end{array}
     \right. \\[6pt]
   \Gamma_C(p) & = & \Gamma_{\mathit{filter}_{\nabla}(C,\, \alpha(p))}(p) \\
   \delta \mid \delta       & = & \delta \\
   \delta \parallel \delta  & = & \delta
   \end{array}

For non-pCRL expressions the alphabet :math:`\alpha(p)` is computed on the fly during the
computation of :math:`push_{\nabla}(A, p)`.

**Example 1**

Let :math:`P = (a+b) \cdot P`. Then :math:`push_{\nabla}(\{a\}, P, \emptyset) = P'`, with
:math:`P' = push_{\nabla}(\{a\}, (a+b) \cdot P, \{(P, \{a\}, P')\}) =
push_{\nabla}(\{a\}, (a+b), \{(P, \{a\}, P')\}) \cdot
push_{\nabla}(\{a\}, P, \{(P, \{a\}, P')\}) = \cdots = a \cdot P'`.

**Example 2**

Let :math:`P = a \cdot \nabla_{\{a\}}(P)`. Then
:math:`push_{\nabla}(\{a\}, P, \emptyset) = P'`, with
:math:`P' = push_{\nabla}(\{a\}, a \cdot \nabla_{\{a\}}(P), \{(P, \{a\}, P')\}) =
push_{\nabla}(\{a\}, a, \{(P, \{a\}, P')\}) \cdot
push_{\nabla}(\{a\}, \nabla_{\{a\}}(P), \{(P, \{a\}, P')\}) = \cdots = a \cdot P'`.

We assume that :math:`P_{B,e}^{\partial}` is a unique name for every
:math:`P \in \{P_1, \ldots, P_n\}`, block set :math:`B` and sequence of data expressions
:math:`e`.

.. math::

   \begin{array}{lll}
   push_{\partial}(B,\, \overline{a}) & = &
     \left\{
     \begin{array}{ll}
       \overline{a} & \text{if } N(\overline{a}) \cap B = \emptyset \\
       \delta & \text{otherwise}
     \end{array}
     \right. \\[6pt]
   push_{\partial}(B,\, P(e)) & = &
     \begin{array}{l}
       P_{B,e}^{\partial}(e), \\
       \text{where } P(d) = p \text{ is the equation of } P, \\
       \text{and where } P_{B,e}^{\partial}(d) = push_{\partial}(B,\, p) \text{ is a new equation}
     \end{array} \\[8pt]
   push_{\partial}(B,\, \delta) & = & \delta \\
   push_{\partial}(B,\, \tau)   & = & \tau \\
   push_{\partial}(B,\, p+q)    & = &
     push_{\partial}(B,\, p) + push_{\partial}(B,\, q) \\
   push_{\partial}(B,\, p \cdot q) & = &
     push_{\partial}(B,\, p) \cdot push_{\partial}(B,\, q) \\
   push_{\partial}(B,\, c \rightarrow p) & = &
     c \rightarrow push_{\partial}(B,\, p) \\
   push_{\partial}(B,\, c \rightarrow p \diamond q) & = &
     c \rightarrow push_{\partial}(B,\, p) \diamond push_{\partial}(B,\, q) \\
   push_{\partial}(B,\, \Sigma_{d:D} p) & = & \Sigma_{d:D}\, push_{\partial}(B,\, p) \\
   push_{\partial}(B,\, p \mathbin{@} t) & = & push_{\partial}(B,\, p) \mathbin{@} t \\
   push_{\partial}(B,\, p \ll q) & = &
     push_{\partial}(B,\, p) \ll push_{\partial}(B,\, q) \\
   push_{\partial}(B,\, p \parallel q) & = &
     push_{\partial}(B,\, p) \parallel push_{\partial}(B,\, q) \\
   push_{\partial}(B,\, p \lmerge q) & = &
     push_{\partial}(B,\, p) \lmerge push_{\partial}(B,\, q) \\
   push_{\partial}(B,\, p \mid q) & = &
     push_{\partial}(B,\, p) \mid push_{\partial}(B,\, q) \\
   push_{\partial}(B,\, \rho_R(p)) & = &
     \rho_R\!\left(push_{\partial}(R^{-1}(B),\, p)\right) \\
   push_{\partial}(B,\, \partial_{B_1}(p)) & = &
     push_{\partial}(B \cup B_1,\, p) \\
   push_{\partial}(B,\, \tau_I(p)) & = &
     \tau_I\!\left(push_{\partial}(B \setminus I,\, p)\right) \\
   push_{\partial}(B,\, \Gamma_C(p)) & = &
     \mathsf{block}(B,\, \Gamma_C(push_{\partial}(B',\, p))) \\
     & & \quad \text{where } B' = B \setminus
       \{b \in B \mid \exists_{\gamma \rightarrow c \in C}.\, b \in \gamma \wedge c \notin B\} \\
   push_{\partial}(B,\, \nabla_V(p)) & = &
     push_{\nabla}(\partial_B(V),\, p,\, \emptyset)
   \end{array}

where

.. math::

   \mathsf{block}(B,\, p) =
     \left\{
     \begin{array}{ll}
       p              & \text{if } B = \emptyset \\
       \partial_B(p)  & \text{otherwise}
     \end{array}
     \right.

**Example 3**

The presence of :math:`R^{-1}(\partial_B(A))` instead of just :math:`R^{-1}(A)` in the
right hand side of the rename operator is explained by the example
:math:`push_{\nabla}(\{b\}, \rho_{\{b \rightarrow c\}} b)`. We see that
:math:`\rho_{\{b \rightarrow c\}} push_{\nabla}(R^{-1}(A), p) =
\rho_{\{b \rightarrow c\}} push_{\nabla}(\{b\}, b) = \rho_{\{b \rightarrow c\}} b = c`,
which is clearly the wrong answer.

Allow sets
~~~~~~~~~~~

There are two rules in the definition of :math:`push_{\nabla}` where the allow set
can/should not be computed explicitly. The computation of
:math:`push_{\nabla}(A, p \parallel q)` involves computation of
:math:`push_{\nabla}(p, A^{\subseteq})`. We want to avoid the computation of
:math:`A^{\subseteq}`, since it can become very large. The computation of
:math:`push_{\nabla}(A, \tau_I(p))` involves computation of
:math:`push_{\nabla}(p, \tau_I^{-1}(A))`. The set :math:`\tau_I^{-1}(A) = AI^*` is
infinite.

In the implementation we use allow sets of the form :math:`A^{\subseteq} I^*`, where
:math:`A` is a set of multi action names and :math:`I` is a set of action names. The
:math:`{}^{\subseteq}` is optional and :math:`I` may be empty. Such an allow set is stored
as two sets :math:`A` and :math:`I`, together with an attribute that tells if
:math:`{}^{\subseteq}` is applicable. We need to show that allow sets are closed under the
operations in :math:`push_{\nabla}`.

.. math::

   \begin{array}{lll}
   \partial_B(A^{\subseteq} I^*) & = & \tau_B(A)^{\subseteq}\, \tau_B(I)^* \\
   \tau_{I_1}^{-1}(A^{\subseteq} I^*) & = & \partial_{I_1}(A^{\subseteq})(I \cup I_1)^* \\
   (A^{\subseteq} I^*) \cap V & = &
     \{\beta \in V \mid \exists_{\alpha \in A}.\, \tau_I(\beta) \sqsubseteq \alpha\} \\
   R^{-1}(A^{\subseteq} I^*) & = & R^{-1}(A^{\subseteq})\, R^{-1}(I)^* \\
   C^{-1}(A^{\subseteq} I^*) & \subseteq &
     C^{-1}(A)^{\subseteq}\, act(C^{-1}(I))^* \\
   (A^{\subseteq} I^*) \leftarrowtail A_1 & = & A^{\subseteq} I^* \\
   (A^{\subseteq} I^*)^{\subseteq} & = & A^{\subseteq} I^* \\[4pt]
   \partial_B(A I^*) & = & \partial_B(A)\, \tau_B(I)^* \\
   \tau_{I_1}^{-1}(A I^*) & = & \partial_{I_1}(A)(I \cup I_1)^* \\
   (A I^*) \cap V & = &
     \{\beta \in V \mid \exists_{\alpha \in A}.\, \tau_I(\beta) = \alpha\} \\
   R^{-1}(A I^*) & = & R^{-1}(A)\, R^{-1}(I)^* \\
   C^{-1}(A I^*) & \subseteq & C^{-1}(A)\, act(C^{-1}(I))^* \\
   (A I^*)^{\subseteq} & = & A^{\subseteq} I^*
   \end{array}

where we used the following properties:

.. math::

   \begin{array}{lll}
   \partial_B(A_1 A_2) & = & \partial_B(A_1)\, \partial_B(A_2) \\
   \partial_B(A^{\subseteq}) & = & \tau_B(A)^{\subseteq} \\
   R^{-1}(A_1 A_2) & = & R^{-1}(A_1)\, R^{-1}(A_2) \\
   R^{-1}(A^*) & = & R^{-1}(A)^* \\
   C^{-1}(A^{\subseteq}) & \subseteq & C^{-1}(A)^{\subseteq} \\
   C^{-1}(A_1 A_2) & = & C^{-1}(A_1)\, C^{-1}(A_2) \\
   C^{-1}(A^*) & = & C^{-1}(A)^* \\
   A^{\subseteq} \leftarrowtail A_1 & = & A^{\subseteq}
   \end{array}

Note that in case of the communication we only have an inclusion relation instead of equality.
This is done to stay within the format :math:`A^{\subseteq} I^*`. As a consequence the
implementation uses an over-approximation of :math:`C^{-1}(A^{\subseteq} I^*)` and
:math:`C^{-1}(A I^*)`. Furthermore note that the property
:math:`R^{-1}(A^{\subseteq}) = R^{-1}(A)^{\subseteq}` does not hold. A counter example is
:math:`R = \{b \rightarrow a\}` and :math:`A = \{a, b \mid c\}`. In that case we have
:math:`R^{-1}(A^{\subseteq}) = \{a, b, c\}^{\subseteq}` and
:math:`R^{-1}(A)^{\subseteq} = \{a, b\}^{\subseteq}`. Another property that was initially
assumed, but that does not hold, is
:math:`(A I^*) \leftarrowtail A_1 = (A \leftarrowtail \tau_I(A_1))\, I^*`.

Optimization for :math:`push_{\nabla}`
----------------------------------------

In some cases the :math:`push_{\nabla}` operator produces expressions that are too large.
This section proposes an optimization for the case :math:`push_{\nabla}(A, \Gamma_C(p))`
that can help to prevent this problem for certain practical cases.

.. math::

   push_{\nabla}(A,\, \Gamma_C(p)) =
     \left\{
     \begin{array}{ll}
       \mathsf{allow}(A,\, \Gamma_{C \setminus C'}(push_{\nabla\Gamma}(A',\, C',\, p)))
         & \text{if } C \neq C' \\
       push_{\nabla\Gamma}(A,\, C,\, p) & \text{otherwise,}
     \end{array}
     \right.

with
:math:`C' = \{\beta \rightarrow b \in C \mid b \notin \bigcup_{\beta' \rightarrow b' \in C} \beta'\}`
and :math:`A' = ((C \setminus C')(A))^{\subseteq}` and

.. math::

   \begin{array}{lll}
   push_{\nabla\Gamma}(A,\, C,\, p \parallel q) & = &
     \mathsf{allow}(A,\, \Gamma_C(\mathsf{allow}(C^{-1}(A),\, p' \parallel q')))
     \ \text{where}
     \left\{
     \begin{array}{lll}
       p'  & = & push_{\nabla\Gamma}(A',\, C,\, p) \\
       q'  & = & push_{\nabla\Gamma}(A'',\, C,\, q) \\
       A'  & = & C^{-1}(A)^{\subseteq} \setminus (C^{-1}(A) \setminus A) \\
       A'' & = & (C^{-1}(A) \leftarrowtail \alpha(p')) \setminus (C^{-1}(A) \setminus A)
     \end{array}
     \right. \\[8pt]
   push_{\nabla\Gamma}(A,\, C,\, p \lmerge q) & = &
     \mathsf{allow}(A,\, \Gamma_C(\mathsf{allow}(C^{-1}(A),\, p' \lmerge q')))
     \ \text{where}
     \left\{
     \begin{array}{lll}
       p'  & = & push_{\nabla\Gamma}(A',\, C,\, p) \\
       q'  & = & push_{\nabla\Gamma}(A'',\, C,\, q) \\
       A'  & = & C^{-1}(A)^{\subseteq} \setminus (C^{-1}(A) \setminus A) \\
       A'' & = & (C^{-1}(A) \leftarrowtail \alpha(p')) \setminus (C^{-1}(A) \setminus A)
     \end{array}
     \right. \\[8pt]
   push_{\nabla\Gamma}(A,\, C,\, p \mid q) & = &
     \mathsf{allow}(A,\, \Gamma_C(\mathsf{allow}(C^{-1}(A),\, p' \mid q')))
     \ \text{where}
     \left\{
     \begin{array}{lll}
       p'  & = & push_{\nabla\Gamma}(A',\, C,\, p) \\
       q'  & = & push_{\nabla\Gamma}(A'',\, C,\, q) \\
       A'  & = & C^{-1}(A)^{\subseteq} \setminus (C^{-1}(A) \setminus A) \\
       A'' & = & (C^{-1}(A) \leftarrowtail \alpha(p')) \setminus (C^{-1}(A) \setminus A)
     \end{array}
     \right. \\[8pt]
   push_{\nabla\Gamma}(A,\, C,\, \partial_B(p)) & = &
     push_{\nabla\Gamma}(\partial_B(A),\, C,\, p) \\
   push_{\nabla\Gamma}(A,\, C,\, \nabla_V(p)) & = &
     push_{\nabla\Gamma}(A \cap V,\, C,\, p) \\
   push_{\nabla\Gamma}(A,\, C,\, p) & = &
     \mathsf{allow}(A,\, \Gamma_C(p'))
     \ \text{where } p' = push_{\nabla}(C^{-1}(A),\, p)
     \text{ for all other cases of } p
   \end{array}

Note that in this case the allow set :math:`A` has the general shape
:math:`(A_1^{\subseteq} \setminus A_2^{\subseteq}) I^*` (?), with the subset operator
:math:`\subseteq` optional, and with :math:`I` possibly empty. To implement this
optimization, it needs to be investigated if such a set :math:`A` is closed under the
operations :math:`\partial_B(A)`, :math:`\tau_{I_1}^{-1}(A)`, :math:`A \cap V`,
:math:`R^{-1}(A)`, :math:`C^{-1}(A)`, :math:`A \leftarrowtail A_1`,
:math:`A^{\subseteq}` and :math:`C(A)`.
