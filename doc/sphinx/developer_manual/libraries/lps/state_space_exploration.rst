.. _lps-state-space-exploration:

State Space Exploration
=======================

| *Author*: Wieger Wesselink

Graph Exploration
-----------------

State space exploration is an instance of graph exploration. Consider a directed graph
and take a node :math:`s_0`. We assume there is a function :math:`successors` that returns the
successor nodes of a vertex. An abstract algorithm for exploring the graph starting from vertex
:math:`s_0` is

**Algorithm: Graph exploration**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreGraph}{$s_0$}
   \State $todo := \{ s_0 \}$
   \State $discovered := \{ s_0 \}$
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ s \in todo$
     \State $todo := todo \setminus \{s\}$
     \State $discovered := discovered \cup \{s\}$
     \For{$s' \in successors(s)$}
       \If{$s' \notin discovered$}
         \State $discovered := discovered \cup \{s' \}$
         \State $todo := todo \cup \{ s' \}$
       \EndIf
     \EndFor
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Event points
~~~~~~~~~~~~

There are many different applications of state space exploration. The Boost Graph Library
[Siek02]_ uses a clever idea to separate such applications from the exploration
itself. It is done by distinguishing *event points* in the algorithm that the user can respond
to by means of callback functions. For our purposes we select the following events:

.. _table-eventpoints:

.. list-table:: Event points
   :widths: 30 70

   * - :math:`\textsf{discover\_state}`
     - is invoked when a state is encountered for the first time
   * - :math:`\textsf{examine\_transition}`
     - is invoked on every transition
   * - :math:`\textsf{start\_state}`
     - is invoked on a state right before its outgoing transitions are being explored
   * - :math:`\textsf{finish\_state}`
     - is invoked on a state after all of its outgoing transitions have been explored

The events are named in terms of states and transitions instead of vertices and edges, since this
is closer to our application domain. The exploration algorithm with event points included looks
like this:

**Algorithm: Graph exploration with event points**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreGraph}{$s_0, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $todo := \{ s_0 \}$
   \State $discovered := \{ s_0 \}$
   \State \colorbox{lightgray}{$\textsf{discover\_state}(s_0)$}
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ s \in todo$
     \State $todo := todo \setminus \{s\}$
     \State \colorbox{lightgray}{$\textsf{start\_state}(s)$}
     \State $discovered := discovered \cup \{s\}$
     \For{$s' \in successors(s)$}
       \If{$s' \notin discovered$}
         \State $discovered := discovered \cup \{s' \}$
         \State \colorbox{lightgray}{$\textsf{discover\_state}(s')$}
         \State $todo := todo \cup \{ s' \}$
       \EndIf
       \State \colorbox{lightgray}{$\textsf{examine\_transition}(s, a, s')$}
     \EndFor
     \State \colorbox{lightgray}{$\textsf{finish\_state}(s)$}
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Applications
------------

Many applications can be easily expressed in terms of the given event points.

Deadlock checking
~~~~~~~~~~~~~~~~~

With deadlock checking we are looking for states that have no outgoing transitions.
By introducing one boolean variable :math:`\textsf{has\_transitions}` we can implement deadlock
checking as follows. The callback functions are printed as comments in gray.

**Algorithm: Deadlock checking implemented using event points**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{FindDeadlock}{$s_0, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State \Comment{\colorbox{lightgray}{bool has\_transitions}}
   \State $todo := \{ s_0 \}$
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ s \in todo$
     \State $todo := todo \setminus \{s\}$
     \State $\textsf{start\_state}(s)$ \Comment{\colorbox{lightgray}{has\_transitions := false}}
     \State $discovered := discovered \cup \{s\}$
     \For{$s' \in successors(s)$}
       \If{$s' \notin discovered$}
         \State $discovered := discovered \cup \{s' \}$
         \State $\textsf{discover\_state}(s')$
         \State $todo := todo \cup \{ s' \}$
       \EndIf
       \State $\textsf{examine\_transition}(s, a, s')$ \Comment{\colorbox{lightgray}{has\_transitions := true}}
     \EndFor
     \State $\textsf{finish\_state}(s)$ \Comment{\colorbox{lightgray}{if (!has\_transitions) report\_deadlock(s)}}
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Search strategies
-----------------

Exploration can be done with different search strategies. We describe three of them: breadth-first,
depth-first and highway. They mainly differ in the order in which the elements of the todo set are
processed. In breadth-first search nodes at the present depth are explored before nodes at a higher
depth. In depth-first search the highest-depth nodes are explored first. Highway search is a variant
that uses a breadth-first search, but it only explores a part of the state space.

In all three cases the :math:`todo` list is stored in a double ended queue. We use the slicing
operator to denote parts of a list. For example, :math:`A[m:n]` corresponds to the sublist
:math:`A[m, \ldots, n-1]`.

Breadth-first search
~~~~~~~~~~~~~~~~~~~~

**Algorithm: Breadth-first search**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreGraphBreadthFirst}{$s_0, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $todo :=$ \colorbox{lightgray}{$[s_0]$}
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \While{\colorbox{lightgray}{$|todo| > 0$}}
     \State \colorbox{lightgray}{$s := todo[0]$}
     \State \colorbox{lightgray}{$todo := todo[1:|todo|]$}
     \State $\textsf{start\_state}(s)$
     \State $discovered := discovered \cup \{s\}$
     \For{$s' \in successors(s)$}
       \If{$s' \notin discovered$}
         \State $discovered := discovered \cup \{s' \}$
         \State $\textsf{discover\_state}(s')$
         \State \colorbox{lightgray}{$todo := todo \concat [s']$}
       \EndIf
       \State $\textsf{examine\_transition}(s, a, s')$
     \EndFor
     \State $\textsf{finish\_state}(s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Depth-first search
~~~~~~~~~~~~~~~~~~

**Algorithm: Depth-first search**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreGraphDepthFirst}{$s_0, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $todo :=$ \colorbox{lightgray}{$[s_0]$}
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \While{\colorbox{lightgray}{$|todo| > 0$}}
     \State \colorbox{lightgray}{$s := todo[|todo| - 1]$}
     \State \colorbox{lightgray}{$todo := todo[0:|todo| - 1]$}
     \State $\textsf{start\_state}(s)$
     \State $discovered := discovered \cup \{s\}$
     \For{$s' \in successors(s)$}
       \If{$s' \notin discovered$}
         \State $discovered := discovered \cup \{s' \}$
         \State $\textsf{discover\_state}(s')$
         \State \colorbox{lightgray}{$todo := todo \concat [s']$}
       \EndIf
       \State $\textsf{examine\_transition}(s, a, s')$
     \EndFor
     \State $\textsf{finish\_state}(s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Highway search
~~~~~~~~~~~~~~

In highway search (see [EGWW09]_) a breadth first search is done, with the restriction that at most
:math:`N` states are put in the todo list for each level. The variable :math:`L` maintains the
number of states in the todo list corresponding to the current level, and the variable :math:`c`
counts how many elements have been added corresponding to the next level. Once :math:`c` reaches
the maximum value :math:`N`, elements are being overwritten randomly.

.. note::

   The specification below deviates from the published version of highway search in the
   sense that overwritten elements are added to the set :math:`discovered`. To avoid this,
   the structure of the algorithm needs to be changed significantly.

**Algorithm: Highway search**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreGraphHighway}{$s_0,$ \colorbox{lightgray}{$N$}$, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $todo := [s_0]$
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \State \colorbox{lightgray}{$L := |todo|$}
   \State \colorbox{lightgray}{$c := 0$}
   \While{$|todo| > 0$}
     \State $s := todo[0]$
     \State $todo := todo[1:|todo|]$
     \State $\textsf{start\_state}(s)$
     \For{$s' \in successors(s)$}
       \If{$s' \notin discovered$}
         \State $discovered := discovered \cup \{s' \}$
         \State $\textsf{discover\_state}(s')$
         \State \colorbox{lightgray}{$c := c + 1$}
         \If{\colorbox{lightgray}{$c \leq N$}}
           \State \colorbox{lightgray}{$todo := todo \concat [s']$}
         \Else
           \State \colorbox{lightgray}{$k := random(\{ 1, \ldots, c\})$}
           \If{\colorbox{lightgray}{$k \leq N$}}
             \State \colorbox{lightgray}{$todo[|todo| - k] := s'$}
           \EndIf
         \EndIf
       \EndIf
       \State $\textsf{examine\_transition}(s, a, s')$
     \EndFor
     \State $\textsf{finish\_state}(s)$
     \State \colorbox{lightgray}{$L := L - 1$}
     \If{\colorbox{lightgray}{$L = 0$}}
       \State \colorbox{lightgray}{$L := |todo|$}
       \State \colorbox{lightgray}{$c := 0$}
     \EndIf
   \EndWhile
   \EndProcedure
   \end{algorithmic}

In the algorithm of [EGWW09]_, the set :math:`Q_d` stores todo elements corresponding to the
current level, and the set :math:`Q_{d+1}` stores todo elements corresponding to the next level.
The algorithm above uses only one list :math:`todo` that stores both of them. At each iteration of
the while loop the first :math:`L` elements of the :math:`todo` list belong to the current level,
and the remaining elements belong to the next level. Furthermore, the algorithm above contains only
one application of a random generator, compared to two applications in the original version. The
element :math:`k` is chosen randomly in the range :math:`[1, \ldots, c]`. There is an :math:`N/c`
probability that this value is in the range :math:`[1, \ldots, N]`. If :math:`k` is inside the
range, the element in the :math:`todo` list with index :math:`k` (counting from the end) is
overwritten. This behaviour matches with the published version.

Cycle detection
---------------

For cycle detection the event points in :ref:`table-eventpoints` are insufficient. In [Siek02]_
the following recursive depth first algorithm is given:

**Algorithm: Recursive cycle detection algorithm as specified in Boost**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{boost\_dfs\_recursive}{$u$}
   \State $color[u] := gray$
   \State $\textsf{discover\_vertex}(u)$
   \For{$(a,v) \in out\_edges(u)$}
     \State $\textsf{examine\_edge}(a, v)$
     \If{$color[v] = white$}
       \State $\textsf{tree\_edge}(a, v)$
       \State \Call{dfs\_recursive}{$v$}
     \ElsIf{$color[v] = gray$}
       \State $\textsf{back\_edge}(a, v)$
     \Else
       \State $\textsf{forward\_or\_cross\_edge}(a, v)$
     \EndIf
     \State $color[u] := black$
     \State $\textsf{finish\_vertex}(u)$
   \EndFor
   \EndProcedure
   \end{algorithmic}

The code in Boost uses an iterative version:

**Algorithm: Iterative cycle detection algorithm as implemented in Boost**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{boost\_dfs\_iterative}{$u$}
   \State $color[u] := gray$
   \State $\textsf{discover\_vertex}(u)$
   \State $stack := [(u, out\_edges(u))]$
   \While{$|stack| > 0$}
     \State $u, E := stack.pop\_back()$
     \While{$|E| > 0$}
       \State $a, v := E[0]$
       \State $\textsf{examine\_edge}(u, a, v)$
       \If{$color[v] = white$}
         \State $\textsf{tree\_edge}(u, a, v)$
         \State $stack.push\_back(u, E[1:])$
         \State $u := v$
         \State $color[u] := gray$
         \State $\textsf{discover\_vertex}(u)$
         \State $E := out\_edges(u)$
       \Else
         \If{$color[v] = gray$}
           \State $\textsf{back\_edge}(u, a, v)$
         \Else
           \State $\textsf{forward\_or\_cross\_edge}(u, a, v)$
         \EndIf
         \State $E := E[1:]$
       \EndIf
     \EndWhile
     \State $color[u] := black$
     \State $\textsf{finish\_vertex}(u)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

For our purposes we rewrite this as:

**Algorithm: Recursive cycle detection**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{dfs\_recursive}{$s_0, gray$}
   \State $gray := gray \cup \{ s_0 \}$
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \For{$(a, s_1) \in out\_edges(s_0)$}
     \State $\textsf{examine\_edge}(s_0, a, s_1)$
     \If{$s_1 \notin discovered$}
       \State $\textsf{tree\_edge}(s_0, a, s_1)$
       \State $discovered := discovered \cup \{ s_1 \}$
       \State \Call{dfs\_recursive}{$s_1, gray$}
     \ElsIf{$s_1 \in todo$}
       \State $\textsf{back\_edge}(s_0, a, s_1)$
     \Else
       \State $\textsf{forward\_or\_cross\_edge}(s_0, a, s_1)$
     \EndIf
     \State $gray := gray \setminus \{ s_0 \}$
     \State $\textsf{finish\_state}(s_0)$
   \EndFor
   \EndProcedure
   \end{algorithmic}

**Algorithm: Iterative cycle detection**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{dfs\_iterative}{$s_0$}
   \State $todo := [(s_0, out\_edges(s_0))]$
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \While{$|todo| > 0$}
     \State $s, E := todo.back()$
     \While{$|E| > 0$}
       \State $a, s_1 := E.pop\_front()$
       \State $\textsf{examine\_edge}(s_0, a, s_1)$
       \If{$s_1 \notin discovered$}
         \State $\textsf{tree\_edge}(s_0, a, s_1)$
         \State $discovered := discovered \cup \{ s_1 \}$
         \State $\textsf{discover\_state}(s_1)$
         \State $todo.back() := (s, E)$
         \State $todo := todo \concat [(s_1, out\_edges(s_1))]$
         \State $s, E := todo.back()$
       \ElsIf{$s_1 \in todo$}
         \State $\textsf{back\_edge}(s_0, a, s_1)$
       \Else
         \State $\textsf{forward\_or\_cross\_edge}(s_0, a, s_1)$
       \EndIf
     \EndWhile
     \State $\textsf{finish\_state}(s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Whenever the :math:`\textsf{back\_edge}` event is triggered, a cycle is found.

Untimed state space exploration
--------------------------------

Consider the following untimed linear process specification :math:`P`, with initial state :math:`d_0`.

.. math::

   \begin{array}{l}
   P(d)=
   \sum\limits_{i\in I}\sum\limits_{e_i}c_i(d, e_i)\rightarrow a_i(f_i(d,e_i)) \cdot P(g_i(d,e_i))
   \end{array}

This linear process is a symbolic representation of a state space, or labeled transition system
(LTS). The previously described graph exploration algorithms can be applied to explore a state space.
Let :math:`rewr` be a rewriter. An algorithm for untimed state space exploration is

**Algorithm: Untimed LPS exploration**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreLPS}{$P(d), d_0, rewr, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $s_0 := rewr(d_0, \emptylist)$
   \State $todo := \{ s_0 \}$
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ s \in todo$
     \State $todo := todo \setminus \{s\}$
     \State $discovered := discovered \cup \{s\}$
     \State $\textsf{start\_state}(s)$
     \For{$i \in I$}
       \State $condition := rewr(c_i(d, e_i), [d := s])$
       \If{$condition = false$}
         \State \textbf{continue}
       \EndIf
       \State $E := \{ e \mid rewr(condition, [e_i := e]) = true \}$
       \For{$e \in E$}
         \State $a := a_i(rewr(f_i(d,e_i), [d:=s,e_i:=e]))$
         \State $s' := rewr(g_i(d,e_i), [d:=s,e_i:=e])$
         \If{$s' \notin discovered$}
           \State $todo := todo \cup \{ s' \}$
           \State $discovered := discovered \cup \{s'\}$
           \State $\textsf{discover\_state}(s')$
         \EndIf
         \State $\textsf{examine\_transition}(s, a, s')$
       \EndFor
     \EndFor
     \State $\textsf{finish\_state}(s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

The set :math:`E` is computed using the Enumerate algorithm. This computation may be expensive.
Hence the condition :math:`c(d,e_i)` is first rewritten, since if it evaluates to :math:`false`
the computation of :math:`E` can be skipped.

Timed state space exploration
------------------------------

Consider the following timed linear process specification :math:`P`, with initial state :math:`d_0`.

.. math::

   \begin{array}{l}
   P(d)=
   \sum\limits_{i\in I}\sum\limits_{e_i}c_i(d, e_i)\rightarrow a_i(f_i(d,e_i))
   \mbox{\colorbox{lightgray}{$\at{t_i(d,e_i)}$}}
   \cdot P(g_i(d,e_i)).
   \end{array}

Note that the time tag :math:`t_i(d,e_i)` is optional. If it is omitted, the corresponding action
may happen at an arbitrary time. In timed state space exploration, care is taken that on every trace
the time tags are increasing. In order to achieve that, a time stamp is recorded for each state in
the state space. We use the notation :math:`t \aftertime s` to denote the state :math:`s` with
associated time stamp :math:`t`. An algorithm for timed state space exploration is

**Algorithm: Timed LPS exploration**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreLPSTimed}{$P(d), d_0, rewr, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $s_0 := rewr(d_0, \emptylist)$
   \State $todo := \{ \hlmath{0 \aftertime s_0} \}$
   \State $discovered := \{ \hlmath{0 \aftertime s_0} \}$
   \State $\textsf{discover\_state}(0 \aftertime s_0)$
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ \hlmath{t \aftertime s} \in todo$
     \State $todo := todo \setminus \{ \hlmath{t \aftertime s} \}$
     \State $discovered := discovered \cup \{ \hlmath{t \aftertime s} \}$
     \State $\textsf{start\_state}(t \aftertime s)$
     \For{$i \in I$}
       \State $condition := rewr(c_i(d, e_i), [d := s])$
       \If{$condition = false$}
         \State \textbf{continue}
       \EndIf
       \State $E := \{ e \mid rewr(condition, [e_i := e]) = true \}$
       \For{$e \in E$}
         \State \colorbox{lightgray}{$t' := rewr(t_i(d,e_i), [d:=s,e_i:=e])$}
         \If{$\hlmath{t' \leq t}$}
           \State \textbf{continue}
         \EndIf
         \State $a := a_i(rewr(f_i(d,e_i), [d:=s,e_i:=e]))$
         \State $s' := rewr(g_i(d,e_i), [d:=s,e_i:=e])$
         \If{$\hlmath{t' \aftertime s'} \notin discovered$}
           \State $todo := todo \cup \{ \hlmath{t' \aftertime s'} \}$
           \State $discovered := discovered \cup \{ \hlmath{t' \aftertime s'} \}$
           \State $\textsf{discover\_state}(t' \aftertime s')$
         \EndIf
         \State $\textsf{examine\_transition}(t \aftertime s, a\at{t'}, t' \aftertime s')$
       \EndFor
     \EndFor
     \State $\textsf{finish\_state}(t \aftertime s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Stochastic state space exploration
------------------------------------

Consider the following stochastic linear process specification :math:`P`, with initial state
:math:`\frac{p(h)}{h} \cdot P(g(h))`.

.. math::

   \begin{array}{l}
   P(d)=
   \sum\limits_{i\in I}\sum\limits_{e_i}c_i(d, e_i)\rightarrow a_i(f_i(d,e_i))
   \mbox{\colorbox{lightgray}{$\dfrac{p_i(d,e_i,h_i)}{h_i}$}}
   \cdot P(g_i(d,e_i,h_i)),
   \end{array}

where :math:`p` and :math:`p_i` are stochastic distributions.
We define a *stochastic state* as a set
:math:`\{(q_1, s_1), \ldots, (q_m, s_m)\}` with :math:`q_j,\ j = 1 \ldots m` a sequence of
non-zero probabilities that sum up to 1, and :math:`s_j,\ j = 1 \ldots m` a sequence of states.
The function :math:`\textsc{ComputeStochasticState}` is used to compute a stochastic state from
its symbolic representation.

**Algorithm: Computation of a stochastic state**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ComputeStochasticState}{$h, p, g, rewr, \sigma$}
   \State $result := \emptyset$
   \State $H := \{ (h',q) \mid q = rewr(p, \sigma[h := h']) \wedge q > 0 \}$
   \For{$(h', q) \in H$}
     \State $s := rewr(g, \sigma[h := h'])$
     \State $result := result \cup \{ (q, s) \}$
   \EndFor
   \State \Return $result$
   \EndProcedure
   \end{algorithmic}

The set :math:`H` is computed using the Enumerate algorithm.

An algorithm for stochastic state space exploration is

**Algorithm: Stochastic LPS exploration**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreLPSStochastic}{$P(d), \frac{p(h)}{h} \cdot P(g(h)), rewr, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}, \textsf{discover\_initial\_state}$}
   \State \colorbox{lightgray}{$\hat{s_0} := \textsc{ComputeStochasticState}(h, p(h), g(h), rewr, \emptylist)$}
   \State \colorbox{lightgray}{$S := \{ s_i \mid (q_i, s_i) \in \hat{s_0} \}$}
   \State \colorbox{lightgray}{$\textsf{discover\_initial\_state}(\hat{s_0})$}
   \For{\colorbox{lightgray}{$s \in S$}}
     \State $todo := todo \cup \{ s \}$
     \State $discovered := discovered \cup \{ s \}$
     \State $\textsf{discover\_state}(s)$
   \EndFor
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ s \in todo$
     \State $todo := todo \setminus \{s\}$
     \State $discovered := discovered \cup \{s\}$
     \State $\textsf{start\_state}(s)$
     \For{$i \in I$}
       \State $condition := rewr(c_i(d, e_i), [d := s])$
       \If{$condition = false$}
         \State \textbf{continue}
       \EndIf
       \State $E := \{ e \mid rewr(condition, [e_i := e]) = true \}$
       \For{$e \in E$}
         \State $a := a_i(rewr(f_i(d,e_i), [d:=s,e_i:=e]))$
         \State \colorbox{lightgray}{$\hat{s'} := \textsc{ComputeStochasticState}(h_i, p_i(d,e_i,h_i), g_i(d,e_i,h_i), rewr, [d:=s, e_i:=e])$}
         \State \colorbox{lightgray}{$S' := \{ s_i \mid (q_i, s_i) \in \hat{s'} \}$}
         \For{\colorbox{lightgray}{$s' \in S'$}}
           \If{$s' \notin discovered$}
             \State $todo := todo \cup \{ s' \}$
             \State $discovered := discovered \cup \{s'\}$
             \State $\textsf{discover\_state}(s')$
           \EndIf
         \EndFor
         \State $\textsf{examine\_transition}(s, a,$ \colorbox{lightgray}{$\hat{s'}$}$)$
       \EndFor
     \EndFor
     \State $\textsf{finish\_state}(s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Caching
-------

The computation of the set of solutions :math:`E` in :math:`\textsc{ExploreLPS}` is expensive.
Therefore it may be a good idea to cache these solutions. Caching can be done locally (i.e. using a
separate cache for each summand), or globally. This leads to the following variants of the
algorithm. We assume that :math:`FV` is a function that computes free variables of an expression.
Let :math:`\mathcal{D}` be the set of process parameters (i.e. the elements of :math:`d`).

Local caching
~~~~~~~~~~~~~

In the local caching algorithm for each summand :math:`i` a mapping :math:`C_i` is maintained.
The cache key is comprised of the actual values of the process parameters that appear in the
condition :math:`c_i(d, e_i)`.

**Algorithm: LPS exploration with local caching**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreLPSLocallyCached}{$P(d), d_0, rewr, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $s_0 := rewr(d_0, \emptylist)$
   \State $todo := \{ s_0 \}$
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \For{\colorbox{lightgray}{$i \in I$}}
     \State \colorbox{lightgray}{$C_i := \emptymap$}
     \State \colorbox{lightgray}{$\gamma_i := FV(c_i(d, e_i)) \cap \mathcal{D}$}
   \EndFor
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ s \in todo$
     \State $todo := todo \setminus \{s\}$
     \State $discovered := discovered \cup \{s\}$
     \State $\textsf{start\_state}(s)$
     \For{$i \in I$}
       \State \colorbox{lightgray}{$key := \gamma_i[d:=s]$}
       \If{\colorbox{lightgray}{$key \in keys(C_i)$}}
         \State \colorbox{lightgray}{$E := C_i[key]$}
       \Else
         \State \colorbox{lightgray}{$E := \{ e \mid rewr(c_i(d, e_i), [d:=s,e_i:=e]) = true \}$}
         \State \colorbox{lightgray}{$C_i := C_i \cup \{(key, E)\}$}
       \EndIf
       \For{$e \in E$}
         \State $a := a_i(rewr(f_i(d,e_i), [d:=s,e_i:=e]))$
         \State $s' := rewr(g_i(d,e_i), [d:=s,e_i:=e])$
         \If{$s' \notin discovered$}
           \State $todo := todo \cup \{ s' \}$
           \State $discovered := discovered \cup \{s'\}$
           \State $\textsf{discover\_state}(s')$
         \EndIf
         \State $\textsf{examine\_transition}(s, a, s')$
       \EndFor
     \EndFor
     \State $\textsf{finish\_state}(s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

Global caching
~~~~~~~~~~~~~~

In the global caching algorithm one mapping :math:`C` is maintained. To achieve this, the condition
of the summands is added to the cache key. If many summands share the same condition, global caching
may be beneficial. In practice this doesn't seem to happen much.

**Algorithm: LPS exploration with global caching**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreLPSGloballyCached}{$P(d), d_0, rewr, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $todo := \{d_0\}$
   \State $discovered := \emptyset$
   \State \colorbox{lightgray}{$C := \emptyset$}
   \For{\colorbox{lightgray}{$i \in I$}}
     \State \colorbox{lightgray}{$\gamma_i := FV(c_i(d, e_i)) \cap \mathcal{D}$}
   \EndFor
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ s \in todo$
     \State $todo := todo \setminus \{s\}$
     \State $discovered := discovered \cup \{s\}$
     \State $\textsf{start\_state}(s)$
     \For{$i \in I$}
       \State \colorbox{lightgray}{$key := c_i(d,e_i) \concat \gamma_i[d:=s]$}
       \If{\colorbox{lightgray}{$key \in keys(C)$}}
         \State \colorbox{lightgray}{$T := C[key]$}
       \Else
         \State \colorbox{lightgray}{$T := \{ t \mid rewr(c_i(d, e_i), [d:=s,e_i:=t]) = true \}$}
         \State \colorbox{lightgray}{$C := C \cup \{(key, T)\}$}
       \EndIf
       \For{$e \in E$}
         \State $a := a_i(rewr(f_i(d,e_i), [d:=s,e_i:=e]))$
         \State $s' := rewr(g_i(d,e_i), [d:=s,e_i:=e])$
         \If{$s' \notin discovered$}
           \State $todo := todo \cup \{ s' \}$
           \State $discovered := discovered \cup \{s'\}$
           \State $\textsf{discover\_state}(s')$
         \EndIf
         \State $\textsf{examine\_transition}(s, a, s')$
       \EndFor
     \EndFor
     \State $\textsf{finish\_state}(s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

In this algorithm :math:`C` is a mapping, with :math:`keys(C) = \{ k \mid \exists_{v}: (k,v) \in C \}`.
We use the notation :math:`C[k]` to denote the unique element :math:`v` such that :math:`(k,v) \in C`.

Confluence Reduction
--------------------

Confluence reduction (see [GP00]_, [Blo01]_ and [BP02]_) is an on-the-fly state space exploration
method that produces a reduced state space. For confluence reduction we assume that the set of
summands :math:`I` is partitioned into a set :math:`I_{regular}` of 'regular' summands, and a set
:math:`I_{confluent}` of confluent :math:`\tau`-summands. The confluent :math:`\tau`-summands are
used to determine a unique representative state that is reachable via confluent :math:`\tau` steps.
This is done using the graph algorithm :math:`\textsc{FindRepresentative}`. This leads to the
following variant of the algorithm:

**Algorithm: LPS exploration with confluence reduction**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{ExploreLPSConfluence}{$P(d), d_0, rewr, \textsf{discover\_state}, \textsf{examine\_transition}, \textsf{start\_state}, \textsf{finish\_state}$}
   \State $s_0 := \hlmath{\textsc{FindRepresentative}(rewr(d_0, \emptylist))}$
   \State $todo := \{ s_0 \}$
   \State $discovered := \{ s_0 \}$
   \State $\textsf{discover\_state}(s_0)$
   \While{$todo \neq \emptyset$}
     \State $\textbf{choose}\ s \in todo$
     \State $todo := todo \setminus \{s\}$
     \State $discovered := discovered \cup \{s\}$
     \State $\textsf{start\_state}(s)$
     \For{$i \in \hlmath{I_{regular}}$}
       \State $condition := rewr(c_i(d, e_i), [d := s])$
       \If{$condition = false$}
         \State \textbf{continue}
       \EndIf
       \State $E := \{ e \mid rewr(condition, [e_i := e]) = true \}$
       \For{$e \in E$}
         \State $a := a_i(rewr(f_i(d,e_i), [d:=s,e_i:=e]))$
         \State $s' := \hlmath{\textsc{FindRepresentative}(rewr(g_i(d,e_i), [d:=s,e_i:=e]))}$
         \If{$s' \notin discovered$}
           \State $todo := todo \cup \{ s' \}$
           \State $discovered := discovered \cup \{s'\}$
           \State $\textsf{discover\_state}(s')$
         \EndIf
         \State $\textsf{examine\_transition}(s, a, s')$
       \EndFor
     \EndFor
     \State $\textsf{finish\_state}(s)$
   \EndWhile
   \EndProcedure
   \end{algorithmic}

As suggested in [BP02]_ Tarjan's strongly connected component (SCC) algorithm (see [Tar72]_) can be
used to compute a unique representative.

Tarjan's SCC algorithm
~~~~~~~~~~~~~~~~~~~~~~

A recursive implementation of Tarjan's strongly connected components algorithm that uses four global
variables :math:`stack`, :math:`low`, :math:`disc` and :math:`result`. The helper function
:math:`\textsc{StrongConnect}` computes the connected component reachable from node :math:`u`. In
this function it is assumed that the function call :math:`\textsf{successors}(u)` returns the
successor states of :math:`u` in a deterministic order.

**Algorithm: Tarjan's Strongly Connected Component Algorithm**

*Input:* :math:`G=(V,E)`: A graph with nodes :math:`V` and edges :math:`E`.

*Output:* :math:`result`: A sequence containing all strongly connected components of the graph.

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{Tarjan}{$G$}
   \State $stack := \emptylist$
   \State $low := \emptymap$ \Comment{the empty mapping is denoted as $\emptymap$}
   \State $disc := \emptymap$ \Comment{$u \in low$ means $u$ is a key of mapping $low$}
   \State $result := \emptylist$
   \For{$u \in V$}
     \If{$u \notin low$}
       \State \Call{StrongConnect}{$u$}
     \EndIf
   \EndFor
   \State \Return $result$
   \EndProcedure
   \end{algorithmic}

**Algorithm: Helper function StrongConnect**

*Input:* :math:`u`: An element of :math:`V`.

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{StrongConnect}{$u$}
   \State $k := |disc|$ \Comment{$k$ is the discovery time assigned to node $u$}
   \State $disc[u] := k$
   \State $low[u] := k$ \Comment{initially $low[u] = disc[u]$}
   \State $stack := stack \concat [u]$
   \For{$v \in \textsf{successors}(u)$}
     \If{$v \notin low$}
       \State \Call{StrongConnect}{$v$}
       \State $low[u] := \textsf{min}(low[u], low[v])$
     \ElsIf{$v \in stack$}
       \State $low[u] := \textsf{min}(low[u], disc[v])$
     \EndIf
   \EndFor
   \If{$low[u] = disc[u]$} \Comment{an SCC has been found}
     \State $comp := \emptylist$
     \While{\textbf{true}}
       \State $v := stack[|stack| - 1]$ \Comment{assign the top of the stack to $v$}
       \State $stack := stack[0:|stack| - 1]$ \Comment{pop an element from the stack}
       \State $comp := comp \concat [v]$
       \If{$v = u$}
         \State \textbf{break}
       \EndIf
     \EndWhile
     \State $result := result \concat [comp]$
   \EndIf
   \EndProcedure
   \end{algorithmic}

A side effect of a call :math:`\textsc{Tarjan}(u)` is that :math:`result` contains the connected
components that have been found.

FindRepresentative
~~~~~~~~~~~~~~~~~~

Due to properties of confluent :math:`\tau`-summands, there is always only one terminal strongly
connected component, i.e. a strongly connected component without outgoing edges. Furthermore, the
first strongly connected component reported by Tarjan's algorithm is always terminating. For our
implementation of :math:`\textsc{FindRepresentative}` we prefer to use an iterative version of
Tarjan's SCC algorithm. The reason for this is that an iterative version can be more easily
interrupted once the first SCC has been found. The algorithm description in [Oqv]_ has been used
as a model for our solution.

**Algorithm: Find a unique representative node in a graph**

.. math::
   :nowrap:
   :class: math-left

   \begin{algorithmic}[1]
   \Procedure{FindRepresentative}{$u$}
   \State $stack := \emptylist$
   \State $low := \emptymap$
   \State $disc := \emptymap$
   \State $work := [(u, 0)]$
   \While{$work \neq \emptylist$}
     \State $(u, i) := work[|work| - 1]$
     \State $work := work[0 : |work| - 1]$
     \If{$i = 0$}
       \State $k := |disc|$
       \State $disc[u] := k$
       \State $low[u] := k$
       \State $stack := stack \concat [u]$
     \EndIf
     \State $recurse := false$
     \For{$j \in [i, \ldots, |successors(u)|]$}
       \State $v := successors(u)[j]$
       \If{$v \notin disc$}
         \State $work := work \concat [(u, j+1)]$
         \State $work := work \concat [(v, 0)]$
         \State $recurse := true$
         \State \textbf{break}
       \ElsIf{$v \in stack$}
         \State $low[u] := min(low[u], disc[v])$
       \EndIf
     \EndFor
     \If{$recurse$}
       \State \textbf{continue}
     \EndIf
     \If{$low[u] = disc[u]$}
       \State $result := u$
       \While{\textbf{true}}
         \State $v := stack[|stack| - 1]$
         \State $stack := stack[0:|stack| - 1]$
         \If{$v = u$}
           \State \textbf{break}
         \EndIf
         \If{$v < result$}
           \State $result := v$
         \EndIf
       \EndWhile
       \State \Return $result$
     \EndIf
     \If{$work \neq \emptylist$}
       \State $v := u$
       \State $(u, z) := work[|work| - 1]$
       \State $low[u] := min(low[u], low[v])$
     \EndIf
   \EndWhile
   \EndProcedure
   \end{algorithmic}

References
----------

.. [Siek02] Jeremy Siek, Lie-Quan Lee and Andrew Lumsdale.
   *The Boost Graph Library: User Guide and Reference Manual*.
   Addison-Wesley, 2002.

.. [EGWW09] Tom A. N. Engels, Jan Friso Groote, Muck van Weerdenburg and Tim A. C. Willemse.
   Search algorithms for automated validation.
   *Journal of Logic and Algebraic Programming*, 78(4):274–287, 2009.

.. [GP00] Jan Friso Groote and Jaco van de Pol.
   State Space Reduction Using Partial tau-Confluence.
   In *MFCS*, LNCS 1893, pages 383–393. Springer, 2000.

.. [Blo01] Stefan Blom.
   Partial t-confluence for Efficient State Space Generation, 2001.

.. [BP02] Stefan Blom and Jaco van de Pol.
   State Space Reduction by Proving Confluence.
   In *CAV*, LNCS 2404, pages 596–609. Springer, 2002.

.. [Tar72] Robert Tarjan.
   Depth first search and linear graph algorithms.
   *SIAM Journal on Computing*, 1(2), 1972.

.. [Oqv] Jesper Öqvist.
   Iterative Tarjan Strongly Connected Components in Python.
   https://llbit.se/?p=3379. Accessed 2019-03-26.
