Gossips
-------

**Contribution of this section**

#. use of sets,
#. exercise with abstract datatypes, and
#. complicated µ-calculus formulae.

**New tools**:
none.

Initially, each of five teenage schoolgirls Ann, Beth, Carol, Daisy, and Elvy,
knows a unique, new, thrilling gossip. They urge to phone each other to exchange
all the gossips they know, thus sharing the same gossips right after the phone
call. We wonder, what is the minimum number of phone calls that have to take
place such that all five girls know all the five new gossips?

In a naive approach for a solution of sharing all gossips, Ann calls Beth,
Carol, Daisy and Elvy in a row, collecting all four other gossips. Note that
gradually Ann has more gossips to share, because she knows the gossips of the
girls she has called so far. For example, finishing her fourth call both Ann and
Elvy know all gossips. Next, Ann calls three others, say Daisy, Carol and Beth,
to bring them up to date with all the latest news. This takes 7 phone calls in
total.

In our modelling below, we choose to index the five gossips by numbers from ``1``
up to ``5``. If one girl knows the gossips ``1`` and ``4``, and another girl
knows gossips ``2``, ``3`` and ``4``, then the two girls both know the gossips
``1``, ``2``, ``3`` and ``4`` after talking to each other on the phone.

A model using sets of positive integers is given below (also available as
:download:`gossip-blanks.mcrl2 <files/gossip-blanks.mcrl2>`. Dealing with sets,
we first introduce constants ``N`` and ``Gossips`` of type ``Pos`` and
``Set(Pos)``, respectively. For the definition of the set ``Gossips`` we use set
comprehension, ``Gossips`` comprises all positive ``k`` such that ``1 <= k <=
N``.

.. literalinclude:: files/gossip-blanks.mcrl2
   :language: mcrl2

We distinguish the actions ``done``, that indicates that a particular girl has
gathered all gossips, and ``all_done``, that indicates that all girls know all
gossips. Furthermore, we have actions ``call``, ``answer`` and ``exchange`` for
initiating a phone call, picking up the phone, and the synchronized execution of
these two. Note the type of the latter three actions: ``Pos # Pos # Set(Pos) #
Set(Pos)``, two identifiers for the caller and callee, two gossip sets to
reflect their respective knowledge at that moment.

For each of the girls, we define a process ``Girl`` carrying an identifier
``myid``, a positive integer, and the current knowledge of gossips ``mygs``, a
set of positive integers, as parameters. Basically, a girl with identifier
``myid`` can either phone another executing the action ``call`` or answer the
phone of another executing the action ``answer``. As the actions have place
holders for the identifier and gossips of the other girl, the summation ``sum
herid:Pos, hergs:Set(Pos)`` quantifies over all possible values, thereby binding
the variables ``herid`` and ``hergs``. After excuting either a ``call`` or
``answer`` action, the knowledge of a girl owning the process gets updates, it now
incorporates the knowledge of the other girl too.

As termination condition for one ``Girl`` process, we will have equality of the
latest gossip knowledge parameter of the process and the set constant
``Gossips``. If all ``Girl`` processes are done, i.e. can perform the ``done``
action, they can synchronize yielding the action ``all_done``. The delta
explicitly terminates the process.

.. admonition:: Exercise

   Fill in the blanks in :download:`gossip-blanks.mcrl2
   <files/gossip-blanks.mcrl2>`. This means you will have to
   
   * determine the type of the function ``Gossips``,
   * finish the definition of the equation for ``Gossips``, and
   * finish the communication definition.
   
.. admonition:: Solution
   :class: collapse
   
   ``Gossips`` is a constant of type ``Set(Pos)``, with equation
   ``Gossips = Gossips = { k:Pos | k >= 1 && k <= N }``. The communication
   enforces the ``done`` actions to communicate to ``all_done``.
   
   When these modifications have been made, the specification looks as follows
   (:download:`gossip.mcrl2 <files/gossip.mcrl2>`)
   
   .. literalinclude:: files/gossip.mcrl2
      :language: mcrl2

.. admonition:: Exercise

   Generate the state space for :download:`gossip.mcrl2 <files/gossip.mcrl2>`
   
.. admonition:: Hint

   If you think verification takes a long time add the ``-v`` flag
   to :ref:`tool-lps2lts`, this will show you progress messages. If you are on a
   platform other that Windows, you can also pass the ``-rjittyc`` flag
   to :ref:`tool-lps2lts` to use the compiling rewriter. This is more efficient
   than the default jitty rewriter.
   
.. admonition:: Solution
   :class: collapse
   
   We run the following commands::
   
     $ mcrl22lps gossip.mcrl2 gossip.lps
     $ lps2lts -v -rjittyc gossip.lps gossip.lts
       [09:02:44.255 verbose] Detected mCRL2 extension.
       [09:02:44.256 verbose] removing unused parts of the data specification.
       [09:02:44.262 verbose] using 'mcrl2compilerewriter' to compile rewriter.
       [09:02:44.274 verbose] compiling rewriter...
       [09:02:46.764 verbose] loading rewriter...
       [09:02:46.765 verbose] writing state space in lts format to 'gossip.lts'.
       [09:02:46.765 verbose] generating state space with 'breadth' strategy...
       [09:02:46.766 verbose] monitor: level 1 done. (1 state, 20 transitions)
       [09:02:46.769 verbose] monitor: level 2 done. (10 states, 200 transitions)
       [09:02:46.795 verbose] monitor: level 3 done. (75 states, 1500 transitions)
       [09:02:46.971 verbose] monitor: level 4 done. (430 states, 8600 transitions)
       [09:02:47.229 verbose] monitor: currently at level 5 with 1000 states and 20000 transitions explored and 3556 states seen.
       [09:02:47.770 verbose] monitor: currently at level 5 with 2000 states and 40000 transitions explored and 5961 states seen.
       [09:02:47.927 verbose] monitor: level 5 done. (1725 states, 34500 transitions)
       [09:02:48.417 verbose] monitor: currently at level 6 with 3000 states and 60000 transitions explored and 6974 states seen.
       [09:02:49.179 verbose] monitor: currently at level 6 with 4000 states and 80000 transitions explored and 7601 states seen.
       [09:02:49.901 verbose] monitor: currently at level 6 with 5000 states and 100000 transitions explored and 8203 states seen.
       [09:02:50.633 verbose] monitor: currently at level 6 with 6000 states and 120000 transitions explored and 8844 states seen.
       [09:02:50.986 verbose] monitor: level 6 done. (4180 states, 83600 transitions)
       [09:02:51.631 verbose] monitor: currently at level 7 with 7000 states and 140001 transitions explored and 9152 states seen.
       [09:02:52.698 verbose] monitor: currently at level 7 with 8000 states and 160001 transitions explored and 9152 states seen.
       [09:02:53.730 verbose] monitor: currently at level 7 with 9000 states and 180001 transitions explored and 9152 states seen.
       [09:02:53.871 verbose] monitor: level 7 done. (2731 states, 54621 transitions)
       [09:02:53.871 verbose] Starting to save file gossip.lts
       [09:03:16.396 verbose] done with state space generation (7 levels, 9152 states and 183041 transitions)

Next we are going to investigate some properties of the gossiping girls.

.. admonition:: Exercise

   Verify that the model is deadlock free.
   
.. admonition:: Solution
   :class: collapse
   
   There are several approaches to verifying absence of deadlock. We discuss
   three of them.
   
   #. The first approach instructs lps2lts to report deadlocks; if no deadlock
      is reported, the system is deadlock free, using the ``-D`` option.
      We also instruct the tool to save a trace to the deadlock,
      and terminate once a deadlock has been found using the option ``-t1``::
      
        $ lps2lts -rjittyc -D -t1 gossip.lps
        
      This does not report anything, so the specification is deadlock free.
      
   #. Assuming we already have the transition system, say in :file:`gossip.lts`,
      we can also use the marking facilities of :ref:`tool-ltsview` to check for
      absence of deadlock. We perform the following steps::
      
        $ ltsview gossip.lts
        
      Now open the :guilabel:`Mark` dialog through :menuselection:`Tools -->
      Mark...`. In this dialog select :guilabel:`Mark deadlocks`. Now any
      cluster in the state space containing a deadlock will be coloured red, and
      the deadlock state is also coloured red.
      
      This method does not colour any states and clusters for the gossiping
      girls, so the system is deadlock free.
      
   #. The third approach uses the µ-calculus to perform the verification.
      Absence of deadlock can be specified as follows
      (:download:`gossip1.mcf <files/gossip1.mcf>`):
      
      .. literalinclude:: files/gossip1.mcf
         :language: mcrl2
         
      We can now verify this using :ref:`tool-lps2pbes` and
      :ref:`tool-pbessolve`::
      
        $ lps2pbes -f gossip1.mcf gossip.lps gossip1.pbes
        $ pbessolve -rjittyc gossip1.pbes
        true        
      
.. admonition:: Exercse

   It is straightforward to see that for ``N`` gossiping girls, there is a path
   of length ``(N - 1) + (N - 2)`` states leading to the situation where all
   girls know all gossips. Simply have the first girl call all other (``N-1``)
   girls; then the first and the last girl know all gossips. Then the first girl
   calls the ``N - 2`` remaining girls.
   
   Verify this claim for ``N = 5``, i.e. show that for 5 gossiping girls, there
   is a path of ``7`` phone calss to the situation where all girls know all
   gossips.
   
.. admonition:: Solution
   :class: collapse
   
   Again we can pursue various different paths in this verification. We will
   illustrate three of them.
   
   #. The first approach uses the µ-calculus to formalise the property. Because
      we are in the finite case, it is tempting to choose the following
      formalisation (:download:`gossip3.mcf <files/gossip3.mcf>`).
      
      .. literalinclude:: files/gossip3.mcf
         :language: mcrl2
         
      This simply says that there is a path of length 7 to a state in which
      an ``all_done`` action can be performed.
      
      However, if we try to verify this property using::
      
        $ lps2pbes -f gossip3.mcf gossip.lps | pbessolve -rjittyc
        
      It seems that :ref:`tool-lps2pbes` is getting stuck. This is caused by the
      translation of µ-calculus formula with an LPS to a PBES, that has to look
      ahead 8 levels in the state space, by recursively evaluating the guards in
      the LPS. This causes a vast blowup in the computations that are performed
      internally. We can restrict the number of levels that the tool needs to
      look ahead by introducing fixed points as is done in
      (:download:`gossip3a.mcf <files/gossip3a.mcf>`):
      
      .. literalinclude:: files/gossip3a.mcf
         :language: mcrl2
         
      Using this formalisation leads to blow up in the number equations in the
      PBES, hence generating the PBES has become more efficient, but the solving
      process may be slower, but does succeed::
      
        $ lps2pbes -f gossip3a.mcf gossip.lps | pbessolve -rjittyc
        true
   
   #. The µ-calculus formulae in the previous approach must be modified if the
      value of ``N`` is changed, which is, in most cases, undesirable. The
      µ-calculus used in mCRL2 supports parameterisation with data. We can use
      this functionality to write a formula for this property for arbitrary
      values of ``N`` as follows (:download:`gossip2.mcf <files/gossip2.mcf>`).
      
      .. literalinclude:: files/gossip2.mcf
         :language: mcrl2
         
      Here the number of actions that has been performed so far is counted by
      the parameter ``n`` of the formula.
      
      We verify the property using::
      
        $ lps2pbes -f gossip2.mcf gossip.lps | pbessolve -rjittyc
        true
        
      Verification of the property using this formula is quick.
   
   #. As a final approach, we see whether we can use the facilities of
      :ref:`tool-lps2lts` to search for an action. We use the knowledge that
      ``all_done`` is performed when all girls know all gossips. We store traces
      to the occurrences of ``all_done``::
      
        $ lps2lts -rjittyc -aall_done -t gossip.lps
        [13:36:38.405 info]    detect: action 'all_done' found and saved to 'gossip.lps_act_0_all_done.trc'(state index: 6571).
        
      A trace has been saved to :file:`gossip.lps_act_0_all_done.trc`; we
      inspect the trace::
      
        $ tracepp gossip.lps_act_0_all_done.trc
        exchange(4, 5, {4}, {5})
        exchange(3, 4, {3}, {4, 5})
        exchange(2, 1, {2}, {1})
        exchange(2, 3, {1, 2}, {3, 4, 5})
        exchange(3, 5, {1, 2, 3, 4, 5}, {4, 5})
        exchange(4, 1, {3, 4, 5}, {1, 2})
        all_done
        
      We see that the only trace that we stored is six exchanges long, instead
      of the seven exchanges that we are looking for. The previous approach
      proved that there is a path of seven exchanges though. This can be
      explained by the way in which :ref:`tool-lps2lts` generates traces. When
      searching for an ``all_done`` action, the tool will save the trace to the
      state that is reached with the ``all_done`` action. It will, however, only
      store one trace per state that is reached, i.e. if there is another path
      to the state doing the ``all_done`` action, it will not store this trace.
      As a result, this approach cannot be used for the verification task at
      hand.
      
      
.. admonition:: Exercise

   As a last exercise, check whether or not there is a path shorter than
   ``(N-1)+(N-2)`` leading to the situation where all girls know all gossips.
   
.. admonition:: Solution
   :class: collapse
   
   Again various approaches are possible here. We discuss two of them.
   
   #. In the previous exercise, we generated a trace to an ``all_done`` action
      using :ref:`tool-lps2lts`. The trace comprised six exchanges, so we easily
      verify that a trace shorter than 7 exchanges to the situation in which all
      girls know all gossips is possible.
      
   #. The second approach formalises this property using the modal µ-calculus.
      The following, propositional µ-calculus formula
      (:download:`gossip4.mcf <files/gossip4.mcf>`) expresses the property.
      
      .. literalinclude:: files/gossip4.mcf
         :language: mcrl2
         
      We verify this using the following commands::
      
        $ lps2pbes -f gossip4.mcf | pbessolve -rjittyc        
        false
        
      Indeed the property does not hold, as we also observed from the
      :ref:`tool-lps2lts` output.
        
.. note::

   In all the verification above, the command :ref:`tool-pbessolve` can be
   replaced by :ref:`tool-pbespgsolve`. The latter uses a different algorithm
   for doing the actual verification; it first translates the pbes into a
   parity game, and then solves the parity game.
   
.. note::

   For ``N >= 4``, the minimal number of calls that needs to be made is
   ``2N - 4``, as was shown in [Hur00]_. 
            
.. [Hur00] C.A.J. Hurkens. Spreading gossip efficiently.
   Nieuw Archief voor Wiskunde, Vol. 5/1, No. 2, p208-. 2000.   
   `(PDF) <http://alexandria.tue.nl/repository/articles/623317.pdf>`_
   


