A Vending Machine
-----------------

**Contribution of this section**

#. specifying processes,
#. linearisation,
#. state space exploration,
#. visualisation of state spaces,
#. comparison/reduction using behavioural equivalences, and
#. verification of modal mu-calculus properties.

**New tools**:
:ref:`tool-mcrl22lps`, :ref:`tool-lps2lts`, :ref:`tool-ltsgraph`, :ref:`tool-ltscompare`,
:ref:`tool-ltsconvert`, :ref:`tool-lps2pbes`, :ref:`tool-pbessolve`.

Our first little step consists of number of variations on the good old
vending machine, a user ``User`` interacting with a machine ``Mach``. By
way of this example we will encounter the basic ingredients of
mCRL2. In the first variation of the vending machine, a very
primitive machine and user, are specified. Some properties are
verified. In the second variation non-determinism is considered and,
additionally, some visualization and comparison tools from the toolset
are illustrated. The third variation comes closer to a rudimentary
prototype specification.

First variation
^^^^^^^^^^^^^^^
After inserting a coin of 10 cents, the user can push the button for
an apple. An apple will then be put in the drawer of the machine. See
Figure :ref:`fig-mach1`.

.. _fig-mach1:

.. figure:: /_static/tutorial/machine/mach1.*
   :align: center

   Vending machine 1

Vending machine 1 can be specified by the following mCRL2, also
included in the file :download:`vm01.mcrl2 <files/vm01.mcrl2>`.

.. literalinclude:: files/vm01.mcrl2
   :language: mcrl2

The specification is split in three sections:

#. ``act``, a declaration of actions of 6 actions,
#. ``proc``, the definition of 2 processes, and
#. ``init``, the initialization of the system.

The process ``User`` is recursively defined as doing an ``ins10``
action, followed by an ``optA`` action, followed by the process ``User``
again. The process ``Mach`` is similar, looping on the action ``acc10``
followed by the action ``putA``. Note, only four actions are used in
the definition of the processes. In particular, the action ``coin`` and
``ready`` are not referred to.

The initialization of the system has a typical form. A number of
parallel processes, in the context of a communication function, with a
limited set of actions allowed. So, ``||`` is the parallel
operator, in this case putting the processes ``User`` and ``Mach`` in
parallel. The communication function is the first argument of the
``comm`` operator. Here, we have that synchronization of an
``ins10`` action and an ``acc10`` action yields the action ``coin``,
whereas synchronization of ``optA`` and ``putA`` yields ``ready``. The
actions of the system that are allowed, are mentioned in the first
argument of the allow operator ``allow``. Thus, for our first
system only ``coin`` and ``ready`` are allowed actions.

We compile the specification in the file :download:`vm01.mcrl2 <files/vm01.mcrl2>`
to a so-called linear process, saved in the file :file:`vm01.lps`. This
can be achieved by running::

  $ mcrl22lps vm01.mcrl2 vm01.lps

on the command line. The linear process is the internal representation
format of mCRL2, and is not meant for human inspection. However, from
:file:`vm01.lps` a labeled transition system, LTS for short, can be
obtained by running::

  $ lps2lts vm01.lps vm01.lts

which can be viewed by the :ref:`tool-ltsgraph` facility, by typing::

   $ ltsgraph vm01.lts

at the prompt.  Some manual beautifying yields the picture in
Figure :ref:`fig-lts-vm01`.

.. _fig-lts-vm01:

.. figure:: /_static/tutorial/machine/vm01.*
   :align: center

   LTS of vending machine 1

Apparently, starting from state 0 the system shuttles between state
0 and 1 alternating the actions ``coin`` and ``ready``. Enforced by
the allow operator, unmatched ``ins10``, ``acc10``,
``optA`` and ``putA`` actions are excluded. The actions synchronize
pairwise, ``ins10`` with ``acc10``, ``optA`` with ``putA``, to produce
``coin`` and ``ready``, respectively.

As a first illustration of model checking in mCRL2, we consider some
simple properties to be checked against the
specification :download:`vm01.mcrl2 <files/vm01.mcrl2>`. Given the LTS of the
system, the properties obviously hold.

a. :download:`vm01a.mcf <files/vm01a.mcf>`:

   .. literalinclude:: files/vm01a.mcf
      :language: mcrl2

   In this property, ``[true*]`` represents all finite sequences of actions
   starting from the initial state. ``<true*.ready>`` expresses the existence
   of a sequence of actions ending with the action ``ready``. The last
   occurence of ``true`` in this property is a logical formula to be evaluated
   in the current state. Thus, if the property is satisfied by the system,
   then after any finite sequence of actions, ``[true*]``, the system can
   continue with some finite sequence of actions ending with ``ready``,
   ``<true*.ready>``, and reaches a state in which the formula ``true`` holds.
   Since ``true`` always holds, this property states that a next ``ready`` is
   always possible.

b. :download:`vm01b.mcf <files/vm01b.mcf>`:

   .. literalinclude:: files/vm01b.mcf
      :language: mcrl2

   This property is less liberal than property (a). Here, ``<ready> true``
   requires a ``ready`` action to be possible for the system, after any finite
   sequence, ``[true*]``. This property does not hold. A ``ready`` action is
   not immediately followed by a ``ready`` action again. Also, ``ready`` is
   not possible in the initial state.

c. :download:`vm01c.mcf <files/vm01c.mcf>`:

   .. literalinclude:: files/vm01c.mcf
      :language: mcrl2

   This property uses the complement construct. ``!coin`` are all actions
   different from ``coin``. So, any sequence of actions with ``ready`` as its
   one but final action and ending with an action different from ``coin``,
   leads to a state where ``false`` holds. Since no such state exists, there
   are no path of the form ``true*.ready.!coin``. Thus, after any ``ready``
   action, any action that follows, if any, will be ``coin``.

d. :download:`vm01d.mcf <files/vm01d.mcf>`:

   .. literalinclude:: files/vm01d.mcf
      :language: mcrl2

   This property is a further variation involving conjunction ``&&``.

Model checking with mCRL2 is done by constructing a so-called
parameterised boolean equation system or PBES from a linear process
specification and a modal :math:`\mu`-calculus formula. For example, to
verify property (a) above, we call the :ref:`tool-lps2pbes`
tool. Assuming property (a) to be in file :download:`vm01a.mcf
<files/vm01a.mcf>`, running::

   $ lps2pbes vm01.lps -f vm01a.mcf vm01a.pbes

creates from the system in linear format and the formula in the file
:download:`vm01.mcrl2 <files/vm01.mcrl2>` right after the ``-f`` switch, a PBES
in the file :file:`vm01a.pbes`. On calling the PBES solver on
:file:`vm01a.pbes`::

   $ pbessolve vm01a.pbes

the mCRL2 tool answers::

   true

So, for vending machine 1 it holds that action ``ready`` is always
possible in the future. Instead of making separate steps explicity,
the verification can also be captured by a single, pipe-line command::

   $ mcrl22lps vm01.mcrl2 | lps2pbes -f vm01a.mcf | pbessolve

Running the other properties yields the expected
results. Properties (c) and (d) do hold, property (b) does not hold, as
indicated by the following snippet::

   $ mcrl22lps vm01.mcrl2 | lps2pbes -f vm01b.mcf | pbessolve
   false
   $ mcrl22lps vm01.mcrl2 | lps2pbes -f vm01c.mcf | pbessolve
   true
   $ mcrl22lps vm01.mcrl2 | lps2pbes -f vm01d.mcf | pbessolve
   true

Second variation
^^^^^^^^^^^^^^^^
Next, we add a chocolate bar to the assortment of the vending machine. A
chocolate bar costs 20 cents, an apple 10 cents. The machine will now accept
coins of 10 and 20 cents. The scenarios allowed are (i) insertion of 10 cent and
purchasing an apple, (ii) insertion of 10 cent twice or 20 cent once and
purchasing a chocolate bar. Additionally, after insertion of money, the user can
push the change button, after which the inserted money is returned. See
Figure :ref:`fig-mach2`.

.. _fig-mach2:

.. figure:: /_static/tutorial/machine/mach2.*
   :align: center

   Vending machine 2

.. admonition:: Exercise

   Extend the following mCRL2 specification (:download:`vm02-holes.mcrl2 <files/vm02-holes.mcrl2>`)
   to describe the vending machine sketched above, and save the resulting
   specification as :file:`vm02.mcrl2`. The actions that are involved,
   and a possible specification of the ``Mach`` process have been given.
   The machine is required to perform a ``prod`` action for administration
   purposes.

   .. literalinclude:: files/vm02-holes.mcrl2
      :language: mcrl2

   Linearise your specification using :ref:`tool-mcrl22lps`, saving the LPS
   as :file:`vm02.lps`.

.. admonition:: Solution
   :class: collapse

   A sample solution is available in :download:`vm02.mcrl2 <files/vm02.mcrl2>`.
   This can be linearised using::

   $ mcrl22lps vm02.mcrl2 vm02.lps

A visualization of the specified system can be obtained by first
converting the linear process into a labeled transition system (in
so-called SVC-format) by::

  $  lps2lts vm02.lps vm02.lts

and next loading the SVC file :file:`vm02.lts` into the ltsgraph tool
by::

  $  ltsgraph vm02.lts

The LTS can be beautified (a bit) using the ``start`` button in the
optimization panel of the user interface. Manual manipulation by
dragging states is also possible. For small examples, increasing
the natural transition length may provide better results.

.. admonition:: Exercise

   Prove that your specification satisfies the following properties:

   a. no three 10ct coins can be inserted in a row,
   b. no chocolate after 10ct only, and
   c. an apple only after 10ct, a chocolate after 20ct.

.. admonition:: Solution
   :class: collapse

   Each of the properties can be expressed as a µ-calculus formula. Possible
   solutions are given as :download:`vm02a.mcf <files/vm02a.mcf>`,
   :download:`vm02b.mcf <files/vm02b.mcf>`, and
   :download:`vm02c.mcf <files/vm02c.mcf>`.

   Each of the properties can be checked using a combination of :ref:`tool-mcrl22lps`,
   :ref:`tool-lps2pbes` and :ref:`tool-pbessolve`. The following is a sample script that
   performs the verification::

     $ mcrl22lps vm02.mcrl2 vm02.lps
     $ lps2pbes vm02.lps -f vm02a.mcf | pbessolve
     true
     $ lps2pbes vm02.lps -f vm02b.mcf | pbessolve
     true
     $ lps2pbes vm02.lps -f vm02c.mcf | pbessolve
     true

   So the conclusion of the verification is that all three properties hold.

The file :download:`vm02-taus.mcrl2 <files/vm02-taus.mcrl2>` contains the
specification of a system performing ``coin10`` and ``coin20`` actions as well
as so-called :math:`\tau`-steps. Using the :ref:`tool-ltscompare` tool you can
compare your model under branching bisimilarity with the LTS of the system
:file:`vm02-taus`, after hiding the actions ``readyA``, ``readyC``, ``out10``,
``out20``, ``prod`` using the following command::

  $ ltscompare -ebranching-bisim --tau=out10,out20,readyA,readyC,prod vm02.lts vm02-taus.lts

.. note::

  You first need to generate the state space of :download:`vm02-taus.mcrl2
  <files/vm02-taus.mcrl2>` using :ref:`tool-mcrl22lps` and :ref:`tool-lps2lts`.

Using :ref:`tool-ltsconvert`, the LTS for :file:`vm02.mcrl2` can be minimized
with respect to branching bisimulation after hiding the readies and returns::

  $ ltsconvert -ebranching-bisim --tau=out10,out20,readyA,readyC,prod vm02.lts vm02min.lts

.. admonition:: Exercise

   Compare the LTSs :file:`vm02min.lts` and `vm02-taus.lts` visually using
   :ref:`tool-ltsgraph`.

Third variation
^^^^^^^^^^^^^^^
A basic version of a vending machine with parametrized actions is available
in the file :download:`vm03-basic.mcrl2 <files/vm03-basic.mcrl2>`.

.. admonition:: Exercise

   Modify this specification such that all coins of denomination
   50ct, 20ct, 10ct and 5ct can be inserted. The machine accumulates
   upto a total of 60 cents. If sufficient credit, an apple or
   chocolate bar is supplied after selection. Money is returned after
   pressing the change button.
