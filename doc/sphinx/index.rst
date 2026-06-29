.. include:: /_includes/toplevel.inc

.. title:: Home

.. _release: https://www.mcrl2.org/release/user_manual/download.html
.. _version: download.html

.. release::

   This website was last updated on: |today|.

#################
The mCRL2 toolset
#################

.. rst-class:: image-no-text-wrap
.. image:: /_static/book.*
          :align: left
          :alt: Modeling and analysis of communicating systems (book)
          :target: https://mitpress.mit.edu/9780262027717/

mCRL2 is a formal specification language with an associated toolset. The
toolset can be used for modelling, validation and verification of
concurrent systems and protocols. It can be run on Windows, Linux, Apple
MacOS and FreeBSD.

The toolset supports a collection of tools for linearisation, simulation,
state-space exploration and generation and tools to optimise and analyse
specifications. Moreover, state spaces can be manipulated, visualised and
analysed.

The mCRL2 toolset is developed at the department of Mathematics and
Computer Science of the `Technische Universiteit Eindhoven
<http://www.tue.nl>`_, in collaboration with
the `University of Twente <http://fmt.cs.utwente.nl>`_.

.. note::

   The latest release, mCRL2 202507.0, is now available from the :ref:`download` page.
   See the `changelog <https://github.com/mCRL2org/mCRL2/blob/master/CHANGES>`_ for a detailed list
   of changes.

.. note::

   The Coursera lectures are now freely available on Youtube, divided in four parts. It is intended
   that the full course including exercises will soon be available on openlearnity.com.

   - `Automata and behavioural equivalences <https://www.youtube.com/watch?v=1brp80vNSy4&list=PLTlBz41EgXBTRm4jcsRwKWoTazWxiieNi>`__.
   - `Model process behaviour <https://www.youtube.com/watch?v=1yP9kMpkv38&list=PLTlBz41EgXBTmbFWpxMsZRT2vcnxxC6fF>`__.
   - `Requirements by modal formulas <https://www.youtube.com/watch?v=-vmhNqO16MA&list=PLTlBz41EgXBQ_vIWcFcGJa8_-a8-o19_F>`__.
   - `Modelling software, protocols, and other behaviour <https://www.youtube.com/watch?v=JV2djZV_dls&list=PLTlBz41EgXBSZPXTzP_MuvbO6MnmTOly6>`__.

.. note::

   Listen to the `Parallel Rhapsody <https://suno.com/song/e90327be-a900-4fc3-82b2-8365cc837cfc>`_
   (thanks to Frederik Gosewehr and Gert Veltink).

.. note::

   .. image:: /_templates/mcrl2/static/understanding.*
      :align: right
      :width: 90px
      :alt: Understanding Behaviour of Distributed Systems using mCRL2 (book)
      :target: https://link.springer.com/book/10.1007/978-3-031-23008-0

   A new introductory book appeared, describing the mCRL2 specification language, the
   modal mu-calculus and the extensive toolset around these. The book is full of examples, varying from
   games to distributed algorithms. It is shown how such games can be solved,
   and how the (in)correctness of such algorithms can be understood.

.. note::

   FORTE 2022 Best Artefact Award: Process Algebra Can Save Lives: Static Analysis of XACML Access
   Control Policies using mCRL2. By Arshad, Horne, Johansen, Owe and Willemse.

.. note::
  .. image:: /_static/tacas-ae.*
          :align: right
          :width: 90px
          :alt: TACAS artefact evaluation badge

  A tool paper explaining the latest developments in mCRL2 was published in TACAS2019.
  The new features include counterexamples, probabilistic processes and efficient algorithms for several types of behavioural relations.
  As part of the review process, mCRL2 was tested and accepted by the TACAS artefact evaluation committee.
  The paper is `available <https://doi.org/10.1007/978-3-030-17465-1_2>`_ under open access.

.. note::

   On Coursera four consecutive courses consisting of approximately 60 lectures provide an introduction into the theory and practical
   use of the mCRL2 toolset. It is based on the book 'Modeling and Analysis of Communicating
   Systems' shown above. Unfortunately, these are now behind a pay wall. The courses are:

   - `Automata and behavioural equivalences <https://www.coursera.org/learn/automata-system-validation>`_.
   - `Model process behaviour <https://www.coursera.org/learn/system-validation-behavior>`_.
   - `Requirements by modal formulas <https://www.coursera.org/learn/system-validation-modal-formulas>`_.
   - `Modelling software, protocols, and other behaviour <https://www.coursera.org/learn/system-validation-software-protocols>`_.

   The slides used to make these courses are available by contacting J.F.Groote@tue.nl. They can also
   be downloaded from `MIT Press <https://mitpress.mit.edu/9780262027717/>`_.

.. note::
   .. image:: /_static/Verum-2-Logo-Vert-Colour.*
          :align: right
          :width: 75px
          :alt: Verum is the industry leader in reliable software design.

   mCRL2 is used by the company `Verum <https://www.verum.com>`_ as their verification engine.
   Verum provides model based software development environments for the languages ASD and Dezyne that
   allow to program proven correct embedded software with much less effort than "classical" programming.
   As it stands Verum is the industry leader in reliable software design.

.. list-table::
   :class: teamcity-blurb

   * - .. image:: /_static/teamcity.*
          :width: 75
          :alt: TeamCity
          :target: https://www.jetbrains.com/teamcity/
     - The mCRL2 project uses `TeamCity <https://www.jetbrains.com/teamcity/>`_ for regression testing during
       development and building the daily snapshots.


.. toctree::
   :caption: Home
   :maxdepth: 2
   :hidden:

   /home/introduction
   /home/download
   /home/support
   /home/publications
   /home/showcases

.. toctree::
   :caption: Using mCRL2
   :maxdepth: 2
   :hidden:
   :glob:

   /user_manual/tutorial/tutorial
   /user_manual/language_reference/index
   /user_manual/tools/tools
   /user_manual/fundamentals/index

.. toctree::
   :caption: Development
   :maxdepth: 2
   :hidden:

   /developer_manual/contributing
   /developer_manual/build_instructions/instructions
   /developer_manual/guidelines
   /developer_manual/library_documentation

   /developer_manual/testing
   /developer_manual/performance

.. toctree::
   :caption: Maintainers
   :maxdepth: 2
   :hidden:

   /maintainer_manual/git_workflow
   /maintainer_manual/teamcity_testing
   /maintainer_manual/packaging
