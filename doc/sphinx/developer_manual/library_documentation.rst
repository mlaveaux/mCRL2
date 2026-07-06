mCRL2 library documentation
###########################

ATermpp
=========

.. toctree::
   :maxdepth: 1

   libraries/atermpp/aterm_library
   libraries/atermpp/aterm_implementation_notes

Core
======

.. toctree::
   :maxdepth: 1

   libraries/core/core_library
   libraries/core/core_parsing
   libraries/core/common_functionality
   libraries/core/mcrl2_syntax_definition

Data
======

.. toctree::
   :maxdepth: 2

   libraries/data/data_library



GUI
=====

.. toctree::
   :maxdepth: 1

   libraries/gui/gui_library


LPS
======

.. toctree::
   :maxdepth: 1

   libraries/lps/lps_library
   libraries/lps/lps_implementation_notes
   libraries/lps/confluence_detection
   libraries/lps/state_space_exploration
   libraries/lps/communication_open_terms


LTS
=====

.. toctree::
   :maxdepth: 1

   libraries/lts/lts_library




Modal formula
================

.. toctree::
   :maxdepth: 1

   libraries/modal_formula/modal_formulas_library

PBES
======

.. toctree::
   :maxdepth: 1

   libraries/pbes/pbes_library
   libraries/pbes/pbesinst
   libraries/pbes/pbes-inst-solving
   libraries/pbes/pbes_rewriters
   libraries/pbes/pbesabsinthe
   libraries/pbes/pbes_generation
   libraries/pbes/pbes_greybox_notes
   libraries/pbes/pbes_implementation_notes

Process
=========

.. toctree::
   :maxdepth: 1

   libraries/process/process_library
   libraries/process/global_variables
   libraries/process/process_implementation_notes

SMT
======

.. toctree::
   :maxdepth: 1

   libraries/smt/smt_library

Utilities
===========

.. toctree::
   :maxdepth: 1

   libraries/utilities/utilities_library

Symbolic
===========

.. toctree::
   :maxdepth: 1

   libraries/symbolic/symbolic_library

Source code reference
=====================

We generate the complete source code reference using Doxygen, which can be viewed `here <../doxygen/index.html>`_.

Code generation
===============
A substantial part of the mCRL2 code is generated using scripts. This section gives an overview of the scripts and their purpose.

* ``build/code_generation/generate_classes.py``               Generates classes that derive from atermpp::aterm_appl. The class specifications that are used as input are located in ``build/code_generation/mcrl2_classes.py``.

* ``build/code_generation/generate_data_types.py``            Generates code for standard data types like bag, bool, int, etc.

* ``build/code_generation/generate_template_overloads.py``    Generates overloads for several template functions (pp, find, normalize_sorts, etc.). These overloads are added to header files to avoid the inclusion of heavy templated header files that are expensive to compile. The implementations of the functions are added to files in the source directory, while the declarations need to be added manually (since it is not always known where to put them).

* ``build/code_generation/generate_term_functions.py``        Generates low level code for creating terms (soundness checks, constructors). The script requires command line arguments (see the help function for details). The description of the internal format in the file ``doc/specs/mcrl2.internal.txt`` is used as input for the script.

* ``build/code_generation/generate_traverser_functions.py``   Generates functions that are implemented using the traverser framework (find/replace/rewrite etc).

* ``build/code_generation/generate_traversers.py``            Generates the contents of traverser.h and builder.h files. The class specifications that are used as input are located in ``build/code_generation/mcrl2_classes.py``.

.. note::

   Most of the scripts print a report of the files that have been updated, and the files that have stayed the same.
