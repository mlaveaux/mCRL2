Symbolic library
================

The symbolic library contains the implementation of LDD helper functions and
classes used by the symbolic instantiation and solving algorithms. The symbolic
instantiation and solving algorithms themselves are located in the
:doc:`PBES library <../pbes/pbes_library>` and
:doc:`LPS library <../lps/lps_library>` respectively since these depend on the
PBES and LPS data
structures respectively. However, since these are all closely related we put the
documentation together here.


Symbolic Reachability
---------------------

See :doc:`symbolic-reachability`.

.. toctree::
	:maxdepth: 1

	list_decision_diagrams
	symbolic-reachability
	zielonka
