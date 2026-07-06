Introduction
============
The Core library provides low-level interfaces to the internally used data
structures. Using this library, internal data structures can be manipulated
directly, parsed from and printed to human-readable textual descriptions.

Concepts
========
The low-level data structures are represented as ATerms (short for Annotated
Terms) [ATerm]_. ATerms represent tree-like data structures with features such
as maximal subterm sharing, automatic garbage collection, and compressed binary
exchange formats. The ATerms implementation used by mCRL2 is provided by the
``atermpp`` library.

Structure
=========
The Core library is a collection of headers and source files in the
``mcrl2::core`` namespace. Its main components are:

- **Identifier strings** (``identifier_string.h``) — a thin wrapper around
  ``atermpp::aterm_string`` used to represent identifiers throughout the mCRL2
  libraries.

- **Parsing infrastructure** (``dparser.h``, ``parse.h``) — C++ wrappers around
  the `dparser <https://github.com/graydon/dparser>`_ GLR parser library. Provides
  ``parse_node``, ``parser_table``, ``parser``, and ``parser_actions`` for
  traversing and collecting information from parse trees. See
  :doc:`core_parsing` for details.

- **Traversal and builder framework** (``traverser.h``, ``builder.h``,
  ``add_binding.h``) — CRTP base classes for read-only traversals and
  term-transforming builders over mCRL2 data structures. Each library in mCRL2
  extends these bases with overloads for its own types. The common patterns are
  described in :doc:`common_functionality`.

- **Pretty printing** (``print.h``, ``print_format.h``) — a ``stream_printer``
  base class and ``pp()`` helper function for converting internal ATerm
  representations to human-readable text.

- **ATerm loading** (``load_aterm.h``) — utilities for reading ATerms from
  binary or text streams.

- **Soundness checks** (``detail/soundness_checks.h``) — auto-generated
  predicate templates used to validate term structure at runtime (enabled when
  ``MCRL2_ENABLE_SOUNDNESS_CHECKS`` is set).

References
----------

.. [ATerm] M.G.J. van den Brand and P. Klint (2007). "ATerms for manipulation and
   exchange of structured data: It's all about sharing." *Information and Software
   Technology* 49:55--64.
