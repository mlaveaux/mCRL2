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

- **Identifier strings** (:mcrl2:`identifier_string.h <identifier_string.h>`) — a thin wrapper around
  ``atermpp::aterm_string`` used to represent identifiers throughout the mCRL2
  libraries.

- **Parsing infrastructure** (:mcrl2:`dparser.h <dparser.h>`, :mcrl2:`parse.h <parse.h>`) — C++ wrappers around
  the `dparser <https://github.com/graydon/dparser>`_ GLR parser library. Provides
  :mcrl2:`parse_node <mcrl2::core::parse_node>`,
  :mcrl2:`parser_table <mcrl2::core::parser_table>`,
  :mcrl2:`parser <mcrl2::core::parser>`, and
  :mcrl2:`parser_actions <mcrl2::core::parser_actions>` for
  traversing and collecting information from parse trees. See
  :doc:`core_parsing` for details.

- **Traversal and builder framework** (:mcrl2:`traverser.h <traverser.h>`, :mcrl2:`builder.h <builder.h>`,
  :mcrl2:`add_binding.h <add_binding.h>`) — CRTP base classes for read-only traversals and
  term-transforming builders over mCRL2 data structures. The base classes are
  :mcrl2:`traverser <mcrl2::core::traverser>`,
  :mcrl2:`builder <mcrl2::core::builder>`, and
  :mcrl2:`add_binding <mcrl2::core::add_binding>`.
  Each library in mCRL2 extends these bases with overloads for its own types.
  The common patterns are described in :doc:`common_functionality`.

- **Pretty printing** (:mcrl2:`print.h <print.h>`, :mcrl2:`print_format.h <print_format.h>`) — a
  :mcrl2:`stream_printer <mcrl2::core::stream_printer>` base class and ``pp()``
  helper function for converting internal ATerm representations to human-readable
  text.

- **ATerm loading** (:mcrl2:`load_aterm.h <load_aterm.h>`) — utilities for reading ATerms from
  binary or text streams.

- **Soundness checks** (:mcrl2:`detail/soundness_checks.h <soundness_checks.h>`) — auto-generated
  predicate templates used to validate term structure at runtime (enabled when
  ``MCRL2_ENABLE_SOUNDNESS_CHECKS`` is set).

References
----------

.. [ATerm] M.G.J. van den Brand and P. Klint (2007). "ATerms for manipulation and
   exchange of structured data: It's all about sharing." *Information and Software
   Technology* 49:55--64.
