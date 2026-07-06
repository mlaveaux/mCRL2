Parsing
=======

The Core library provides C++ wrappers around the
`dparser <https://github.com/graydon/dparser>`_ GLR parser library. Three
grammars are compiled into the library and exposed as ``extern "C"`` symbols:

.. code-block:: c++

   extern D_ParserTables parser_tables_mcrl2;  // full mCRL2 language
   extern D_ParserTables parser_tables_fsm;    // FSM file format
   extern D_ParserTables parser_tables_dot;    // DOT graph format

The grammar source for the mCRL2 language is
``libraries/core/source/mcrl2_syntax.g`` and is compiled to C by
``make_dparser`` during the build. For a description of the formal syntax see
:doc:`mcrl2_syntax_definition`.

Key classes
-----------

.. table:: Classes in ``mcrl2/core/dparser.h`` and ``mcrl2/core/parse.h``
   :widths: 30 20 50

   ========================================  ===============  ============================================================
   Class                                     Header           Purpose
   ========================================  ===============  ============================================================
   ``parse_node``                            ``dparser.h``    Wraps a ``D_ParseNode``. Gives access to the matched
                                                              symbol (``symbol()``), source text (``string()``),
                                                              children (``child(i)``, ``child_count()``), and
                                                              source location (``line()``, ``column()``).
   ``parser_table``                          ``dparser.h``    Wraps ``D_ParserTables``. Supports symbol lookup by index
                                                              (``symbol_name()``) and finding start symbols by name
                                                              (``start_symbol_index()``).
   ``parser``                                ``dparser.h``    Wraps ``D_Parser``. Drives a parse run via ``parse()``,
                                                              holds the ``parser_table``, and provides
                                                              ``destroy_parse_node()`` for cleanup.
   ``parser_actions``                        ``parse.h``      Base class for parse-tree traversal. Provides
                                                              ``traverse(node, f)`` that walks the tree calling ``f``
                                                              on each node, and factory functions for ``visitor`` and
                                                              ``collector`` callbacks (see below).
   ``default_parser_actions``                ``parse.h``      Extends ``parser_actions`` with ``parse_list()``,
                                                              ``parse_vector()``, ``parse_Id()``, and
                                                              ``parse_Number()`` helpers.
   ``parse_node_exception``                  ``parse.h``      Exception thrown when a parse node does not have the
                                                              expected structure. Includes line/column context.
   ``parse_node_unexpected_exception``       ``parse.h``      Specialisation of ``parse_node_exception`` for nodes
                                                              of an entirely unexpected type; adds a dump of the
                                                              node to the error message.
   ========================================  ===============  ============================================================

Parsing an expression
---------------------

The following pattern is used throughout the mCRL2 libraries to parse a
single expression or specification:

.. code-block:: c++

   #include "mcrl2/core/parse.h"
   #include "mcrl2/core/detail/dparser_functions.h"

   // 1. Create a parser bound to a grammar and the standard callbacks.
   core::parser p(parser_tables_mcrl2,
                  core::detail::ambiguity_fn,
                  core::detail::syntax_error_fn);

   // 2. Look up the start symbol for the non-terminal you want to parse.
   unsigned int start = p.start_symbol_index("DataExpr");

   // 3. Parse a string. The caller must destroy the returned node.
   core::parse_node node = p.parse(text, start, /*partial_parses=*/false);

   // 4. Walk the tree using parser_actions.
   core::default_parser_actions actions(p);
   MyResult result = actions.parse_MyNonTerminal(node);

   // 5. Release the parse tree.
   p.destroy_parse_node(node);

The ``parse_identifier`` function in ``parse.h`` is the simplest
self-contained example of this pattern.

Parse-tree traversal helpers
-----------------------------

``parser_actions`` provides two inner callback types for use with
``traverse()``:

``visitor``
   Calls a function ``f(node)`` whenever a node's symbol name matches a given
   type string. Traversal stops descending into a matched node.

``collector``
   Like ``visitor``, but appends ``f(node)`` to a container instead of calling
   a side-effecting function.

``set_collector``
   Like ``collector``, but inserts into a set container.

Factory functions ``make_visitor``, ``make_collector``, and
``make_set_collector`` are provided so that the callback type can be deduced
automatically.

Error handling
--------------

Two dparser callback functions in ``mcrl2/core/detail/dparser_functions.h``
are registered with the parser at construction time:

``ambiguity_fn``
   Resolves syntactic ambiguities in the ``_ -> _ <> _`` (if-then-else)
   operator for process expressions by preferring the rule with higher
   priority.

``syntax_error_fn``
   A custom syntax-error handler that prints both line number and column to
   give more helpful diagnostics than dparser's default output.

The error counter in ``mcrl2/core/dparser.h`` (``detail::dparser_error_message_count``)
is used internally to cap the number of error messages emitted during a single
parse run. The relevant free functions are:

.. code-block:: c++

   void  reset_dparser_error_message_count();
   void  increment_dparser_error_message_count();
   std::size_t get_dparser_error_message_count();
   void  set_dparser_max_error_message_count(std::size_t n);  // default: 1

Parser utility helpers
-----------------------

``mcrl2/core/parser_utility.h`` provides predicates and warning functions
for detecting operator-precedence patterns that changed between mCRL2
versions:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Helper
     - Purpose
   * - ``detail::is_and_or_node``
     - Detects ``x && (y || z)`` in the parse tree.
   * - ``detail::is_left_merge_merge``
     - Detects ``x ||_ (y || z)`` in the parse tree.
   * - ``warn_and_or(node)``
     - Prints a warning for each ``x && y || z`` occurrence (currently a no-op stub).
   * - ``warn_left_merge_merge(node)``
     - Prints a warning for each ``x ||_ y || z`` occurrence (currently a no-op stub).
   * - ``detail::foreach_parse_node(x, f)``
     - Calls ``f`` on every node in the subtree rooted at ``x``.
