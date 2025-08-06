// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/dparser.cpp
/// \brief add your file description here.

#include "mcrl2/core/dparser.h"
#include "mcrl2/core/detail/dparser_functions.h"
#include "mcrl2/utilities/logger.h"

#include <d.h>

#include <iostream>
#include <iomanip>
#include <locale>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else  
  import utilities;
#endif

// This corresponds to the D_ParseNode_to_PNode macro in parse.c.
PNode* D_ParseNode_to_PNode(D_ParseNode* _apn)
{
  if (_apn == nullptr)
  {
    return nullptr;
  }

  return reinterpret_cast<PNode*>(reinterpret_cast<char*>(_apn) - (sizeof(PNode) - sizeof(D_ParseNode)));
}


extern "C"
{
  extern D_ParserTables parser_tables_mcrl2;
}



namespace mcrl2::core {

std::string parse_node::add_context(const std::string& message) const
{
  return detail::add_context(&node->start_loc, message);
}

int parse_node::symbol() const
{
  return node->symbol;
}

int parse_node::child_count() const
{
  return d_get_number_of_children(node);
}

int parse_node::priority() const
{
  return D_ParseNode_to_PNode(node)->priority;
}

// 0 <= i < child_count()
parse_node parse_node::child(int i) const
{
  return parse_node(d_get_child(node, i));
}

parse_node parse_node::find_in_tree(int symbol) const
{
  return parse_node(d_find_in_tree(node, symbol));
}

std::string parse_node::string() const
{
  return std::string(node->start_loc.s, node->end - node->start_loc.s);
}

std::string parse_node::tree() const
{
  if (child_count() < 2)
  {
    return this->string();
  }
  std::stringstream result;
  result << "(" << child(0).tree();
  for (int i = 1; i < child_count(); ++i)
  {
    result << " " << child(i).tree();
  }
  result << ")";
  return result.str();
}

int parse_node::column() const
{
  return node->start_loc.col;
}

int parse_node::line() const
{
  return node->start_loc.line;
}

std::string parse_node::pathname() const
{
  return std::string(node->start_loc.pathname);
}

parse_node::~parse_node()
{
  if (parser)
  {
    free_D_ParseNode(parser, node);
  }
}

// Prints a tree of
std::string parser_table::tree(const core::parse_node& node) const
{
  std::stringstream result;
  result << symbol_name(node) << "(";
  if (node.child_count() == 0)
  {
    result << '"' << node.string() << '"';
  }
  else
  {
    result << tree(node.child(0));
  }
  for (int i = 1; i < node.child_count(); ++i)
  {
    result << " " << tree(node.child(i));
  }
  result << ")";
  return result.str();
}

// Returns the number of symbols in the table
unsigned int parser_table::symbol_count() const
{
  return m_table.nsymbols;
}

// Returns the name of the i-th symbol
std::string parser_table::symbol_name(unsigned int i) const
{
  if (i >= m_table.nsymbols)
  {
    print();
    std::ostringstream out;
    out << "parser_table::symbol_name: index " << i << " out of bounds!";
    throw std::runtime_error(out.str());
  }
  const char* name = m_table.symbols[i].name;
  if (!name)
  {
    return "";
  }
  return std::string(name);
}

std::string parser_table::symbol_name(const parse_node& node) const
{
  return symbol_name(node.symbol());
}

// Returns the 'start symbol' of the i-th symbol
int parser_table::start_symbol(unsigned int i) const
{
  return m_table.symbols[i].start_symbol;
}

// Returns true if the i-th symbol is of type D_SYMBOL_NTERM
bool parser_table::is_term_symbol(unsigned int i) const
{
  return m_table.symbols[i].kind == D_SYMBOL_NTERM;
}

unsigned int parser_table::start_symbol_index(const std::string& name) const
{
  for (unsigned int i = 0; i < symbol_count(); i++)
  {
    if (is_term_symbol(i) && symbol_name(i) == name)
    {
      return start_symbol(i);
    }
  }
  throw mcrl2::runtime_error("unknown start symbol '" + name + "'");
  return 0;
}

void parser_table::print() const
{
  std::clog << "--------------------" << std::endl;
  std::clog << "-   symbol table   -" << std::endl;
  std::clog << "--------------------" << std::endl;
  for (unsigned int i = 0; i < symbol_count(); i++)
  {
    std::clog << std::setw(3) << i << " " << symbol_name(i) << std::endl;
  }
  std::clog << "--------------------" << std::endl;
}

parser::parser(D_ParserTables& tables, D_AmbiguityFn ambiguity_fn, D_SyntaxErrorFn syntax_error_fn, std::size_t max_error_message_count)
  : m_table(tables)
{
  detail::set_dparser_max_error_message_count(max_error_message_count);
  m_parser = new_D_Parser(&tables, 0);
  m_parser->initial_globals = this;
  m_parser->save_parse_tree = 1;
  m_parser->initial_scope = nullptr;
  m_parser->dont_compare_stacks = 1;
  m_parser->dont_use_deep_priorities_for_disambiguation = 1;
  m_parser->dont_use_greediness_for_disambiguation = 1;
  m_parser->dont_use_height_for_disambiguation = 1;
  if (ambiguity_fn)
  {
    m_parser->ambiguity_fn = ambiguity_fn;
  }
  if (syntax_error_fn)
  {
    m_parser->syntax_error_fn = syntax_error_fn;
  }
}

parser::~parser()
{
  free_D_Parser(m_parser);
}

const parser_table& parser::symbol_table() const
{
  return m_table;
}

unsigned int parser::start_symbol_index(const std::string& name) const
{
  return m_table.start_symbol_index(name);
}

parse_node parser::parse(const std::string& text, unsigned int start_symbol_index, bool partial_parses)
{
  detail::reset_dparser_error_message_count();
  m_parser->start_state = start_symbol_index;
  m_parser->partial_parses = partial_parses ? 1 : 0;
  D_ParseNode* result = dparse(m_parser, const_cast<char*>(text.c_str()), static_cast<int>(text.size()));
  if (!result || m_parser->syntax_errors)
  {    
    if (result != nullptr) {
      free_D_ParseNode(m_parser, result);
    }
    throw mcrl2::runtime_error("syntax error");
  }
  return parse_node(result, m_parser);
}

void parser::print_symbol_table() const
{
  m_table.print();
}

std::string parser::indent(unsigned int count) const
{
  return std::string(count, ' ');
}

std::string parser::truncate(const std::string& s, unsigned int max_size) const
{
  std::string result = s.substr(0, max_size);

  // truncate at newline
  std::string::size_type pos = result.find('\n');
  if (pos != std::string::npos)
  {
    result = result.substr(0, pos);
  }

  return result;
}

void parser::print_tree(const parse_node& node, unsigned int level) const
{
  if (node)
  {
    std::string symbol = m_table.symbol_name(node.symbol());
    std::string prefix = indent(2 * level);
    std::cout << prefix << "--- " << symbol << " \"" << truncate(node.string()) << "\"" << std::endl;
    for (int i = 0; i <= node.child_count(); i++)
    {
      print_tree(node.child(i), level + 1);
    }
  }
}

void parser::destroy_parse_node(const parse_node& node)
{
  free_D_ParseNode(m_parser, node.node);
}

/// \brief Callback function for nodes in the parse tree
void parser::announce(D_ParseNode& node_ref)
{
  parse_node node(&node_ref);
  std::cout << "parsed " << m_table.symbol_name(node.symbol()) << " " << node.string() << std::endl;
}

namespace detail {

std::string add_context(const d_loc_t* loc, const std::string& message)
{
  std::stringstream s;
  s << "Line " << loc->line << ", column " << loc->col << ": "
    << message << std::endl;
  char* beg = loc->s - loc->col;
  char* end = loc->s;
  while (*end != '\0' && *end != '\n' && *end != '\r')
  {
    ++end;
  }
  std::string line(beg, end);
  s << "  " << line << std::endl;
  for (int i = 0; i < loc->col + 2; ++i)
  {
    s << ' ';
  }
  s << '^';
  return s.str();
}

inline
bool is_all_of_type(D_ParseNode* nodes[], int n, const char* type, const core::parser_table& table)
{
  for (int i = 0; i < n; i++)
  {
    core::parse_node node(nodes[i]);
    if (table.symbol_name(node) != type)
    {
      return false;
    }
  }
  return true;
}

inline
void print_ambiguous_nodes(D_ParseNode* nodes[], int n, const char* type, const core::parser_table& table)
{
  mCRL2log(log::verbose) << "--- " << type << " ambiguity" << std::endl;
  for (int i = 0; i < n; ++i)
  {
    core::parse_node vi(nodes[i]);
    // mCRL2log(log::verbose) << vi.tree() << " " << table.tree(vi) << std::endl;
    mCRL2log(log::verbose) << "ALT " << table.tree(vi) << std::endl;
  }
}

inline
void print_chosen_node(D_ParseNode* node, const core::parser_table& table)
{
  core::parse_node vi(node);
  mCRL2log(log::verbose) << "CHOOSE " << table.tree(vi) << std::endl;
}

/// \brief Function for resolving parser ambiguities.
D_ParseNode* ambiguity_fn(struct D_Parser * /*p*/, int n, struct D_ParseNode **v)
{
  core::parser_table table(parser_tables_mcrl2);

  // Maurice:
  //
  // Our grammar is ambiguous for RegFrm and ActFrm with brackets, an RegFrm expression (a && b) can be parsed as both RegFrm(ActFrm) and ActFrm(ActFrm) using either rule
  // RegFrm :: '(' RegFrm ')' | ActFrm or rule ActFrm :: '(' ActFrm ')' | Id.
  //
  // Introducing a ActFrmNoBrackets doesn't work since rule priorities are only applied to recursive applications of the rule, i.e., T :: T && T, and not T :: U && U, and U :: T.
  //
  // Instead we take the parse tree that starts with RegFrm(ActFrm).  
  for (int i = 1; i < n; ++i)
  {
    core::parse_node vi(v[i]);
    if (table.symbol_name(vi) == "RegFrm" && table.symbol_name(vi.child(0)) == "ActFrm")
    {
      mCRL2log(log::trace) << "Prioritising the bracket rule for RegFrm." << std::endl;
      return vi.node;
    }
  }

  // Maurice: We should take the parse tree with the lowest priority at the root and not the highest priority.
  // in the expression a * b + c, the parse tree (a * b) + c has "+" as root, and not "*", which would be a * (b + c), whereas * has a higher priority than +.
  // There are at least two nodes.
  core::parse_node candidate(v[0]);
  mCRL2log(log::trace) << "First candidate " << candidate.tree() << " with priority " << candidate.priority() << std::endl;

  // Indicates that that the candidate is actually lower priority than another node.
  bool chosen_candidate = true;
  for (int i = 1; i < n; ++i)
  {
    core::parse_node vi(v[i]);
    if (table.symbol_name(vi) == table.symbol_name(candidate))
    {    
      if (vi.priority() == candidate.priority()) 
      {
        // There are two nodes with the same priority.
        mCRL2log(log::trace) << "Two candidates with the same priority." << std::endl;
        chosen_candidate = false;
      }

      if (vi.priority() < candidate.priority()) 
      {
        mCRL2log(log::trace) << "Selecting " << vi.tree() << " as the parse tree with lower priority " << vi.priority() <<  std::endl;
        candidate = vi;
      }
    }
  }

  if (chosen_candidate) {
    // The candidate has a lower priority than the other nodes.
    mCRL2log(log::trace) << "The parse tree with the lowest priority: " << candidate.tree() << std::endl;
    return candidate.node;
  }

  // If we reach this point, then the ambiguity is unresolved. We print all
  // ambiguities on the debug output, then throw an exception.
  for (int i = 0; i < n; ++i)
  {
    core::parse_node vi(v[i]);
    mCRL2log(log::info) << "Candidate: " << vi.tree() << std::endl;
    mCRL2log(log::debug) << "Candidate (verbose): " << table.tree(vi) << std::endl;
  }
  
  throw mcrl2::runtime_error("Failed to parse the input. Several ambiguous parse trees where found, see the candidates above.");
}

static void log_location(struct D_Parser *ap)
{
  // We recover information about the last parsed node by casting D_Parser to Parser, which
  // is the structure that the dparser library internally uses to keep its administration in.
  std::string after;
  SNode *s = ((Parser*)ap)->snode_hash.last_all;
  ZNode *z = s != nullptr ? s->zns.v[0] : nullptr;
  while (z != nullptr && z->pn->parse_node.start_loc.s == z->pn->parse_node.end)
  {
    z = (z->sns.v && z->sns.v[0]->zns.v) ? z->sns.v[0]->zns.v[0] : nullptr;
  }
  if (z && z->pn->parse_node.start_loc.s != z->pn->parse_node.end)
  {
    after = std::string(z->pn->parse_node.start_loc.s, z->pn->parse_node.end);
  }

  std::string message = "syntax error";
  if (!after.empty())
  {
    message = message + " after '" + after + "'";
  }
  mCRL2log(log::error) << add_context(&ap->loc, message) << std::endl;
}

void syntax_error_fn(struct D_Parser *ap)
{
  core::detail::increment_dparser_error_message_count();
  if (core::detail::get_dparser_error_message_count() > core::detail::get_dparser_max_error_message_count())
  {
    return;
  }
  log_location(ap);
  if (ap->loc.s == nullptr)
  {
    mCRL2log(log::error) << "Unexpected end of input." << std::endl;
  }
  else
  {
    // Dive into the internals of dparser to recover some extra diagnostics.
    Parser* p = (Parser*)ap;
    if (p->pnode_hash.all && p->pnode_hash.all->latest)
    {
      core::parse_node n(&p->pnode_hash.all->latest->parse_node);
      D_Symbol &s = p->t->symbols[n.symbol()];
      if (s.kind == D_SYMBOL_INTERNAL)
      {
        /* DParser stores production rules in order: search for the corresponding nonterminal. */
        int parentsym = n.symbol() - 1;
        while (p->t->symbols[parentsym].kind == D_SYMBOL_INTERNAL)
        {
          --parentsym;
        }
        s = p->t->symbols[parentsym];
      }

      switch (s.kind)
      {
      case D_SYMBOL_STRING:
      case D_SYMBOL_TOKEN:
        {
          std::locale loc;
          mCRL2log(log::error) << "Unexpected "
                                         << (std::isalpha(n.string()[0], loc) ? "keyword " : "")
                                         << "'" << n.string() << "'" << std::endl;
        }
        break;
      case D_SYMBOL_NTERM:
        mCRL2log(log::error) << "Unexpected " << s.name << " '" << n.string() << "'" << std::endl;
        break;
      default:
        // TODO: check if we can give more sensible output in the remaining cases.
        break;
      }
    }
  }
}

} // namespace detail

void parser::custom_parse_error(const std::string& message) const
{
  core::detail::increment_dparser_error_message_count();
  if (core::detail::get_dparser_error_message_count() > core::detail::get_dparser_max_error_message_count())
  {
    return;
  }
  detail::log_location(m_parser);
  mCRL2log(log::error) << message << std::endl;
}

} // namespace mcrl2::core



