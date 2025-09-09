#!/usr/bin/env python3
# Author(s): Jeroen Keiren
# Copyright: see the accompanying file COPYING or copy at
# https://github.com/mCRL2org/mCRL2/blob/master/COPYING
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

# Classes for data expressions and sort expressions, used for code generation
# Basically reflects the parse tree of a specification, but allows for more
# flexibility.

import string
from string import Template
import copy
import re
import logging

# The following lists identifiers that must be
# escaped by a namespace when the corresponding function
# is called, to work around a bug in GCC 4.4
IDS_WITH_NAMESPACE = ['set_comprehension', 'bag_comprehension']

LOG = logging.getLogger("data")
# Uncomment to add lots of debug information
#LOG.setLevel(logging.DEBUG)

# create console handler and set level to debug
CH = logging.StreamHandler()
LOG.addHandler(CH)

def indent(n):
  result = ""
  i = 0
  while i < n:
    result += " "
    i += 1
  return result

def add_namespace(function_symbol_name, function_symbol_namespace, other_namespace = "undefined"):
  global IDS_WITH_NAMESPACE
  if (function_symbol_namespace == other_namespace or is_standard_function(function_symbol_name) or function_symbol_namespace == "undefined") and not (str(function_symbol_name) in IDS_WITH_NAMESPACE):
    return function_symbol_name
  else:
    return "sort_%s::%s" % (function_symbol_namespace, function_symbol_name)

# Remove trailing _ from a string
def remove_underscore(s):
  return s[:-1] if s.endswith('_') else s

# Escape an initial @ sign is present. Needed for doxygen code extraction
def escape(x):
  x = str(x)
  return '\\' + x if x.startswith('@') else x

def fcode(x, spec):
  if isinstance(x, str):
    return x
  else:
    return x.code(spec)

# Merge elements of a list of lists into a single list.
# Elements in the result are sorted and unique.
def merge(args):
  result = []
  for arg in args:
    for x in arg:
      if x not in result:
        result.append(x)
  return sorted(result)

# The elements of l1 that do not occur in l2
def difference(l1, l2):
  result = []
  for x in l1:
    if x not in l2:
      result.append(x)
  return result

def is_standard_function(s):
  return str(s) in ["equal_to", "not_equal_to", "if_", "less", "less_equal", "greater", "greater_equal"]

def target_sort(sort_expr):
  if isinstance(sort_expr, sort_arrow):
    return sort_expr.codomain
  else:
    return sort_expr

class identifier():
  def __init__(self, string):
    self.string = string

  def __str__(self):
    return self.string

  def __eq__(self, other):
    if hasattr(other, "string"):
      return self.string == other.string
    else:
      return self.string == other

  def __lt__(self, other):
    if hasattr(other, "string"):
      return self.string < other.string
    else:
      return self.string == other

  def __hash__(self):
    return hash(self.string)

class label():
  def __init__(self, id):
    assert(isinstance(id, identifier))
    self.identifier = id

  def __str__(self):
    return self.identifier.__str__()

  def __eq__(self, other):
    if hasattr(other, "identifier"):
      return self.identifier == other.identifier
    else:
      return NotImplemented

  def __hash__(self):
    return hash(self.identifier)

class variable_declaration():
  def __init__(self, id, sort_expr):
    assert(isinstance(id, identifier))
    assert(isinstance(sort_expr, sort_expression))
    self.id = id
    self.sort_expression = sort_expr

  def __str__(self):
    return "{0} : {1}".format(self.id, self.sort_expression)

  def __eq__(self, other):
    if hasattr(other, "id") and hasattr(other, "sort_expression"):
      return self.id == other.id and self.sort_expression == other.sort_expression
    else:
      return NotImplemented

  def sort_parameters(self, spec):
    return self.sort_expression.sort_parameters(spec)

  def code(self, spec):
    return "        variable v{0}(\"{0}\",{1})".format(self.id, self.sort_expression.code(spec))

class variable_declaration_list():
  def __init__(self, elements):
    assert(all([isinstance(x, variable_declaration) for x in elements]))
    self.elements = elements

  def __str__(self):
    return ", ".join([str(e) for e in self.elements])

  def __eq__(self, other):
    if hasattr(other, "elements"):
      return self.elements == other.elements
    else:
      return NotImplemented

  def push_back(self, element):
    if element not in self.elements:
      self.elements.append(element)

  def find_variable(self, variable):
    for e in self.elements:
      if e.id == variable.id:
        return e

  def code(self, spec):
    return "".join(["{0};\n".format(e.code(spec)) for e in self.elements])

class function_declaration():
  def __init__(self, id, sort_expr, l, internextern, definedby):
    assert(isinstance(id, identifier))
    assert(isinstance(sort_expr, sort_expression))
    assert(isinstance(l, label))
    assert(isinstance(internextern, internal_external))
    assert(isinstance(definedby, defined_by))
    self.id = id
    self.sort_expression = sort_expr
    self.label = l
    self.internextern = internextern
    self.definedby = definedby
    self.namespace = ""
    self.original_namespace = ""

  def __eq__(self, other):
    if (hasattr(other, "id") and hasattr(other, "label") and hasattr(other, "namespace") and
            hasattr(other, "original_namespace") and hasattr(other, "internextern") and hasattr(other, "definedby")):
      return (self.id == other.id and self.label == other.label and self.namespace == other.namespace and self.original_namespace == other.original_namespace and
              self.internextern == other.internextern and self.definedby == other.definedby)
    else:
      return NotImplemented

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string

  def projection_arguments(self, projection_arguments):
    if isinstance(self.sort_expression, sort_arrow):
      arguments = self.sort_expression.domain.projection_arguments()
      for (a, idx) in arguments:
        value = (self.id, self.label, idx)
        try:
          projection_arguments[a].append(value)
        except:
          projection_arguments[a] = [value]

    return projection_arguments

  def sort_parameters(self, spec):
    return self.sort_expression.sort_parameters(spec)

  def __str__(self):
    return "{0} : {1}".format(self.id, self.sort_expression)

  def generator_code(self, spec, function_declarations):
    if self.namespace != self.original_namespace or self.namespace != function_declarations.namespace or self.namespace != spec.namespace:
      return ""

    sort_params = self.sort_expression.sort_parameters(spec)
    # Determine whether this function is overloaded
    functions = [x for x in function_declarations.elements if (x.id == self.id) and (x.label == self.label) and (x.namespace == self.namespace)]
    extra_parameters = []
    if len(functions) > 1:
      try:
        extra_parameters.append(self.sort_expression.domain.code(spec))
      except:
        pass # in case sort_expression has no domain

    return "        result.push_back({0}({1}));\n".format(add_namespace(self.label, self.namespace), ", ".join([s.code(spec) for s in sort_params] + extra_parameters))

  def mCRL2_usable_functions(self, spec, function_declarations):
    if self.namespace != self.original_namespace or self.namespace != function_declarations.namespace or self.namespace != spec.namespace: # or self.internextern.is_internal(): Also functions declared internal are externally usable, as they can occur in lps's en pbes's that are sometimes typechecked, for instance when merging two specifications.
      return ""

    sort_params = self.sort_expression.sort_parameters(spec)
    # Determine whether this function is overloaded
    functions = [x for x in function_declarations.elements if (x.id == self.id) and (x.label == self.label) and (x.namespace == self.namespace)]
    extra_parameters = []
    if len(functions) > 1:
      try:
        extra_parameters.append(self.sort_expression.domain.code(spec))
      except:
        pass # in case sort_expression has no domain

    return "        result.push_back({0}({1}));\n".format(add_namespace(self.label, self.namespace), ", ".join([s.code(spec) for s in sort_params] + extra_parameters))

  def cplusplus_implementable_functions(self, spec, data_parameters, function_declarations):
    if self.namespace != self.original_namespace or self.namespace != function_declarations.namespace or self.namespace != spec.namespace or not self.definedby.is_defined_by_code():
      return ""

    sort_params = self.sort_expression.sort_parameters(spec)
    # Determine whether this function is overloaded
    functions = [x for x in function_declarations.elements if (x.id == self.id) and (x.label == self.label) and (x.namespace == self.namespace) and (x.definedby.is_defined_by_code())]
    extra_parameters = []
    if len(functions) > 1:
      try:
        extra_parameters.append(self.sort_expression.domain.code(spec))
      except:
        pass # in case sort_expression has no domain
    return "        result[{0}({1})]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>({0}_application,".format(add_namespace(self.label, self.namespace), ", ".join([s.code(spec) for s in sort_params] + extra_parameters), data_parameters) + \
                    "\"{0}_manual_implementation\");\n".format(add_namespace(self.label, self.namespace))


class internal_external():
  def __init__(self, status):
    assert(status=="internal" or status=="external")
    self.status = status

  def __eq__(self, other):
    assert(self.status=="internal" or self.status=="external")
    assert(other.status=="internal" or other.status=="external")
    return self.status == other.status

  def is_internal(self):
    assert(self.status=="internal" or self.status=="external")
    return self.status=="internal";

  def is_external(self):
    assert(self.status=="internal" or self.status=="external")
    return self.status=="external";

  def string(self):
    return self.status

class defined_by():
  def __init__(self, status):
    assert(status=="defined_by_rewrite_rules" or status=="defined_by_code")
    self.status = status

  def __eq__(self, other):
    assert(self.status=="defined_by_rewrite_rules" or self.status=="defined_by_code")
    assert(other.status=="defined_by_rewrite_rules" or other.status=="defined_by_code")
    return self.status == other.status

  def is_defined_by_code(self):
    assert(self.status=="defined_by_rewrite_rules" or self.status=="defined_by_code")
    return self.status == "defined_by_code"

  def is_defined_by_rewrite_rules(self):
    assert(self.status=="defined_by_rewrite_rules" or self.status=="defined_by_code")
    return self.status == "defined_by_rewrite_rules"

  def string(self):
    return self.status

class function_declaration_list():
  def __init__(self, elements):
    assert(all([isinstance(x, function_declaration) for x in elements]))
    self.elements = elements
    self.namespace = ""
    self.original_namespace = ""

  def __eq__(self, other):
    return self.elements == other.elements and self.namespace == other.namespace and self.original_namespace == other.original_namespace

  def has_cplusplus_implementable_code(self):
    for e in self.elements:
      if e.namespace == self.namespace and e.definedby.is_defined_by_code():
        return True
    return False

  def empty(self):
    return len(self.elements) == 0

  def push_back(self, element):
    assert(isinstance(element, function_declaration))
    self.elements += [element]

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    for e in self.elements:
      e.set_namespace(string)

  def sort_parameters(self, spec):
    return merge([e.sort_parameters(spec) for e in self.elements])

  def find_function(self, function, argumentcount):
    if function.id == "==":
      return function_declaration(function.id, sort_expression(), label(identifier("equal_to")),
                    internal_external("internal"), defined_by("defined_by_rewrite_rules"))
    elif function.id == "!=":
      return function_declaration(function.id, sort_expression(), label(identifier("not_equal_to")),
                    internal_external("internal"), defined_by("defined_by_rewrite_rules"))
    elif function.id == "if":
      return function_declaration(function.id, sort_expression(), label(identifier("if_")),
                    internal_external("internal"), defined_by("defined_by_rewrite_rules"))
    elif function.id == "<":
      return function_declaration(function.id, sort_expression(), label(identifier("less")),
                    internal_external("internal"), defined_by("defined_by_rewrite_rules"))
    elif function.id == ">":
      return function_declaration(function.id, sort_expression(), label(identifier("greater")),
                    internal_external("internal"), defined_by("defined_by_rewrite_rules"))
    elif function.id == "<=":
      return function_declaration(function.id, sort_expression(), label(identifier("less_equal")),
                    internal_external("internal"), defined_by("defined_by_rewrite_rules"))
    elif function.id == ">=":
      return function_declaration(function.id, sort_expression(), label(identifier("greater_equal")),
                    internal_external("internal"), defined_by("defined_by_rewrite_rules"))
    else:
      for e in self.elements:
        if e.id == function.id:
          if argumentcount == -1:
            return e
          elif argumentcount == 0 and not isinstance(e.sort_expression, sort_arrow):
            return e
          elif isinstance(e.sort_expression, sort_arrow) and len(e.sort_expression.domain.elements) == argumentcount:
            return e
    return None

  def merge_declarations(self, declarations, is_supertype):
    for d in declarations.elements:
      if is_supertype and any([(d.id == x.id) and (d.label == x.label) and (x.namespace == self.namespace) for x in self.elements]):
        d.set_namespace(self.namespace)
      self.elements.append(d)

  def __str__(self):
    return "".join(["{0}\n".format(e) for e in self.elements])

  def generator_code(self, spec, is_constructor_declaration = False):
    code = ""
    for e in self.elements:
      if not is_constructor_declaration or not spec.sort_specification.is_alias(target_sort(e.sort_expression)):
        code += e.generator_code(spec, self)
    return code

  def mCRL2_usable_functions(self, spec, is_constructor_declaration = False):
    code = ""
    for e in self.elements:
      if not is_constructor_declaration or not spec.sort_specification.is_alias(target_sort(e.sort_expression)):
        code += e.mCRL2_usable_functions(spec, self)
    return code

  def cplusplus_implementable_functions(self, spec, data_parameters, is_constructor_declaration = False):
    code = ""
    for e in self.elements:
      if not is_constructor_declaration or not spec.sort_specification.is_alias(target_sort(e.sort_expression)):
        code += e.cplusplus_implementable_functions(spec, data_parameters, self)
    return code

  def projection_code(self, spec):
    projection_arguments = {}
    for e in self.elements:
      if (e.namespace == e.original_namespace and e.namespace == spec.namespace):
        projection_arguments = e.projection_arguments(projection_arguments)

    code = []
    for p in projection_arguments:
      # Filter duplicates
      case_arguments = projection_arguments[p]
      unique_case_arguments = []
      for (id, label, idx) in case_arguments:
        if not (id, label, idx) in unique_case_arguments:
          unique_case_arguments.append((id, label, idx))

      assertions = ["is_{0}_application(e)".format(label) for (id, label, idx) in unique_case_arguments]

# First group by index
      index_table = {}
      for (id, label, idx) in unique_case_arguments:
        try:
          index_table[idx].append((id, label))
        except:
          index_table[idx] = [(id, label)]

      # If all projection project out the argument with the same index, do not generate a condition;
      # otherwise we generate an if statement for the different cases.
      projection = []

      if len(index_table) == 1:
        # The same position is projected in all cases, just take this position, no need to generate a condition.
        projection.append("return atermpp::down_cast<application>(e)[{0}];".format(list(index_table.keys())[0]))
      else:
        # A projection function of the same name is used to project different positions. We need to generate code for
        # each of the positions.
        for (index, id_label_pairs) in list(index_table.items()):
          conditions = ["is_{}_application(e)".format(label) for (_, label) in id_label_pairs]
          condition = " || ".join(conditions)
          projection_case = '''if ({0})
        {{
          return atermpp::down_cast<application>(e)[{1}];
        }}'''.format(condition, index)
          projection.append(projection_case)

      function = '''      ///\\brief Function for projecting out argument.
      ///        {0} from an application.
      /// \\param e A data expression.
      /// \\pre {0} is defined for e.
      /// \\return The argument of e that corresponds to {0}.
      inline
      const data_expression& {1}(const data_expression& e)
      {{
        assert({2});
        {3}
      }}'''.format(escape(p), p, " || ".join(assertions), "\n".join(projection))


      code.append(function)

    return "{0}\n\n".format("\n\n".join(code))

  def code(self, spec, is_constructor_declaration = False):
    assert(isinstance(spec, specification))
    # First we merge function declarations with the same function symbol/label,
    # hence we get function symbols with multiple sorts. Then for each of these
    # we generate the code
    class multi_function_declaration():
      def __init__(self, id, namespace, sort_expressions, l, internextern, definedby):
        assert(isinstance(id, identifier))
        assert(isinstance(sort_expressions, sort_expression_list))
        assert(isinstance(l, label))
        assert(isinstance(internextern, internal_external))
        assert(isinstance(definedby, defined_by))
        self.id = id
        self.namespace = namespace
        self.sort_expression_list = sort_expressions
        self.label = l
        self.internextern = internextern
        self.definedby = definedby

      def __eq__(self, other):
        return self.id == other.id and self.namespace == other.namespace and self.sort_expression_list == other.sort_expression_list and self.label == other.label

      def __str__(self):
        return "{0} : { {1} }".format(self.id, self.sort_expression_list)

      def function_name(self, fullname, name):
        CODE_TEMPLATE = Template('''
      /// \\brief Generate identifier ${namestring}.
      /// \\return Identifier ${namestring}.
      inline
      const core::identifier_string& ${namecode}_name()
      {
        static core::identifier_string ${namecode}_name = core::identifier_string("${name}");
        return ${namecode}_name;
      }
''')
        return CODE_TEMPLATE.substitute(namestring=escape(fullname), namecode=name, name=fullname)

      def function_constructor(self, fullname, name, sortparams, sort, number_of_data_parameters):
        CODE_TEMPLATE = Template('''
      /// \\brief Constructor for function symbol ${namestring}.
      ${sortparameterstring}
      /// \\return Function symbol ${functionname}.
      inline
      ${const}function_symbol${reference} ${functionname}(${sortparameters})
      {
        ${static}function_symbol ${functionname}(${functionname}_name(), ${sortname});
        return ${functionname};
      }
''')

        return self.function_name(fullname, name) + CODE_TEMPLATE.substitute(
          namestring=escape(fullname),
          sortparameterstring='\n      '.join(['/// \\param {0} A sort expression.'.format(escape(x.code())) for x in sortparams]),
          functionname=escape(name),
          sortparameters = ', '.join(['const sort_expression& {0}'.format(x.code()) for x in sortparams]),
          sortname=sort,
          const='const ' if sortparams == [] else '',
          reference='&' if sortparams == [] else '',
          static='static ' if sortparams == [] else '',
          parameters=", ".join(["a[%s]" %e for e in range(0,number_of_data_parameters)]))

      def polymorphic_function_constructor(self, fullname, name, sortparams):
        CODE_TEMPLATE = Template('''
      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol ${functionname}(${parameters})
      {
        ${targetsort}
        function_symbol ${functionname}(${functionname}_name(), ${sortname});
        return ${functionname};
      }
''')

        # There is polymorphism at play, hence we need to compute the
        # allowed sorts first.
        # That is, we have parameters for all domain elements
        # based on those, we compute the corresponding target sort
        domain_sort_ids = [sort_identifier(identifier("s%s" % i)) for i in range(len(self.sort_expression_list.elements[0].domain.elements))]
        target_sort_id = sort_identifier(identifier("target_sort"))
        new_sort = sort_arrow(domain(False, domain_sort_ids), target_sort_id)

        # If all codomains are equal we have a unique target sort, otherwise
        # this is detemined from the domain_sort_ids
        first_codomain = self.sort_expression_list.elements[0].codomain
        simple = all([x.codomain == first_codomain for x in self.sort_expression_list.elements])
        if simple:
          target_sort = 'sort_expression {0}({1});'.format(target_sort_id, first_codomain.code(spec))
        else:
          CASE_TEMPLATE = Template('''        ${elsestr}if (${condition})
        {
          target_sort = ${sort};
        }''')

          TARGET_SORT_TEMPLATE = Template('''sort_expression ${target_sort_id};
${cases}
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for ${functionname} with domain sorts \" + ${sortmsg} + \". \");
        }
''')

          cases = []
          for (i,sort) in enumerate(self.sort_expression_list.elements):
            cases.append(CASE_TEMPLATE.substitute(
              elsestr = '' if i == 0 else 'else ',
              condition = ' && '.join(['{0} == {1}'.format(domain_sort_ids[j_domsort[0]].code(spec), j_domsort[1].code(spec)) for j_domsort in enumerate(sort.domain.sorts())]),
              sort = sort.codomain.code(spec)
            ))

          target_sort = TARGET_SORT_TEMPLATE.substitute(
            target_sort_id = target_sort_id.code(spec),
            cases = '\n'.join(cases),
            functionname = self.label,
            sortmsg = " + \", \" + ".join(['pp({0})'.format(domain_sort_ids[j].code(spec)) for j in range(len(sort.domain.elements))])
            )

        parameters = []
        if domain_sort_ids != [] and sortparams != [] and simple:
          parameters += ['const sort_expression& ' for x in sortparams]
        else:
          parameters += ['const sort_expression& {0}'.format(x.code()) for x in sortparams]
        parameters += ['const sort_expression& {0}'.format(x.code()) for x in domain_sort_ids]

        return self.function_name(fullname, name) + CODE_TEMPLATE.substitute(
          namestring = escape(fullname),
          sortparameterstring = '\n      '.join(['/// \\param {0} A sort expression'.format(escape(x.code())) for x in sortparams + domain_sort_ids]),
          functionname = name,
          parameters = ', '.join(parameters),
          targetsort = target_sort,
          sortname = new_sort.code(spec)
          )

      def function_recogniser(self, fullname, name, sortparams):
        CODE_TEMPLATE = Template('''
      /// \\brief Recogniser for function ${namestring}.
      /// \\param e A data expression.
      /// \\return true iff e is the function symbol matching ${namestring}.
      inline
      bool is_${functionname}_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e)${getnamef} == ${functionname}${getname}();
        }
        return false;
      }
''')

        return CODE_TEMPLATE.substitute(
          namestring = escape(fullname),
          functionname = name,
          getnamef = '' if sortparams == [] else '.name()',
          getname = '' if sortparams == [] else '_name'
        )

      def polymorphic_function_recogniser(self, fullname, name, sortparams):
        CODE_TEMPLATE = Template('''
      /// \\brief Recogniser for function ${namestring}.
      /// \\param e A data expression.
      /// \\return true iff e is the function symbol matching ${namestring}.
      inline
      bool is_${functionname}_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == ${functionname}_name()${condition};
        }
        return false;
      }
''')

        domain_size = len(self.sort_expression_list.elements[0].domain.elements)
        condition = ' && atermpp::down_cast<function_sort>(f.sort()).domain().size() == {0}'.format(domain_size)

        if sortparams == '':
          cases = []
          for s in self.sort_expression_list.elements:
            d = s.domain
            assert(len(d.elements) == domain_size)
            cases.append('f == {0}({1})'.format(name, d.code(spec)))
          condition += ' && ({0})'.format(' || '.join(cases))

        return CODE_TEMPLATE.substitute(
          namestring = escape(fullname),
          functionname = name,
          condition = condition
        )

      def function_application(self, fullname, name, sort_params, data_params, polymorphic):
        CODE_TEMPLATE = Template('''
      /// \\brief Application of function symbol ${namestring}.
      ${sortparameterstring}
      ${dataparameterstring}
      /// \\return Application of ${namestring} to a number of arguments.
      inline
      application ${functionname}(${parameters})
      {
        return ${nsfunctionname}(${actsortparameters})(${actdataparameters});
      }

      /// \\brief Make an application of function symbol ${namestring}.
      /// \\param result The data expression where the ${namestring} expression is put.
      ${sortparameterstring}
      ${dataparameterstring}
      inline
      void make_${functionname}(data_expression& result, ${parameters})
      {
        make_application(result, ${nsfunctionname}(${actsortparameters}),${actdataparameters});
      }
''')

        formal_sort_params = ['const sort_expression& {0}'.format(fcode(x, spec)) for x in sort_params]
        formal_data_params = ['const data_expression& {0}'.format(fcode(x, spec)) for x in data_params]
        domain_params = ['{0}.sort()'.format(x) for x in data_params] if polymorphic else []
        return CODE_TEMPLATE.substitute(
          namestring = escape(fullname),
          sortparameterstring = '' if sort_params == '' else '\n      '.join(['/// \\param {0} A sort expression.'.format(fcode(x, spec)) for x in sort_params]),
          dataparameterstring = '' if data_params == '' else '\n      '.join(['/// \\param {0} A data expression.'.format(fcode(x, spec)) for x in data_params]),
          functionname = name,
          parameters = ', '.join(formal_sort_params + formal_data_params),
          nsfunctionname = add_namespace(name, self.namespace),
          actsortparameters = ', '.join([fcode(x, spec) for x in sort_params + domain_params]),
          actdataparameters = ', '.join([fcode(x, spec) for x in data_params])
        )

      def constant_implementation_code(self, fullname, name, sort_params, polymorphic):
        if not self.definedby.is_defined_by_code():
          return ""
        CODE_TEMPLATE = Template('''
      /// \\brief The data expression of an application of the constant symbol ${namestring}.
      /// \\details This function is to be implemented manually. \
      /// \\return The data expression corresponding to an application of ${namestring} to a number of arguments.
      inline
      void ${functionname}_manual_implementation(data_expression& result);

      /// \\brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void ${functionname}_application(data_expression& result, [[maybe_unused]] const data_expression& a)
      {
        assert(is_function_symbol(a));
        // assert(a==${functionname}());
        ${functionname}_manual_implementation(result${domain_parameters});
      }\n
''')

        formal_sort_params = ['const sort_expression& {0}'.format(fcode(x, spec)) for x in sort_params]
        # In the code  below data_params are used, but it is not declared nor defined... Long live Python.
        domain_params = ['{0}.sort()'.format(x) for x in data_params] if polymorphic else []
        return CODE_TEMPLATE.substitute(
          namestring = escape(fullname),
          sortparameterstring = '' if sort_params == '' else '\n      '.join(['/// \\param {0} A sort expression.'.format(fcode(x, spec)) for x in sort_params]),
          functionname = name,
          actsortparameters = ', '.join([fcode(x, spec) for x in sort_params + domain_params]),
          domain_parameters = ', '.join([fcode(x, spec) + ', ' for x in domain_params]),
          parameters = ', '.join(formal_sort_params),
        )

      def function_implementation_code(self, fullname, name, sort_params, data_params, polymorphic):
        if not self.definedby.is_defined_by_code():
          return ""
        CODE_TEMPLATE = Template('''
      /// \\brief The data expression of an application of the function symbol ${namestring}.
      /// \\details This function is to be implemented manually.
      ${dataparameterstring}
      /// \\return The data expression corresponding to an application of ${namestring} to a number of arguments.
      inline
      void ${functionname}_manual_implementation(data_expression& result, ${parameters});\n

      /// \\brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void ${functionname}_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==${functionname}());
        ${functionname}_manual_implementation(result, ${domain_parameters}${aaparameters});
      }\n
''')

        formal_sort_params = ['const sort_expression& {0}'.format(fcode(x, spec)) for x in sort_params]
        formal_data_params = ['const data_expression& {0}'.format(fcode(x, spec)) for x in data_params]
        domain_params = ['{0}.sort()'.format(x) for x in data_params] if polymorphic else []
        return CODE_TEMPLATE.substitute(
          namestring = escape(fullname),
          sortparameterstring = '' if sort_params == '' else '\n      '.join(['/// \\param {0} A sort expression.'.format(fcode(x, spec)) for x in sort_params]),
          dataparameterstring = '' if data_params == '' else '\n      '.join(['/// \\param {0} A data expression.'.format(fcode(x, spec)) for x in data_params]),
          functionname = name,
          parameters = ', '.join(formal_data_params),
          domain_parameters = ''.join([fcode(x, spec) + ', ' for x in domain_params]),
          actsortparameters = ', '.join([fcode(x, spec) for x in sort_params + domain_params]),
          aaparameters=", ".join(["a[%s]" %e for e in range(0,len(data_params))])
        )

      def function_application_recogniser(self, fullname, name):
        CODE_TEMPLATE = Template('''
      /// \\brief Recogniser for application of ${namestring}.
      /// \\param e A data expression.
      /// \\return true iff e is an application of function symbol ${functionname} to a
      ///     number of arguments.
      inline
      bool is_${functionname}_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_${functionname}_function_symbol(atermpp::down_cast<application>(e).head());
      }
''')

        return CODE_TEMPLATE.substitute(
          namestring = escape(fullname),
          functionname = name,
        )

      def function_application_code(self, sort, polymorphic = False):
        if not isinstance(sort, sort_arrow):
          return self.constant_implementation_code(self.id, self.label, self.sort_expression_list.sort_parameters(spec), polymorphic)

        return self.function_application(self.id, self.label, self.sort_expression_list.sort_parameters(spec), sort.domain.parameters(spec), polymorphic) + \
               self.function_application_recogniser(self.id, self.label) + \
               self.function_implementation_code(self.id, self.label, self.sort_expression_list.sort_parameters(spec), sort.domain.parameters(spec), polymorphic)


      def code(self, spec):
        assert(isinstance(spec, specification))

        if self.namespace != spec.namespace:
          return ''

        if len(self.sort_expression_list.elements) == 1:
          sort = self.sort_expression_list.elements[0] # as len is 1
          parameter_list_length = (len(sort.domain.parameters(spec)) if isinstance(sort,sort_arrow) else 0)
          return self.function_constructor(self.id, self.label, self.sort_expression_list.sort_parameters(spec), sort.code(spec), parameter_list_length) + \
                 self.function_recogniser(self.id, self.label, self.sort_expression_list.sort_parameters(spec)) + \
                 self.function_application_code(sort)

        else:
          assert len(self.sort_expression_list.elements) > 1

          return self.polymorphic_function_constructor(self.id, self.label, self.sort_expression_list.sort_parameters(spec)) + \
                 self.polymorphic_function_recogniser(self.id, self.label, self.sort_expression_list.formal_parameters_code(spec)) + \
                 self.function_application_code(self.sort_expression_list.elements[0], True)

    class multi_function_declaration_list():
      def __init__(self, elements):
        assert(all([isinstance(x, multi_function_declaration) for x in elements]))
        self.elements = elements

      def __eq__(self, other):
        return self.elements == other.elements

      def push_back(self, element):
        assert(isinstance(element, multi_function_declaration))
        self.elements += [element]

      def push_back_function_declaration(self, element):
        assert(isinstance(element, function_declaration))
        for e in self.elements:
          if (e.id == element.id) and (e.label == element.label):
            if all([str(x) != str(element.sort_expression) for x in e.sort_expression_list.elements]):
              e.sort_expression_list.push_back(element.sort_expression)
            return

        self.elements += [multi_function_declaration(element.id, element.namespace, sort_expression_list([element.sort_expression]), \
                                                     element.label, element.internextern, element.definedby)]

      def __str__(self):
        return "".join(["{0}\n".format(e) for e in self.elements])

      def code(self, spec):
        assert(isinstance(spec, specification))
        code = ""
        for e in self.elements:
          code += "%s" % (e.code(spec))
        return code

    # Add all the elements to a multi_function_declaration_list,
    # and let it do the code generation
    multi_function_declarations = multi_function_declaration_list([])
    for e in self.elements:
      if (e.namespace == self.namespace):
        multi_function_declarations.push_back_function_declaration(e)

    code = multi_function_declarations.code(spec)
    return code

class equation_declaration():
  def __init__(self, lhs, rhs, condition = None):
    assert(isinstance(lhs, data_expression))
    assert(isinstance(rhs, data_expression))
    self.namespace = ""
    self.original_namespace = ""
    if condition != None: assert(isinstance(condition, data_expression))
    self.lhs = lhs
    self.rhs = rhs
    self.condition = condition

  def __eq__(self, other):
    return self.condition == other.condition and self.lhs == other.lhs and self.rhs == other.rhs and self.condition == other.condition and self.namespace == other.namespace and self.original_namespace == other.original_namespace

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string

  def sort_parameters(self, spec, function_spec, variable_spec):
    condition_params = []
    if self.condition:
      condition_params = self.condition.sort_parameters(spec, function_spec, variable_spec)
    return merge([condition_params, self.lhs.sort_parameters(spec, function_spec, variable_spec), self.rhs.sort_parameters(spec, function_spec, variable_spec)])

  def has_lambda(self):
    if self.condition and self.condition.has_lambda():
      return True
    return self.lhs.has_lambda() or self.rhs.has_lambda()

  def has_forall(self):
    if self.condition and self.condition.has_forall():
      return True
    return self.lhs.has_forall() or self.rhs.has_forall()

  def has_exists(self):
    if self.condition and self.condition.has_exists():
      return True
    return self.lhs.has_exists() or self.rhs.has_exists()

  def __str__(self):
    if self.condition:
      return "{0} = {1}".format(self.lhs, self.rhs)
    else:
      return "{0} -> {1} = {2}".format(self.condition, self.lhs, self.rhs)

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    condition = None
    if self.condition != None:
      condition = self.condition.determinise_variable_or_function_symbol(function_spec, variable_spec)
    lhs = self.lhs.determinise_variable_or_function_symbol(function_spec, variable_spec)
    rhs = self.rhs.determinise_variable_or_function_symbol(function_spec, variable_spec)

    return equation_declaration(lhs, rhs, condition)

  def code(self, spec, function_spec, variable_spec):
    variables = []
    if self.condition:
      variables = self.condition.free_variables()
    variables = merge([variables, self.lhs.free_variables(), self.rhs.free_variables()])

    if len(variables) == 0:
      variables_string = "variable_list()"
    else:
      variables_string = "variable_list({" + "{0}".format(", ".join(sorted([v.code(spec, function_spec, variable_spec) for v in variables]))) + "})"

    if self.condition:
      return "result.push_back(data_equation({0}, {1}, {2}, {3}));".format(variables_string, self.condition.code(spec, function_spec, variable_spec), self.lhs.code(spec, function_spec, variable_spec), self.rhs.code(spec, function_spec, variable_spec))
    else:
      return "result.push_back(data_equation({0}, {1}, {2}));".format(variables_string, self.lhs.code(spec, function_spec, variable_spec), self.rhs.code(spec, function_spec, variable_spec))

class equation_declaration_list():
  def __init__(self, elements):
    assert(all([isinstance(x, equation_declaration) for x in elements]))
    self.elements = elements
    self.namespace = ""
    self.original_namespace = ""

  def __eq__(self, other):
    return self.elements == other.elements and self.namespace == other.namespace and self.original_namespace == other.original_namespace

  def push_back(self, element):
    assert(isinstance(element, equation_declaration))
    self.elements.append(element)

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    for e in self.elements:
      e.set_namespace(string)

  def sort_parameters(self, spec, function_spec, variable_spec):
    return merge([e.sort_parameters(spec, function_spec, variable_spec) for e in self.elements])

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    result = equation_declaration_list([])
    result.namespace = self.namespace
    result.original_namespace = self.original_namespace
    for e in self.elements:
      result.push_back(e.determinise_variable_or_function_symbol(function_spec, variable_spec))
    return result

  def has_lambda(self):
    return any([x.has_lambda() for x in self.elements])

  def has_forall(self):
    return any([x.has_forall() for x in self.elements])

  def has_exists(self):
    return any([x.has_exists() for x in self.elements])

  def __str__(self):
    return "".join(["{0}\n".format(e) for e in self.elements])

  def code(self, spec, function_spec, variable_spec):
    self = self.determinise_variable_or_function_symbol(function_spec, variable_spec)

    sort_parameters = self.sort_parameters(spec, function_spec, variable_spec)
    formal_parameters_code = [x.formal_parameter_code() for x in sort_parameters]

    code  = "      /// \\brief Give all system defined equations for %s\n" % (escape(spec.get_namespace()))
    for s in sort_parameters:
      code += "      /// \\param %s A sort expression\n" % (escape(s.code(spec)))
    code += "      /// \\return All system defined equations for sort %s\n" % (escape(spec.get_namespace()))
    code += "      inline\n"
    code += "      data_equation_vector %s_generate_equations_code(%s)\n" % (spec.get_namespace(), ", ".join(formal_parameters_code))
    code += "      {\n"
    code += "%s\n" % (variable_spec.code(spec))
    code += "        data_equation_vector result;\n"
    code += "%s" % (spec.sort_specification.structured_sort_equation_code())
    for e in self.elements:
      code += "        %s\n" % (e.code(spec, function_spec, variable_spec))
    code += "        return result;\n"
    code += "      }\n"

    return code


class data_expression():
  def __init__(self):
    pass

  def __str__(self):
    raise Exception("String operation for data_expression is not implemented")

  def has_lambda(self):
    return False

  def has_forall(self):
    return False

  def has_exists(self):
    return False

class data_application(data_expression):
  def __init__(self, head, arguments):
    assert(isinstance(head, data_expression))
    assert(isinstance(arguments, data_expression_list))
    self.head = head
    self.arguments = arguments

  def __eq__(self, other):
    if other==None:
       return False
    return self.head == other.head and self.arguments == other.arguments

  def sort_parameters(self, spec, function_spec, variable_spec):
    return merge([self.head.sort_parameters(spec, function_spec, variable_spec), self.arguments.sort_parameters(spec, function_spec, variable_spec)])

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    return data_application(self.head.determinise_variable_or_function_symbol(function_spec, variable_spec), self.arguments.determinise_variable_or_function_symbol(function_spec, variable_spec))

  def free_variables(self):
    return merge([self.head.free_variables(), self.arguments.free_variables()])

  def has_lambda(self):
    return self.head.has_lambda() or self.arguments.has_lambda()

  def has_forall(self):
    return self.head.has_forall() or self.arguments.has_forall()

  def has_exists(self):
    return self.head.has_exists() or self.arguments.has_exists()

  def __str__(self):
    return "{0}({1})".format(self.head, self.arguments)

  def code(self, spec, function_spec, variable_spec, indentlog = 0):
    LOG.debug(indent(indentlog) + "Generating code for {0}".format(self))
    if isinstance(self.head, function_symbol):
      # FIXME: Extremely ugly workaround for set complement
      if self.head.id == "!" and len(self.arguments.elements) == 1 and isinstance(self.arguments.elements[0], data_application) and isinstance(self.arguments.elements[0].head, data_variable):
        head_code = "sort_bool_::not_"
      # FIXME: Extremely ugly workaround for bag join
      elif self.head.id == "+" and len(self.arguments.elements) == 2 and isinstance(self.arguments.elements[0], data_application) and isinstance(self.arguments.elements[0].head, data_variable):
        head_code = "sort_nat::plus"
      else:
        head_code = self.head.code(spec, function_spec, variable_spec, indentlog + 2, len(self.arguments.elements))
        # make generated code a bit more readable
        head_code = head_code[:head_code.find("(")]
    else:
      head_code = self.head.code(spec, function_spec, variable_spec, indentlog + 2)

    sort_parameters_code = []

    # FIXME: workaround for complement and join
    if isinstance(self.head, function_symbol) and (head_code != "sort_bool_::not_" and head_code != "sort_nat::plus"):
      sort_parameters = self.head.sort_parameters(spec, function_spec, variable_spec)
      sort_parameters_code += [x.code(spec) for x in sort_parameters]

    sort_parameters_code += [self.arguments.code(spec, function_spec, variable_spec, indentlog + 2)]

    result = "%s(%s)" % (head_code, ", ".join(sort_parameters_code))
    LOG.debug(indent(indentlog) + "Yields {0}".format(result))
    return result

class data_variable_or_function_symbol(data_expression):
  def __init__(self, identifier):
    self.id = identifier

  def __eq__(self, other):
    return self.id == other.id

  def __str__(self):
    return str(self.id)

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    result = function_spec.find_function(self, -1)
    if result:
      return function_symbol(self.id)
    else:
      return data_variable(self.id)

class data_variable(data_expression):
  def __init__(self, id):
    assert(isinstance(id, identifier))
    self.id = id

  def __eq__(self, other):
    if hasattr(other, "id"):
      return self.id == other.id
    else:
      return NotImplemented

  def __lt__(self, other):
    if hasattr(other, "id"):
      return self.id < other.id
    else:
      return NotImplemented

  def sort_parameters(self, spec, function_spec, variable_spec):
    variable = variable_spec.find_variable(self)
    if variable == None:
      return []

    return variable.sort_expression.sort_parameters(spec)

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    return self

  def free_variables(self):
    return [self]

  def __str__(self):
    return str(self.id)

  def code(self, spec, function_spec, variable_spec, indentlog = 0):
    LOG.debug(indent(indentlog) + "Generating code for {0}".format(self))
    result = "v{0}".format(self.id)
    LOG.debug(indent(indentlog) + "Yields {0}".format(result))
    return result

class function_symbol(data_expression):
  def __init__(self, id):
    assert(isinstance(id, identifier))
    self.id = id

  def __eq__(self, other):
    return self.id == other.id

  def sort_parameters(self, spec, function_spec, variable_spec):
    return function_spec.find_function(self, -1).sort_expression.sort_parameters(spec)

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    return self

  def free_variables(self):
    result = set()
    return result

  def __str__(self):
    return str(self.id)

  def code(self, spec, function_spec, variable_spec, indentlog=0, argumentcount = -1):
    LOG.debug(indent(indentlog) + "Generating code for {0} with argumentcount {1} and indent {2}".format(self, argumentcount, indentlog))
    f = function_spec.find_function(self, argumentcount)
    sort_parameters_code = [s.code(spec) for s in f.sort_parameters(spec)]
    result = "%s(%s)" % (add_namespace(f.label, f.namespace, function_spec.namespace), ", ".join(sort_parameters_code))
    LOG.debug(indent(indentlog) + "Yields {0}".format(result))
    return result

class lambda_abstraction(data_expression):
  def __init__(self, var_declaration, expression):
    assert(isinstance(var_declaration, variable_declaration))
    assert(isinstance(expression, data_expression))
    self.variable_declaration = var_declaration
    self.expression = expression

  def __eq__(self, other):
    return self.variable_declaration == other.variable_declaration and self.expression == other.expression

  def sort_parameters(self, spec, function_spec, variable_spec):
    return merge(self.variable_declaration.sort_parameters(spec), self.expression.sort_parameters(spec, function_spec, variable_spec))

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    variable_spec.declarations.push_back(self.variable_declaration)
    return lambda_abstraction(self.variable_declaration, self.expression.determinise_variable_or_function_symbol(function_spec, variable_spec))

  def free_variables(self):
    return difference(self.expression.free_variables(), [data_variable(self.variable_declaration.id)])

  def has_lambda(self):
    return True

  def has_forall(self):
    return self.expression.has_forall()

  def has_exists(self):
    return self.expression.has_exists()

  def __str__(self):
    return "lambda({0}, {1})".format(self.variable_declaration, self.expression)

  def code(self, spec, function_spec, variable_spec, indentlog = 0):
    LOG.debug(indent(indentlog) + "Generating code for {0}".format(self))
    result = "lambda(variable_list({%s}), %s)" % (data_variable(self.variable_declaration.id).code(spec, function_spec, variable_spec, indentlog+2), self.expression.code(spec, function_spec, variable_spec, indentlog+2))
    LOG.debug(indent(indentlog) + "Yields {0}".format(result))
    return result

class forall(data_expression):
  def __init__(self, var_declaration, expression):
    assert(isinstance(var_declaration, variable_declaration))
    assert(isinstance(expression, data_expression))
    self.variable_declaration = var_declaration
    self.expression = expression

  def __eq__(self, other):
    return self.variable_declaration == other.variable_declaration and self.expression == other.expression

  def sort_parameters(self, spec, function_spec, variable_spec):
    return merge([self.variable_declaration.sort_parameters(spec), self.expression.sort_parameters(spec, function_spec, variable_spec)])

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    variable_spec.declarations.push_back(self.variable_declaration)
    return forall(self.variable_declaration, self.expression.determinise_variable_or_function_symbol(function_spec, variable_spec))

  def free_variables(self):
    return difference(self.expression.free_variables(), [data_variable(self.variable_declaration.id)])

  def has_lambda(self):
    return self.expression.has_lambda()

  def has_forall(self):
    return True

  def has_exists(self):
    return self.expression.has_exists()

  def __str__(self):
    return "forall({0}, {1})".format(self.variable_declaration, self.expression)

  def code(self, spec, function_spec, variable_spec, indentlog = 0):
    LOG.debug(indent(indentlog) + "Generating code for {0}".format(self))
    result = "forall(variable_list({%s}), %s)" % (data_variable(self.variable_declaration.id).code(spec, function_spec, variable_spec, indentlog+2), self.expression.code(spec, function_spec, variable_spec, indentlog+2))
    LOG.debug(indent(indentlog) + "Yields {0}".format(result))
    return result

class exists(data_expression):
  def __init__(self, var_declaration, expression):
    assert(isinstance(var_declaration, variable_declaration))
    assert(isinstance(expression, data_expression))
    self.variable_declaration = var_declaration
    self.expression = expression

  def __eq__(self, other):
    return self.variable_declaration == other.variable_declaration and self.expression == other.expression

  def sort_parameters(self, spec, function_spec, variable_spec):
    return merge(self.variable_declaration.sort_parameters(spec), self.expression.sort_parameters(spec, function_spec, variable_spec))

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    variable_spec.declarations.push_back(self.variable_declaration)
    return exists(self.variable_declaration, self.expression.determinise_variable_or_function_symbol(function_spec, variable_spec))

  def free_variables(self):
    return difference(self.expression.free_variables(), [data_variable(self.variable_declaration.id)])

  def has_lambda(self):
    return self.expression.has_lambda()

  def has_forall(self):
    return self.expression.has_forall()

  def has_exists(self):
    return True

  def __str__(self):
    return "exists({0}, {1})".format(self.variable_declaration, self.expression)

  def code(self, spec, function_spec, variable_spec, indentlog=0):
    LOG.debug(indent(indentlog) + "Generating code for {0}".format(self))
    result = "exists(variable_list({%s}), %s)" % (data_variable(self.variable_declaration.id).code(spec, function_spec, variable_spec, indentlog+2), self.expression.code(spec, function_spec, variable_spec, indentlog+2))
    LOG.debug(indent(indentlog) + "Yields {0}".format(result))
    return result

class data_expression_list():
  def __init__(self, elements):
    assert(all([isinstance(x, data_expression) for x in elements]))
    self.elements = elements

  def __eq__(self, other):
    return self.elements == other.elements

  def push_back(self, element):
    assert(isinstance(element, data_expression))
    self.elements += [element]

  def determinise_variable_or_function_symbol(self, function_spec, variable_spec):
    result = data_expression_list([])
    for e in self.elements:
      result.push_back(e.determinise_variable_or_function_symbol(function_spec, variable_spec))
    return result

  def sort_parameters(self, spec, function_spec, variable_spec):
    return merge([e.sort_parameters(spec, function_spec, variable_spec) for e in self.elements])

  def free_variables(self):
    return merge([e.free_variables() for e in self.elements])

  def has_lambda(self):
    return any([x.has_lambda() for x in self.elements])

  def has_forall(self):
    return any([x.has_forall() for x in self.elements])

  def has_exists(self):
    return any([x.has_exists() for x in self.elements])

  def __str__(self):
    return ", ".join([str(e) for e in self.elements])

  def code(self, spec, function_spec, variable_spec, indentlog):
    return ", ".join([e.code(spec, function_spec, variable_spec, indentlog+2) for e in self.elements])

class sort_expression_list():
  def __init__(self, elements):
    assert(all([isinstance(x, sort_expression) for x in elements]))
    self.elements = elements

  def __eq__(self, other):
    return self.elements == other.elements

  def push_back(self, element):
    assert(isinstance(element, sort_expression))
    self.elements += [element]

  def sort_parameters(self, spec):
    return merge([e.sort_parameters(spec) for e in self.elements])

  def formal_parameters_code(self, spec):
    return ", ".join([p.formal_parameter_code() for p in self.sort_parameters(spec)])

  def actual_parameters_code(self, spec):
    return ", ".join([p.actual_parameter_code() for p in self.sort_parameters(spec)])

  def __str__(self):
    return ", ".join(self.elements)

  def code(self):
    return "variable_list({" + "{0}".format(", ".join([e.code() for e in self.elements()])) + "})"

class sort_expression():
  def __init__(self):
    pass

  def sort_parameters(self, spec):
    return []

class sort_identifier(sort_expression):
  def __init__(self, id):
    assert(isinstance(id, identifier))
    self.name = id

  def __eq__(self, other):
    if hasattr(other, "name"):
      return self.name == other.name
    else:
      return NotImplemented

  def __lt__(self, other):
    if hasattr(other, "name"):
      return self.name < other.name
    else:
      return NotImplemented

  def sort_parameters(self, spec):
    if self in spec.sort_parameters():
      return [self]
    else:
      return []

  def __str__(self):
    return str(self.name)

  def formal_parameter_code(self):
    return "const sort_expression& {0}".format(self.name).lower()

  def actual_parameter_code(self):
    return str(self.name).lower()

  def code(self, spec = None):
    #assert(isinstance(spec, specification))
    if (spec != None and spec.sort_specification.has_sort(self)):
      return spec.sort_specification.get_sort(self).inline_code(spec)
    else:
      return str(self.name).lower()

class domain(sort_expression):
  def __init__(self, labelled, elements):
    assert(isinstance(labelled, bool))
    self.labelled = labelled
    self.elements = elements

  def __eq__(self, other):
    if hasattr(other, "labelled") and hasattr(other, "elements"):
      return self.labelled == other.labelled and self.elements == other.elements
    else:
      return NotImplemented

  def sorts(self):
    if self.labelled:
      return [e[0] for e in self.elements]
    else:
      return self.elements

  def __labels(self):
    return [e[1] for e in self.elements]

  def push_back(self, element):
    self.elements.append(element)

  def projection_arguments(self):
    if not self.labelled:
      return []
    else:
      return [(e, i) for (i, e) in enumerate(self.__labels())]

  def sort_parameters(self, spec):
    return merge([p.sort_parameters(spec) for p in self.sorts()])

  def parameters(self, spec):
    return ['arg{0}'.format(i) for i in range(0,len(self.elements))]

  def formal_parameters_code(self, spec):
    return ', '.join(['const data_expression& {0}'.format(str(x)) for x in self.parameters(spec)])

  def actual_parameters_code(self, spec):
    return ', '.join(self.parameters(spec))

  def assertions_code(self):
    return ""

  def __str__(self):
    if self.labelled:
      return " # ".join(["{0} <\"{1}\">".format(e[0], e[1]) for e in self.elements])
    else:
      return " # ".join(map(str, self.elements))

  def code(self, spec):
    return ", ".join([e.code(spec) for e in self.sorts()])

class sort_arrow(sort_expression):
  def __init__(self, dom, codom):
    assert(isinstance(dom, domain))
    assert(isinstance(codom, sort_expression))
    self.domain = dom
    self.codomain = codom

  def __eq__(self, other):
    return self.domain == other.domain and self.codomain == other.codomain

  def sort_parameters(self, spec):
    return merge([self.codomain.sort_parameters(spec), self.domain.sort_parameters(spec)])

  def __str__(self):
    return "{0} -> {1}".format(self.domain, self.codomain)

  def code(self, spec):
    assert(isinstance(spec, specification))
    return "make_function_sort_(%s, %s)" % (self.domain.code(spec), self.codomain.code(spec))

class sort_container(sort_expression):
  def __init__(self, container, element_sort):
    assert(isinstance(container, identifier))
    assert(isinstance(element_sort, sort_expression))
    self.container = container
    self.element_sort = element_sort

  def __eq__(self, other):
    if hasattr(other, "container") and hasattr(other, "element_sort"):
      return self.container == other.container and self.element_sort == other.element_sort
    else:
      return NotImplemented

  def sort_parameters(self, spec):
    return self.element_sort.sort_parameters(spec)

  def __str__(self):
    return "{0}({1})".format(self.container, self.element_sort)

  def code(self, spec):
    assert(isinstance(spec, specification))
    if (spec.sort_specification.has_sort(self)):
      sort = spec.sort_specification.get_sort(self)
      result = "%s(%s)" % (sort.label, self.element_sort.code(spec))
      if sort.original_namespace != spec.namespace and sort.original_namespace != "undefined":
        result = "sort_%s::%s" % (sort.original_namespace, result)
      return result
    else:
      return "%s(%s)" % (str(self.container).lower(), self.element_sort.code(spec))

class structured_sort_declaration():
  def __init__(self, id, l, arguments = None):
    assert(isinstance(id, identifier))
    assert(isinstance(l, label))
    assert(arguments == None or isinstance(arguments, domain))
    self.id = id
    self.label = l
    self.arguments = arguments
    self.namespace = ""
    self.original_namespace = ""

  def __eq__(self, other):
    return self.id == other.id and self.label == other.label and self.arguments == other.arguments and self.namespace == other.namespace and self.original_namespace == other.original_namespace

  def merge_structured_sort(self, function_spec, sort_expr):
    if self.arguments == None:
      s = sort_expr
    else:
      s = sort_arrow(self.arguments, sort_expr)
    f = function_declaration(self.id, s, self.label,internal_external("external"), defined_by("defined_by_rewrite_rules"))
    # f.set_namespace should not be used here, as this function may be
    # introduced in an included specification, and hence we need to preserve
    # original_namespace
    f.original_namespace = self.original_namespace
    f.namespace = self.namespace
    constructor_spec = function_spec.constructor_specification
    constructor_spec.declarations.push_back(f)
    function_spec.constructor_specification = constructor_spec
    return function_spec

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string

  def __str__(self):
    if self.arguments:
      return "{0}({1})".format(self.id, self.arguments)
    else:
      return str(self.id)

  def struct_constructor_arguments(self, spec):
    return "variable_list({" + "{0}".format(", ".join(["structured_sort_constructor_argument(\"%s\", %s)" % (a[1], a[0].code(spec)) for a in self.arguments.elements])) + "})"

  def code(self, spec):
    if self.arguments == None:
      return "structured_sort_constructor(\"%s\", \"%s\")" % (self.id, self.label)
    else:
      return "structured_sort_constructor(\"%s\", %s, \"%s\")" % (self.id, self.struct_constructor_arguments(spec), self.label)

class structured_sort_declaration_list():
  def __init__(self, elements):
    assert(all([isinstance(x, structured_sort_declaration) for x in elements]))
    self.elements = elements
    self.namespace = ""
    self.original_namespace = ""

  def __eq__(self, other):
    return self.elements == other.elements and self.namespace == other.namespace and self.original_namespace == other.original_namespace

  def merge_structured_sort(self, function_spec, sort_expr):
    for e in self.elements:
      function_spec = e.merge_structured_sort(function_spec, sort_expr)
    return function_spec

  def push_back(self, element):
    assert(isinstance(element, structured_sort_declaration))
    self.elements += [element]

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    for e in self.elements:
      e.set_namespace(string)

  def __str__(self):
    return " | ".join(map(str, self.elements))

  def code(self, spec):
    code = ""
    for e in self.elements:
      code += "          constructors.push_back(%s);\n" % (e.code(spec))
    return code

class structured_sort_specification():
  def __init__(self, elements):
    assert(isinstance(elements, structured_sort_declaration_list))
    self.elements = elements
    self.namespace = ""
    self.original_namespace = ""

  def __eq__(self, other):
    if hasattr(other, "elements") and hasattr(other, "namespace") and hasattr(other, "original_namespace"):
      return self.elements == other.elements and self.namespace == other.namespace and self.original_namespace == other.original_namespace
    else:
      return NotImplemented

  def merge_structured_sort(self, function_spec, sort_expr):
    return self.elements.merge_structured_sort(function_spec, sort_expr)

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    self.elements.set_namespace(string)

  def __str__(self):
    return "struct {0}".format(self.elements)

  def code(self, spec):
    return self.elements.code(spec)

class sort_declaration():
  def __init__(self, sort_expr, l, alias = None):
    assert(isinstance(sort_expr, sort_expression))
    assert(alias == None or isinstance(alias, structured_sort_specification))
    assert(isinstance(l, label))
    self.sort_expression = sort_expr
    self.label = l
    self.alias = alias
    self.namespace = ""
    self.original_namespace = ""

  def merge_structured_sort(self, function_spec):
    if self.alias != None:
      function_spec = self.alias.merge_structured_sort(function_spec, self.sort_expression)
    return function_spec

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    if self.alias != None:
      self.alias.set_namespace(string)

  def name(self):
    return self.sort_expression

  def defines_container(self):
    return isinstance(self.sort_expression, sort_container)

  def defines_struct(self):
    return self.alias != None

  def __str__(self):
    if self.alias:
      return "{0} = {1}".format(self.sort_expression, self.alias)
    else:
      return str(self.sort_expression)

  def inline_code(self, spec):
    if (self.namespace == spec.namespace) or (self.namespace == "") or self.namespace == "undefined":
      return "%s()" % (self.label)
    else:
      return "sort_%s::%s()" % (self.namespace, self.label)

  def structured_sort_constructor_code(self):
    if self.alias == None:
      return ""
    else:
      param = ""
      if isinstance(self.sort_expression, sort_container):
        param = str(self.sort_expression.element_sort).lower()
      sort = "%s(%s)" % (self.label, param)
      result = "        function_symbol_vector %s_constructors = detail::%s_struct(%s).constructor_functions(%s);\n" % (self.label, self.label, param, sort)
      return result + "        result.insert(result.end(), %s_constructors.begin(), %s_constructors.end());\n" % (self.label, self.label)

  def structured_sort_mapping_code(self):
    if self.alias == None:
      return ""
    else:
      param = ""
      if isinstance(self.sort_expression, sort_container):
        param = str(self.sort_expression.element_sort).lower()
      sort = "%s(%s)" % (self.label, param)
      result = "        function_symbol_vector %s_mappings = detail::%s_struct(%s).comparison_functions(%s);\n" % (self.label, self.label, param, sort)
      result += "        result.insert(result.end(), %s_mappings.begin(), %s_mappings.end());\n" % (self.label, self.label)
#      result += "        %s_mappings = detail::%s_struct(%s).projection_functions(%s);\n" % (self.label, self.label, param, sort)
#      result += "        result.insert(result.end(), %s_mappings.begin(), %s_mappings.end());\n" % (self.label, self.label)
#      result += "        %s_mappings = detail::%s_struct(%s).recogniser_functions(%s);\n" % (self.label, self.label, param, sort)
#      result += "        result.insert(result.end(), %s_mappings.begin(), %s_mappings.end());\n" % (self.label, self.label)
      return result

  def structured_sort_equation_code(self):
    if self.alias == None:
      return ""
    else:
      param = ""
      if isinstance(self.sort_expression, sort_container):
        param = str(self.sort_expression.element_sort).lower()
      sort = "%s(%s)" % (self.label, param)
      result = "        data_equation_vector %s_equations = detail::%s_struct(%s).constructor_equations(%s);\n" % (self.label, self.label, param, sort)
      result += "        result.insert(result.end(), %s_equations.begin(), %s_equations.end());\n" % (self.label, self.label)
      result += "        %s_equations = detail::%s_struct(%s).comparison_equations(%s);\n" % (self.label, self.label, param, sort)
      result += "        result.insert(result.end(), %s_equations.begin(), %s_equations.end());\n" % (self.label, self.label)
#      result += "        %s_equations = detail::%s_struct(%s).projection_equations(%s);\n" % (self.label, self.label, param, sort)
#      result += "        result.insert(result.end(), %s_equations.begin(), %s_equations.end());\n" % (self.label, self.label)
#      result += "        %s_equations = detail::%s_struct(%s).recogniser_equations(%s);\n" % (self.label, self.label, param, sort)
#      result += "        result.insert(result.end(), %s_equations.begin(), %s_equations.end());\n" % (self.label, self.label)
      return result

  def sort_name(self, id, label):
    code = ""
    code += "      inline\n"
    code += "      const core::identifier_string& %s_name()\n" % (label)
    code += "      {\n"
    code += "        static core::identifier_string %s_name = core::identifier_string(\"%s\");\n" % (label, id)
    code += "        return %s_name;\n" % (label)
    code += "      }\n\n"
    return code

  def container_name(self, id, label):
    code = ""
    code += "      inline\n"
    code += "      const core::identifier_string& %s_name()\n" % (label)
    code += "      {\n"
    code += "        static core::identifier_string %s_name = core::identifier_string(\"%s\");\n" % (label, label)
    code += "        return %s_name;\n" % (label)
    code += "      }\n\n"
    return code

  def sort_expression_constructors(self, id, label):
    code = ""
    code += self.sort_name(id, label)
    code += "      /// \\brief Constructor for sort expression %s.\n" % (escape(id))
    code += "      /// \\return Sort expression %s.\n" % (escape(id))
    code += "      inline\n"
    code += "      const basic_sort& %s()\n" % (label)
    code += "      {\n"
    code += "        static basic_sort %s = basic_sort(%s_name());\n" % (label, label)
    code += "        return %s;\n" % (label)
    code += "      }\n\n"

    code += "      /// \\brief Recogniser for sort expression %s\n" % (escape(str(id)))
    code += "      /// \\param e A sort expression\n"
    code += "      /// \\return true iff e == %s()\n" % (escape(label))
    code += "      inline\n"
    code += "      bool is_%s(const sort_expression& e)\n" % (label)
    code += "      {\n"
    code += "        if (is_basic_sort(e))\n"
    code += "        {\n"
    code += "          return basic_sort(e) == %s();\n" % (label)
    code += "        }\n"
    code += "        return false;\n"
    code += "      }\n"
    return code

  def sort_expression_constructors_container_sort(self, id, label, parameter):
    code = ""
#    code += self.container_name(id, label)
    code += "      /// \\brief Constructor for sort expression %s(%s)\n" % (escape(id), escape(parameter))
    code += "      /// \\param %s A sort expression\n" % (escape(str(parameter).lower()))
    code += "      /// \\return Sort expression %s(%s)\n" % (escape(label), escape(str(parameter).lower()))
    code += "      inline\n"
    code += "      container_sort %s(const sort_expression& %s)\n" % (label, str(parameter).lower())
    code += "      {\n"
    code += "        container_sort %s(%s_container(), %s);\n" % (label, label, str(parameter).lower())
    code += "        return %s;\n" % (label)
    code += "      }\n\n"

    code += "      /// \\brief Recogniser for sort expression %s(%s)\n" % (escape(id), escape(str(parameter).lower()))
    code += "      /// \\param e A sort expression\n"
    code += "      /// \\return true iff e is a container sort of which the name matches\n"
    code += "      ///      %s\n" % (escape(label))
    code += "      inline\n"
    code += "      bool is_%s(const sort_expression& e)\n" % (label)
    code += "      {\n"
    code += "        if (is_container_sort(e))\n"
    code += "        {\n"
    code += "          return container_sort(e).container_name() == %s_container();\n" % (label)
    code += "        }\n"
    code += "        return false;\n"
    code += "      }\n"
    return code

  def code(self, spec):
    if isinstance(self.sort_expression, sort_identifier):
      code = self.sort_expression_constructors(self.sort_expression, self.label)
    else:
      assert(isinstance(self.sort_expression, sort_container))
      code = self.sort_expression_constructors_container_sort(self.sort_expression.container, self.label, self.sort_expression.element_sort)

    if self.alias != None:
      code += "\n"
      code += "      namespace detail {\n\n"
      code += "        /// \\brief Declaration for sort %s as structured sort\n" % (escape(self.label))
      if isinstance(self.sort_expression, sort_container):
        param = str(self.sort_expression.element_sort).lower()
        code += "        /// \\param %s A sort expression\n" % (escape(param))
      else:
        param = ""
      code += "        /// \\return The structured sort representing %s\n" % (escape(self.label))
      code += "        inline\n"
      code += "        structured_sort %s_struct(const sort_expression& %s)\n" % (self.label, param)
      code += "        {\n"
      code += "          structured_sort_constructor_vector constructors;\n"
      code += "%s" % (self.alias.code(spec))
      code += "          return structured_sort(constructors);\n"
      code += "        }\n\n"
      code += "      } // namespace detail\n"
    return code

class sort_declaration_list():
  def __init__(self, elements):
    assert(all([isinstance(x, sort_declaration) for x in elements]))
    self.elements = elements
    self.namespace = ""
    self.original_namespace = ""

  def merge_structured_sorts(self, function_spec):
    for e in self.elements:
      function_spec = e.merge_structured_sort(function_spec)
    return function_spec

  def push_back(self, element):
    assert(isinstance(element, sort_declaration))
    self.elements += [element]

  def is_alias(self, sort_expr):
    assert(isinstance(sort_expr, sort_expression))
    return any([x.sort_expression == sort_expr and x.alias != None for x in self.elements])

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    for e in self.elements:
      e.set_namespace(string)

  def merge_declarations(self, declarations, is_supertype):
    self.elements += declarations.elements

  def defines_container(self):
    return any([x.defines_container() for x in self.elements])

  def defines_struct(self):
    return any([x.defines_struct() for x in self.elements])

  def __str__(self):
    return "".join(["{0}\n".format(e) for e in self.elements])

  def structured_sort_constructor_code(self):
    code = ""
    for e in self.elements:
      if e.namespace == self.namespace:
        code += "%s" % (e.structured_sort_constructor_code())
    return code

  def structured_sort_mapping_code(self):
    code = ""
    for e in self.elements:
      if e.namespace == self.namespace:
        code += "%s" % (e.structured_sort_mapping_code())
    return code

  def structured_sort_equation_code(self):
    code = ""
    for e in self.elements:
      if e.namespace == self.namespace:
        code += "%s" % (e.structured_sort_equation_code())
    return code

  def code(self, spec):
    code = ""
    for e in self.elements:
      if e.namespace == self.namespace:
        code += "%s\n" % (e.code(spec))
    return code

class equation_specification():
  def __init__(self, declarations):
    assert(isinstance(declarations, equation_declaration_list))
    self.declarations = declarations
    self.namespace = ""
    self.original_namespace = ""

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    self.declarations.set_namespace(string)

  def has_lambda(self):
    return self.declarations.has_lambda()

  def has_forall(self):
    return self.declarations.has_forall()

  def has_exists(self):
    return self.declarations.has_exists()

  def __str__(self):
    return "eqn {0}".format(self.declarations)

  def code(self, spec, function_spec, variable_spec):
    return self.declarations.code(spec, function_spec, variable_spec)

class variable_specification():
  def __init__(self, declarations):
    assert(isinstance(declarations, variable_declaration_list))
    self.declarations = declarations

  def find_variable(self, variable):
    return self.declarations.find_variable(variable)

  def __str__(self):
    return "var {0}".format(self.declarations)

  def code(self, spec):
    return self.declarations.code(spec)

class mapping_specification():
  def __init__(self, declarations):
    assert(isinstance(declarations, function_declaration_list))
    self.declarations = declarations
    self.namespace = ""
    self.original_namespace = ""

  def has_cplusplus_implementable_code(self):
    return self.declarations.has_cplusplus_implementable_code()

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    self.declarations.set_namespace(string)

  def find_function(self, function, argumentcount):
    return self.declarations.find_function(function, argumentcount)

  def __str__(self):
    return "map {0}".format(self.declarations)

  def merge_specification(self, spec, is_supertype):
    self.declarations.merge_declarations(spec.declarations, is_supertype)

  def code(self, spec):
    if not self.declarations.empty():
      assert(isinstance(spec, specification))
      sort_parameters = ", ".join(["const sort_expression& %s" % (str(x).lower()) for x in self.declarations.sort_parameters(spec)])
      data_parameters = ", ".join(["%s" % (str(x).lower()) for x in self.declarations.sort_parameters(spec)])
      constructor_data_parameters = ", ".join(["%s" % (str(x).lower()) for x in spec.function_specification.constructor_specification.declarations.sort_parameters(spec)])
      assert(self.namespace == spec.namespace)
      namespace_string = self.namespace
      if namespace_string == "undefined":
        namespace_string = spec.get_namespace()
      code = ""
      code += self.declarations.code(spec)
      code += "      /// \\brief Give all system defined mappings for %s\n" % (escape(namespace_string))
      for s in self.declarations.sort_parameters(spec):
        code += "      /// \\param %s A sort expression\n" % (escape(str(s).lower()))
      code += "      /// \\return All system defined mappings for %s\n" % (escape(namespace_string))
      code += "      inline\n"
      code += "      function_symbol_vector %s_generate_functions_code(%s)\n" % (namespace_string, sort_parameters)
      code += "      {\n"
      code += "        function_symbol_vector result;\n"
      code += self.declarations.generator_code(spec) + (spec.sort_specification.structured_sort_mapping_code())
      code += "        return result;\n"
      code += "      }\n"
      code += "      \n"

      code += "      /// \\brief Give all system defined mappings and constructors for %s\n" % (escape(namespace_string))
      for s in self.declarations.sort_parameters(spec):
        code += "      /// \\param %s A sort expression\n" % (escape(str(s).lower()))
      code += "      /// \\return All system defined mappings for %s\n" % (escape(namespace_string))
      code += "      inline\n"
      code += "      function_symbol_vector %s_generate_constructors_and_functions_code(%s)\n" % (namespace_string, sort_parameters)
      code += "      {\n"
      code += "        function_symbol_vector result=%s_generate_functions_code(%s);\n" % (namespace_string, data_parameters)
      if not spec.function_specification.constructor_specification.declarations.empty():
        code += "        for(const function_symbol& f: %s_generate_constructors_code(%s))\n" % (namespace_string, constructor_data_parameters)
        code += "        {\n"
        code += "          result.push_back(f);\n"
        code += "        }\n"
      code += "        return result;\n"
      code += "      }\n"
      code += "      \n"

      code += "      /// \\brief Give all system defined mappings that can be used in mCRL2 specs for %s\n" % (escape(namespace_string))
      for s in self.declarations.sort_parameters(spec):
        code += "      /// \\param %s A sort expression\n" % (escape(str(s).lower()))
      code += "      /// \\return All system defined mappings for that can be used in mCRL2 specificationis %s\n" % (escape(namespace_string))
      code += "      inline\n"
      code += "      function_symbol_vector %s_mCRL2_usable_mappings(%s)\n" % (namespace_string, sort_parameters)
      code += "      {\n"
      code += "        function_symbol_vector result;\n"
      add_mappings_code = self.declarations.mCRL2_usable_functions(spec) + (spec.sort_specification.structured_sort_mapping_code())
      if add_mappings_code != "":
        code += add_mappings_code
      code += "        return result;\n"
      code += "      }\n"
      code += "\n\n"
      code += "      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it\n"
      code += "      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;\n"
      code += "      /// \\brief Give all system defined mappings that are to be implemented in C++ code for %s\n" % (escape(namespace_string))
      for s in self.declarations.sort_parameters(spec):
        code += "      /// \\param %s A sort expression\n" % (escape(str(s).lower()))
      code += "      /// \\return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for %s\n" % (escape(namespace_string))
      code += "      inline\n"
      code += "      implementation_map %s_cpp_implementable_mappings(%s)\n" % (namespace_string, sort_parameters)
      code += "      {\n"
      code += "        implementation_map result;\n"
      add_mappings_code = self.declarations.cplusplus_implementable_functions(spec, data_parameters)
      if add_mappings_code != "":
        code += add_mappings_code
      code += "        return result;\n"
      code += "      }\n"
      return code
    else:
      return ""

class constructor_specification():
  def __init__(self, declarations):
    assert(isinstance(declarations, function_declaration_list))
    self.declarations = declarations
    self.namespace = ""
    self.original_namespace = ""

  def has_cplusplus_implementable_code(self):
    return self.declarations.has_cplusplus_implementable_code()

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    self.declarations.set_namespace(string)

  def find_function(self, function, argumentcount):
    return self.declarations.find_function(function, argumentcount)

  def __str__(self):
    return "cons {0}".format(self.declarations)

  def merge_specification(self, spec, is_supertype):
    self.declarations.merge_declarations(spec.declarations, is_supertype)

  def code(self, spec):
    assert(isinstance(spec, specification))
    sort_parameters = ", ".join(["const sort_expression& %s" % (str(x).lower()) for x in self.declarations.sort_parameters(spec)])
    data_parameters = ", ".join(["%s" % (str(x).lower()) for x in self.declarations.sort_parameters(spec)])
    assert(self.namespace == spec.namespace)
    namespace_string = self.namespace
    if namespace_string == "undefined":
      namespace_string = spec.get_namespace()
    code  = self.declarations.code(spec)

    code += "      /// \\brief Give all system defined constructors for %s.\n" % (escape(namespace_string))
    for s in self.declarations.sort_parameters(spec):
      code += "      /// \\param %s A sort expression.\n" % (escape(str(s).lower()))
    code += "      /// \\return All system defined constructors for %s.\n" % (escape(namespace_string))
    code += "      inline\n"
    code += "      function_symbol_vector %s_generate_constructors_code(%s)\n" % (namespace_string,sort_parameters)
    code += "      {\n"
    code += "        function_symbol_vector result;\n"
    add_constructors_code = self.declarations.generator_code(spec) + (spec.sort_specification.structured_sort_constructor_code())
    if add_constructors_code != "":
      code += "%s" % (spec.sort_specification.structured_sort_constructor_code())
      code += "%s\n" % (self.declarations.generator_code(spec, True))
    code += "        return result;\n"
    code += "      }\n"

    code += "      /// \\brief Give all defined constructors which can be used in mCRL2 specs for %s.\n" % (escape(namespace_string))
    for s in self.declarations.sort_parameters(spec):
      code += "      /// \\param %s A sort expression.\n" % (escape(str(s).lower()))
    code += "      /// \\return All system defined constructors that can be used in an mCRL2 specification for %s.\n" % (escape(namespace_string))
    code += "      inline\n"
    code += "      function_symbol_vector %s_mCRL2_usable_constructors(%s)\n" % (namespace_string,sort_parameters)
    code += "      {\n"
    code += "        function_symbol_vector result;\n"
    add_constructors_code = self.declarations.mCRL2_usable_functions(spec) + (spec.sort_specification.structured_sort_constructor_code())
    if add_constructors_code != "":
      code += "%s" % (spec.sort_specification.structured_sort_constructor_code())
      code += "%s\n" % (self.declarations.mCRL2_usable_functions(spec, True))
    code += "        return result;\n"
    code += "      }\n"

    code += "      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it\n"
    code += "      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;\n"
    code += "      /// \\brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for %s.\n" % (escape(namespace_string))
    for s in self.declarations.sort_parameters(spec):
      code += "      /// \\param %s A sort expression.\n" % (escape(str(s).lower()))
    code += "      /// \\return All system defined constructors that are to be implemented in C++ for %s.\n" % (escape(namespace_string))
    code += "      inline\n"
    code += "      implementation_map %s_cpp_implementable_constructors(%s)\n" % (namespace_string,sort_parameters)
    code += "      {\n"
    code += "        implementation_map result;\n"
    add_constructors_code = self.declarations.cplusplus_implementable_functions(spec,data_parameters) + (spec.sort_specification.structured_sort_constructor_code())
    if add_constructors_code != "":
      code += "%s" % (spec.sort_specification.structured_sort_constructor_code())
      code += "%s\n" % (self.declarations.cplusplus_implementable_functions(spec, data_parameters, True))
    code += "        return result;\n"
    code += "      }\n"
    return code

class function_specification():
  def __init__(self, mapping_spec, constructor_spec):
    assert(isinstance(mapping_spec, mapping_specification))
    assert(isinstance(constructor_spec, constructor_specification))
    self.constructor_specification = constructor_spec
    self.mapping_specification = mapping_spec
    self.namespace = ""
    self.original_namespace = ""

  def has_cplusplus_implementable_code(self):
    return self.constructor_specification.has_cplusplus_implementable_code() or self.mapping_specification.has_cplusplus_implementable_code()

  def merge_specification(self, spec, is_supertype):
    self.constructor_specification.merge_specification(spec.constructor_specification, is_supertype)
    self.mapping_specification.merge_specification(spec.mapping_specification, is_supertype)

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    if self.constructor_specification != None:
      self.constructor_specification.set_namespace(string)
    self.mapping_specification.set_namespace(string)

  def find_function(self, function, argumentcount):
    result = self.constructor_specification.find_function(function, argumentcount)
    if not result:
      result = self.mapping_specification.find_function(function, argumentcount)
    return result

  def __str__(self):
    return "{0}\n{1}".format(self.constructor_specification, self.mapping_specification)

  def code(self, spec):
    assert(isinstance(spec, specification))
    # This assumes that no functions occur with the same name, but one is
    # constructor and one is mapping.
    code = ""
    code += self.constructor_specification.code(spec)
    code += self.mapping_specification.code(spec)
# We need deepcopy here to prevent undesired duplication
    merged_declarations = copy.deepcopy(self.constructor_specification.declarations)
    for e in self.mapping_specification.declarations.elements:
      merged_declarations.push_back(copy.deepcopy(e))
    code += merged_declarations.projection_code(spec)
    return code

class sort_specification():
  def __init__(self, declarations):
    assert(isinstance(declarations, sort_declaration_list))
    self.declarations = declarations
    self.namespace = ""
    self.original_namespace = ""

  def is_alias(self, sort_expr):
    assert(isinstance(sort_expr, sort_expression))
    return self.declarations.is_alias(sort_expr)

  def merge_structured_sorts(self, function_spec):
    return self.declarations.merge_structured_sorts(function_spec)

  def set_namespace(self, string):
    self.namespace = string
    if self.original_namespace == "":
      self.original_namespace = string
    self.declarations.set_namespace(string)

  def has_sort(self, sort):
    assert(isinstance(sort, sort_expression))
    for d in self.declarations.elements:
      if str(d.sort_expression) == str(sort):
        return True
    return False

  def get_sort(self, sort):
    assert(self.has_sort(sort))
    for d in self.declarations.elements:
      if str(d.sort_expression) == str(sort):
        return d

  def defines_container(self):
    return self.declarations.defines_container()

  def defines_struct(self):
    return self.declarations.defines_struct()

  def merge_specification(self, spec, is_supertype):
    self.declarations.merge_declarations(spec.declarations, is_supertype)

  def __str__(self):
    return "sort {0}".format(self.declarations)

  def structured_sort_constructor_code(self):
    return self.declarations.structured_sort_constructor_code()

  def structured_sort_mapping_code(self):
    return self.declarations.structured_sort_mapping_code()

  def structured_sort_equation_code(self):
    return self.declarations.structured_sort_equation_code()

  def code(self,spec):
    return self.declarations.code(spec)

class specification():
  def __init__(self, sort_spec, function_spec, variable_spec, equation_spec):
    assert(isinstance(sort_spec, sort_specification))
    assert(isinstance(function_spec, function_specification))
    assert(isinstance(variable_spec, variable_specification))
    assert(isinstance(equation_spec, equation_specification))
    self.sort_specification = sort_spec
    self.function_specification = function_spec
    self.variable_specification = variable_spec
    self.equation_specification = equation_spec
    self.includes = None
    self.uses = None
    self.subtypes = None
    self.namespace = ""
    self.original_namespace = ""
    self.originfile = ""

  def has_cplusplus_implementable_code(self):
    return self.function_specification.has_cplusplus_implementable_code()

  def set_includes(self, includes):
    assert(isinstance(includes, include_list))
    self.includes = includes

  def set_using(self, uses):
    assert(isinstance(uses, using_list))
    self.uses = uses

  def set_subtypes(self, subtypes):
    self.subtypes = subtypes

  def sort_parameters(self):
    if self.uses != None:
      return self.uses.sort_parameters()
    else:
      return []

  def has_lambda(self):
    return self.equation_specification.has_lambda()

  def has_forall(self):
    return self.equation_specification.has_forall()

  def has_exists(self):
    return self.equation_specification.has_exists()

  def defines_container(self):
    return self.sort_specification.defines_container()

  def defines_struct(self):
    return self.sort_specification.defines_struct()

  def set_origin_file(self, f):
    self.origin_file = f

  def get_namespace(self):
    if self.namespace != "undefined" and self.namespace != "":
      return self.namespace
    else:
      return self.origin_file.stem

  def set_namespace(self):
    if (len(self.sort_specification.declarations.elements) == 0):
      s = "undefined"
    else:
      s = "{0}".format(self.sort_specification.declarations.elements[0].label)
    self.namespace = s
    if self.original_namespace == "":
      self.original_namespace = s
    self.sort_specification.set_namespace(s)
    self.function_specification.set_namespace(s)
    self.equation_specification.set_namespace(s)

  def merge_specification(self, spec):
    assert(isinstance(spec, specification))
    self.sort_specification.merge_specification(spec.sort_specification, self.subtypes != None)
    self.function_specification.merge_specification(spec.function_specification, self.subtypes != None)

  def __str__(self):
    res = []
    if self.includes:
      res.append(str(self.includes))
    res.append(str(self.sort_specification))
    res.append(str(self.function_specification))
    res.append(str(self.variable_specification))
    res.append(str(self.equation_specification))
    return "\n".join(res)

  def code(self, infilename):
    # Add structured sorts to constructor declarations
    self.function_specification = self.sort_specification.merge_structured_sorts(self.function_specification)
    code  = ""
    code += "// Author(s): Jeroen Keiren\n"
    code += "// Copyright: see the accompanying file COPYING or copy at\n"
    code += "// https://github.com/mCRL2org/mCRL2/blob/master/COPYING\n"
    code += "//\n"
    code += "// Distributed under the Boost Software License, Version 1.0.\n"
    code += "// (See accompanying file LICENSE_1_0.txt or copy at\n"
    code += "// http://www.boost.org/LICENSE_1_0.txt)\n"
    code += "//\n"
    code += "/// \\file mcrl2/data/%s.h\n" % (remove_underscore(self.get_namespace()))
    code += "/// \\brief The standard sort %s.\n" % (self.get_namespace())
    code += "///\n"
    code += "/// This file was generated from the data sort specification\n"
    code += "/// mcrl2/data/build/%s.spec.\n" % (remove_underscore(self.get_namespace()))
    code += "\n"
    code += "#ifndef MCRL2_DATA_%s_H\n" % (infilename.stem.upper())
    code += "#define MCRL2_DATA_%s_H\n\n" % (infilename.stem.upper())
    code += "#include \"mcrl2/data/basic_sort.h\"\n"
    code += "#include \"mcrl2/data/function_sort.h\"\n"
    code += "#include \"mcrl2/data/function_symbol.h\"\n"
    code += "#include \"mcrl2/data/application.h\"\n"
    code += "#include \"mcrl2/data/data_equation.h\"\n"
    code += "#include \"mcrl2/data/standard.h\"\n"

    if self.has_lambda():
      code += "#include \"mcrl2/data/lambda.h\"\n"
    if self.has_forall():
      code += "#include \"mcrl2/data/forall.h\"\n"
    if self.has_exists():
      code += "#include \"mcrl2/data/exists.h\"\n"
    if self.defines_container():
      code += "#include \"mcrl2/data/container_sort.h\"\n"
    if self.defines_struct():
      code += "#include \"mcrl2/data/structured_sort.h\"\n"
    if self.includes != None:
      code += self.includes.code()

    code += "\n";
    code += "#include <functional>\n"
    code += "#include <map>\n"
    code += "\n";
    code += "#ifndef MCRL2_ENABLE_MODULES\n"
    code += "  #include \"mcrl2/utilities/exception.cxx\"\n"
    code += "#else\n"
    code += " import utilities;\n"
    code += "#endif\n"

    code += "\n"
    code += "namespace mcrl2 {\n\n"
    code += "  namespace data {\n\n"
    if self.namespace != "undefined":
      code += "    /// \\brief Namespace for system defined sort %s.\n" % (escape(self.get_namespace()))
      code += "    namespace sort_%s {\n\n" % (self.get_namespace())
    code += self.sort_specification.code(self)
    code += self.function_specification.code(self)
    code += self.equation_specification.code(self, self.function_specification, self.variable_specification)
    code += "\n"
    dependent_sorts = set([])
    auxiliary_sorts = set([])
    for e in self.sort_specification.declarations.elements:
      if str(e).startswith("@") and e.original_namespace == self.namespace:
        auxiliary_sorts.add("         specification.add_system_defined_sort(%s);\n" % (e.inline_code(self)))
    if self.namespace != "undefined":
      code += "    } // namespace sort_%s\n\n" % (self.get_namespace())
    code += "  } // namespace data\n\n"
    code += "} // namespace mcrl2\n\n"
    code += ("#include \"mcrl2/data/detail/%s.h\" // This file contains the manual implementations of rewrite functions.\n" % (remove_underscore(self.get_namespace()))
                 if self.has_cplusplus_implementable_code() else "")
    code += "#endif // MCRL2_DATA_%s_H\n" % (infilename.stem.upper())
    code = code.replace("__", "_")
    p = re.compile('sort_([A-Za-z0-9]*)_([ ]|:)')
    code = p.sub(r'sort_\1\2', code)
    p = re.compile('is_([A-Za-z0-9]*)_\(')
    code = p.sub(r'is_\1(',code)

    return code

class include_list():
  def __init__(self, elements):
    assert(all([isinstance(x, include) for x in elements]))
    self.elements = elements

  def push_back(self, element):
    assert(isinstance(element, include))
    self.elements += [element]

  def __str__(self):
    return "".join(["{0}\n".format(e) for e in self.elements])

  def code(self):
    s = ""
    for e in self.elements:
      s += "%s\n" % (e.code())
    return s

class include():
  def __init__(self, id):
    assert(isinstance(id, identifier))
    self.identifier = id

  def __str__(self):
    return "#include {0}".format(self.identifier)

  def code(self):
    return "#include \"mcrl2/data/{0}.h\"".format(str(self.identifier)[:-5])

class subtype_list():
  def __init__(self, elements):
    assert(all([isinstance(x, subtype) for x in elements]))
    self.elements = elements

  def push_back(self, element):
    assert(isinstance(element, subtype))
    self.elements += [element]

  def __str__(self):
    return "".join(["{0}\n".format(e) for e in self.elements])

  def code(self):
    s = ""
    for e in self.elements:
      s += "%s\n" % (e.code())
    return s

class subtype():
  def __init__(self, id):
    assert(isinstance(id, identifier))
    self.identifier = id

  def __str__(self):
    return "#supertypeof {0}".format(self.identifier)

class using_list():
  def __init__(self, elements):
    assert(all([isinstance(x, using) for x in elements]))
    self.elements = elements

  def push_back(self, element):
    assert(isinstance(element, using))
    self.elements += [element]

  def sort_parameters(self):
    return [x.sort_parameters() for x in self.elements]

  def __str__(self):
    return ", ".join(map(str, self.elements))

  def code(self):
    code = ""
    for e in self.elements:
      code += "%s\n" % (e.code())
    return code

class using():
  def __init__(self, identifier):
    self.identifier = identifier

  def __str__(self):
    return str(self.identifier)

  def sort_parameters(self):
    return sort_identifier(self.identifier)

  def code(self):
    return str(self)
