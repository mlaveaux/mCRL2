#!/usr/bin/env python3

import sys
from mcrl2_classes import *
from mcrl2_utility import *

MCRL2_ROOT = '../../'

class_map = mcrl2_class_map()
all_classes = parse_class_map(class_map)
modifiability_map = make_modifiability_map(all_classes)

file_map = {
  'action_formulas' : MCRL2_ROOT + 'libraries/modal_formula/source/modal_formula.cpp',
  'core' : MCRL2_ROOT + 'libraries/core/source/core.cpp',
  'data' : MCRL2_ROOT + 'libraries/data/source/data.cpp',
  'lps' : MCRL2_ROOT + 'libraries/lps/source/lps.cpp',
  'pbes_system' : MCRL2_ROOT + 'libraries/pbes/source/pbes.cpp',
  'process' : MCRL2_ROOT + 'libraries/process/source/process.cpp',
  'regular_formulas' : MCRL2_ROOT + 'libraries/modal_formula/source/modal_formula.cpp',
  'state_formulas' : MCRL2_ROOT + 'libraries/modal_formula/source/modal_formula.cpp',
}

PP_CLASSNAMES = '''
data::sort_expression_list
data::sort_expression_vector
data::data_expression_list
data::data_expression_vector
data::assignment_list
data::assignment_vector
data::variable_list
data::variable_vector
data::function_symbol_list
data::function_symbol_vector
data::structured_sort_constructor_list
data::structured_sort_constructor_vector
data::data_equation_list
data::data_equation_vector
pbes_system::pbes_equation_vector
pbes_system::pbes_expression_list
pbes_system::pbes_expression_vector
pbes_system::propositional_variable_list
pbes_system::propositional_variable_vector
pbes_system::propositional_variable_instantiation_list
pbes_system::propositional_variable_instantiation_vector
process::action_list
process::action_vector
process::action_label_list
process::action_label_vector
process::process_identifier_list
process::process_identifier_vector
process::process_expression_list
process::process_expression_vector
process::process_equation_list
process::process_equation_vector
'''

NORMALIZE_SORTS_CLASSNAMES = '''
data::data_equation
data::data_equation_list
data::data_equation_vector
data::data_expression
data::sort_expression
data::variable_list
lps::multi_action
lps::specification
lps::stochastic_specification
process::action
process::action_label_list
process::process_equation_vector
process::process_specification
pbes_system::pbes_equation_vector
pbes_system::pbes
pbes_system::pbes_expression
state_formulas::state_formula
'''

TRANSLATE_USER_NOTATION_CLASSNAMES = '''
data::data_expression
data::data_equation
lps::multi_action
pbes_system::pbes
pbes_system::pbes_expression
process::action
process::process_expression
process::process_specification
state_formulas::state_formula
'''

FIND_SORT_EXPRESSIONS_CLASSNAMES = '''
data::data_equation
data::data_expression
data::sort_expression
lps::specification
lps::stochastic_specification
pbes_system::pbes
process::action_label_list
process::process_equation_vector
process::process_expression
process::process_specification
state_formulas::state_formula
'''

FIND_VARIABLES_CLASSNAMES = '''
action_formulas::action_formula
data::data_expression
data::data_expression_list
data::function_symbol
data::variable
data::variable_list
lps::linear_process
lps::stochastic_linear_process
lps::specification
lps::stochastic_specification
lps::deadlock
lps::multi_action
pbes_system::pbes
process::action
state_formulas::state_formula
'''

FIND_FREE_VARIABLES_CLASSNAMES = '''
data::data_expression
data::data_expression_list
lps::linear_process
lps::stochastic_linear_process
lps::specification
lps::stochastic_specification
lps::deadlock
lps::multi_action
lps::process_initializer
lps::stochastic_process_initializer
pbes_system::pbes
pbes_system::pbes_expression
pbes_system::pbes_equation
process::action
process::process_specification
state_formulas::state_formula
'''

FIND_FUNCTION_SYMBOLS_CLASSNAMES = '''
data::data_equation
lps::specification
lps::stochastic_specification
pbes_system::pbes
'''

FIND_PROPOSITIONAL_VARIABLE_INSTANTIATIONS_CLASSNAMES = '''
pbes_system::pbes_expression
'''

FIND_IDENTIFIERS_CLASSNAMES = '''
data::variable_list
lps::specification
lps::stochastic_specification
process::process_specification
pbes_system::pbes_expression
state_formulas::state_formula
'''

FIND_ACTION_LABELS_CLASSNAMES = '''
lps::linear_process
lps::process_initializer
lps::specification
lps::stochastic_specification
state_formulas::state_formula
'''

SEARCH_VARIABLE_CLASSNAMES = '''
data::data_expression
pbes_system::pbes_expression
'''

def has_specification(type):
    return type.endswith('specification') or type.endswith('pbes')

def is_modifiable(type):
    if type in modifiability_map:
        return modifiability_map[type]
    elif type.endswith('_list'):
        return False
    elif type.endswith('_vector'):
        return True
    elif type.endswith('pbes'):
        return True
    elif type.endswith('vector<pbes_equation>'):
        return True
    raise Exception('Unknown type %s!' % type)

def extract_namespace(classname):
    return re.sub('::.*', '', classname)

def generate_traverser_overloads(classnames, function, return_type, parameters: list[str], code_map):
    for classname in classnames:
        namespace = extract_namespace(classname)
        params = [(f'const {classname}&', 'x')]
        params.extend([(p, f'arg{i}') for i, p in enumerate(parameters)])

        text = re.sub('>>', '> >', '%s %s(%s) { return %s::%s< %s >(%s); }\n' % (return_type, function, ", ".join(map(lambda x : x[0] + " " + x[1], params)), namespace, function, classname, 
                                                                                 ", ".join(map(lambda x : x[1], params))))
        if namespace in code_map:
            code_map[namespace].append(text)

def generate_builder_overloads(classnames, function, code_map):
    for classname in classnames:
        namespace = extract_namespace(classname)
        if is_modifiable(classname):
            text = 'void %s(%s& x) { %s::%s< %s >(x); }\n' % (function, classname, namespace, function, classname)
        else:
            text = '%s %s(const %s& x) { return %s::%s< %s >(x); }\n' % (classname, function, classname, namespace, function, classname)
        if namespace in code_map:
            code_map[namespace].append(text)

# special because of additional variable argument
def generate_search_variable_overloads(classnames, function, return_type, code_map):
    for classname in classnames:
        namespace = extract_namespace(classname)
        text = re.sub('>>', '> >', '%s %s(const %s& x, const data::variable& v) { return %s::%s< %s >(x, v); }\n' % (return_type, function, classname, namespace, function, classname))
        if namespace in code_map:
            code_map[namespace].append(text)

# special because of additional data_specification argument
def generate_normalize_sorts_overloads(classnames, code_map):
    for classname in classnames:
        namespace = extract_namespace(classname)
        if is_modifiable(classname):
            text = 'void normalize_sorts(%s& x, const data::sort_specification& sortspec) { %s::normalize_sorts< %s >(x, sortspec); }\n' % (classname, namespace, classname)
        else:
            text = '%s normalize_sorts(const %s& x, const data::sort_specification& sortspec) { return %s::normalize_sorts< %s >(x, sortspec); }\n' % (classname, classname, namespace, classname)
        if has_specification(classname):
            text = re.sub('x, sortspec', 'x, x.data()', text)
            text = re.sub('& sortspec', '& /* sortspec */', text)
        if namespace in code_map:
            code_map[namespace].append(text)

code_map = {}
for namespace in file_map:
    code_map[namespace] = []

# add pp overloads for all known classes
for name in sorted(all_classes):
    c = all_classes[name]
    PP_CLASSNAMES = PP_CLASSNAMES + '\n%s::%s' % (c.namespace(), c.classname())

classnames = PP_CLASSNAMES.strip().split()
generate_traverser_overloads(classnames, 'pp', 'std::string', ['bool'], code_map)

classnames = NORMALIZE_SORTS_CLASSNAMES.strip().split()
generate_normalize_sorts_overloads(classnames, code_map)

classnames = TRANSLATE_USER_NOTATION_CLASSNAMES.strip().split()
generate_builder_overloads(classnames, 'translate_user_notation', code_map)

classnames = FIND_SORT_EXPRESSIONS_CLASSNAMES.strip().split()
generate_traverser_overloads(classnames, 'find_sort_expressions', 'std::set<data::sort_expression>', [], code_map)

classnames = FIND_VARIABLES_CLASSNAMES.strip().split()
generate_traverser_overloads(classnames, 'find_all_variables', 'std::set<data::variable>', [], code_map)

classnames = FIND_FREE_VARIABLES_CLASSNAMES.strip().split()
generate_traverser_overloads(classnames, 'find_free_variables', 'std::set<data::variable>', [], code_map)

classnames = FIND_FUNCTION_SYMBOLS_CLASSNAMES.strip().split()
generate_traverser_overloads(classnames, 'find_function_symbols', 'std::set<data::function_symbol>', [], code_map)

classnames = FIND_PROPOSITIONAL_VARIABLE_INSTANTIATIONS_CLASSNAMES.strip().split()
generate_traverser_overloads(classnames, 'find_propositional_variable_instantiations', 'std::set<pbes_system::propositional_variable_instantiation>', [], code_map)

classnames = FIND_IDENTIFIERS_CLASSNAMES.strip().split()
generate_traverser_overloads(classnames, 'find_identifiers', 'std::set<core::identifier_string>', [], code_map)

classnames = FIND_ACTION_LABELS_CLASSNAMES.strip().split()
generate_traverser_overloads(classnames, 'find_action_labels', 'std::set<process::action_label>', [], code_map)

classnames = SEARCH_VARIABLE_CLASSNAMES.strip().split()
generate_search_variable_overloads(classnames, 'search_variable', 'bool', code_map)

result = True
for namespace in code_map:
    filename = file_map[namespace]
    text = ''.join(code_map[namespace])
    label = 'generated %s overloads' % namespace
    result = insert_text_in_file(filename, text, label) and result
sys.exit(not result) # 0 result indicates successful execution
