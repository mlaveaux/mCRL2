#!/usr/bin/env python3

#~ Copyright 2009, 2010 Wieger Wesselink.
#~ Distributed under the Boost Software License, Version 1.0.
#~ (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

# This file contains tables that are used to generate classes and traversal functions
# for these classes. A prerequisite is that each class has a corresponding aterm
# representation (the first column of each table). The second column contains the
# constructor of the classes.

import re
import string

# Modifiers in the second column of class specifications:
#
# C = generate container typedefs (vector + list)
# E = it is an expression derived class
# I = generate is-function
# M = the class is modifiable
# O = generate constructor overloads
# s = skip constructor generation
# S = skip class generation
# U = generate user section
# X = it is an expression super class
# W = do not generate swap overload
# N = term has an additional index as last argument  TODO: Constructors must still be adapted. 
# i = this class can be cast to an aterm_int.

CORE_CLASSES = r'''
identifier_string() : public atermpp::aterm_string | SC | String | An identifier
'''

CONTAINER_TYPES = r'''
container_type() : public atermpp::aterm       | XCO | SortConsType | Container type
list_container() : public data::container_type | EIO | SortList     | Container type for lists
set_container()  : public data::container_type | EIO | SortSet      | Container type for sets
bag_container()  : public data::container_type | EIO | SortBag      | Container type for bags
fset_container() : public data::container_type | EIO | SortFSet     | Container type for finite sets
fbag_container() : public data::container_type | EIO | SortFBag     | Container type for finite bags
'''

STRUCTURED_SORT_ELEMENTS = r'''
structured_sort_constructor_argument(const core::identifier_string& name, const sort_expression& sort)                                                            : public atermpp::aterm | ICUO | StructProj | An argument of a constructor of a structured sort
structured_sort_constructor(const core::identifier_string& name, const structured_sort_constructor_argument_list& arguments, core::identifier_string& recogniser) : public atermpp::aterm | ICUO | StructCons | A constructor for a structured sort
'''

SORT_EXPRESSION_CLASSES = r'''
sort_expression()                                                                         : public atermpp::aterm        | XOCU | SortExpr             | A sort expression
basic_sort(const core::identifier_string& name)                                           : public data::sort_expression | EO   | SortId               | A basic sort
container_sort(const container_type& container_name, const sort_expression& element_sort) : public data::sort_expression | EO   | SortCons             | A container sort
structured_sort(const structured_sort_constructor_list& constructors)                     : public data::sort_expression | CEOU | SortStruct           | A structured sort
function_sort(const sort_expression_list& domain, const sort_expression& codomain)        : public data::sort_expression | EO   | SortArrow            | A function sort
untyped_sort()                                                                            : public data::sort_expression | EO   | UntypedSortUnknown   | Unknown sort expression
untyped_possible_sorts(const sort_expression_list& sorts)                                 : public data::sort_expression | EO   | UntypedSortsPossible | Multiple possible sorts
untyped_sort_variable(const atermpp::aterm_int& value)                                    : public data::sort_expression | EIOU | UntypedSortVariable  | Untyped sort variable
'''

BINDER_TYPES = r'''
binder_type()                             : public atermpp::aterm       | XCO | BindingOperator   | Binder
untyped_set_or_bag_comprehension_binder() : public data::binder_type    | EIO | UntypedSetBagComp | Binder for untyped set or bag comprehension
set_comprehension_binder()                : public data::binder_type    | EIO | SetComp           | Binder for set comprehension
bag_comprehension_binder()                : public data::binder_type    | EIO | BagComp           | Binder for bag comprehension
forall_binder()                           : public data::binder_type    | EIO | Forall            | Binder for universal quantification
exists_binder()                           : public data::binder_type    | EIO | Exists            | Binder for existential quantification
lambda_binder()                           : public data::binder_type    | EIO | Lambda            | Binder for lambda abstraction
'''

ASSIGNMENT_EXPRESSION_CLASSES = r'''
assignment_expression()                                                               : public atermpp::aterm                      | XIOC  | WhrDecl                     | Assignment expression
assignment(const variable& lhs, const data_expression& rhs)                           : public data::assignment_expression         | EIOUC | DataVarIdInit               | Assignment of a data expression to a variable
untyped_identifier_assignment(const core::identifier_string& lhs, const data_expression& rhs) : public data::assignment_expression | EIOUC | UntypedIdentifierAssignment | Assignment of a data expression to a string
'''

DATA_EXPRESSION_CLASSES = r'''
data_expression()                                                                                             : public atermpp::aterm        | XCU    | DataExpr          | A data expression
variable(const core::identifier_string& name, const sort_expression& sort)                                    : public data::data_expression | EOCUs  | DataVarId         | A data variable
function_symbol(const core::identifier_string& name, const sort_expression& sort)                             : public data::data_expression | EOCUsN | OpId              | A function symbol
application(const data_expression& head, data_expression_list const& arguments)                               : public data::data_expression | EOUSW  | DataAppl          | An application of a data expression to a number of arguments
where_clause(const data_expression& body, const assignment_expression_list& declarations)                     : public data::data_expression | EOU    | Whr               | A where expression
abstraction(const binder_type& binding_operator, const variable_list& variables, const data_expression& body) : public data::data_expression | EO     | Binder            | An abstraction expression
machine_number(std::size_t value)                                                                             : public data::data_expression | EOCi   | Number            | A machine number
untyped_identifier(const core::identifier_string& name)                                                       : public data::data_expression | EO     | UntypedIdentifier | An untyped identifier
'''

ABSTRACTION_EXPRESSION_CLASSES = r'''
abstraction()                                                                                 : public data::data_expression | XES | None | Abstraction
forall(const variable_list& variables, const data_expression& body)                           : public data::abstraction     | ES  | None | Universal quantification
exists(const variable_list& variables, const data_expression& body)                           : public data::abstraction     | ES  | None | Existential quantification
lambda(const variable_list& variables, const data_expression& body)                           : public data::abstraction     | ES  | None | Lambda abstraction
set_comprehension(const variable_list& variables, const data_expression& body)                : public data::abstraction     | ES  | None | Set comprehension
bag_comprehension(const variable_list& variables, const data_expression& body)                : public data::abstraction     | ES  | None | Bag comprehension
untyped_set_or_bag_comprehension(const variable_list& variables, const data_expression& body) : public data::abstraction     | ES  | None | Untyped Set/Bag comprehension
'''

DATA_CLASSES = r'''
alias(const basic_sort& name, const sort_expression& reference)                                                                         : public atermpp::aterm | CIO | SortRef              | A sort alias
data_equation(const variable_list& variables, const data_expression& condition, const data_expression& lhs, const data_expression& rhs) : public atermpp::aterm | COU | DataEqn              | A data equation
data_specification()                                                                                                                    : public atermpp::aterm | SW  | DataSpec             | A data specification
untyped_data_parameter(const core::identifier_string& name, const data_expression_list& arguments)                                      : public atermpp::aterm | CI  | UntypedDataParameter | An untyped parameter
'''

STATE_FORMULA_CLASSES = r'''
state_formula()                                                                                                 : public atermpp::aterm                | XICU | StateFrm                 | A state formula
true_()                                                                                                         : public state_formulas::state_formula | EI   | StateTrue                | The value true for state formulas
false_()                                                                                                        : public state_formulas::state_formula | EI   | StateFalse               | The value false for state formulas
not_(const state_formula& operand)                                                                              : public state_formulas::state_formula | EI   | StateNot                 | The not operator for state formulas
minus(const state_formula& operand)                                                                             : public state_formulas::state_formula | EI   | StateMinus               | The minus operator for state formulas
and_(const state_formula& left, const state_formula& right)                                                     : public state_formulas::state_formula | EI   | StateAnd                 | The and operator for state formulas
or_(const state_formula& left, const state_formula& right)                                                      : public state_formulas::state_formula | EI   | StateOr                  | The or operator for state formulas
imp(const state_formula& left, const state_formula& right)                                                      : public state_formulas::state_formula | EI   | StateImp                 | The implication operator for state formulas
plus(const state_formula& left, const state_formula& right)                                                     : public state_formulas::state_formula | EI   | StatePlus                | The plus operator for state formulas with values
const_multiply(const data::data_expression& left, const state_formula& right)                                   : public state_formulas::state_formula | EI   | StateConstantMultiply    | The multiply operator for state formulas with values
const_multiply_alt(const state_formula& left, const data::data_expression& right)                               : public state_formulas::state_formula | EI   | StateConstantMultiplyAlt | The multiply operator for state formulas with values
forall(const data::variable_list& variables, const state_formula& body)                                         : public state_formulas::state_formula | EI   | StateForall              | The universal quantification operator for state formulas
exists(const data::variable_list& variables, const state_formula& body)                                         : public state_formulas::state_formula | EI   | StateExists              | The existential quantification operator for state formulas
infimum(const data::variable_list& variables, const state_formula& body)                                        : public state_formulas::state_formula | EI   | StateInfimum             | The infimum over a data type for state formulas
supremum(const data::variable_list& variables, const state_formula& body)                                       : public state_formulas::state_formula | EI   | StateSupremum            | The supremum over a data type for state formulas
sum(const data::variable_list& variables, const state_formula& body)                                            : public state_formulas::state_formula | EI   | StateSum                 | The sum over a data type for state formulas
must(const regular_formulas::regular_formula& formula, const state_formula& operand)                            : public state_formulas::state_formula | EI   | StateMust                | The must operator for state formulas
may(const regular_formulas::regular_formula& formula, const state_formula& operand)                             : public state_formulas::state_formula | EI   | StateMay                 | The may operator for state formulas
yaled()                                                                                                         : public state_formulas::state_formula | EI   | StateYaled               | The yaled operator for state formulas
yaled_timed(const data::data_expression& time_stamp)                                                            : public state_formulas::state_formula | EI   | StateYaledTimed          | The timed yaled operator for state formulas
delay()                                                                                                         : public state_formulas::state_formula | EI   | StateDelay               | The delay operator for state formulas
delay_timed(const data::data_expression& time_stamp)                                                            : public state_formulas::state_formula | EI   | StateDelayTimed          | The timed delay operator for state formulas
variable(const core::identifier_string& name, const data::data_expression_list& arguments)                      : public state_formulas::state_formula | EI   | StateVar                 | The state formula variable
nu(const core::identifier_string& name, const data::assignment_list& assignments, const state_formula& operand) : public state_formulas::state_formula | EI   | StateNu                  | The nu operator for state formulas
mu(const core::identifier_string& name, const data::assignment_list& assignments, const state_formula& operand) : public state_formulas::state_formula | EI   | StateMu                  | The mu operator for state formulas
'''

REGULAR_FORMULA_CLASSES = r'''
regular_formula()                                                                                                       : public atermpp::aterm                    | XIC | RegFrm        | A regular formula
seq(const regular_formula& left, const regular_formula& right)                                                          : public regular_formulas::regular_formula | EI  | RegSeq        | The seq operator for regular formulas
alt(const regular_formula& left, const regular_formula& right)                                                          : public regular_formulas::regular_formula | EI  | RegAlt        | The alt operator for regular formulas
trans(const regular_formula& operand)                                                                                   : public regular_formulas::regular_formula | EI  | RegTrans      | The trans operator for regular formulas
trans_or_nil(const regular_formula& operand)                                                                            : public regular_formulas::regular_formula | EI  | RegTransOrNil | The 'trans or nil' operator for regular formulas
untyped_regular_formula(const core::identifier_string& name, const regular_formula& left, const regular_formula& right) : public regular_formulas::regular_formula | EI  | UntypedRegFrm | An untyped regular formula or action formula
'''

ACTION_FORMULA_CLASSES = r'''
action_formula()                                                            : public atermpp::aterm                  | XIC  | ActFrm            | An action formula
true_()                                                                     : public action_formulas::action_formula | EI   | ActTrue           | The value true for action formulas
false_()                                                                    : public action_formulas::action_formula | EI   | ActFalse          | The value false for action formulas
not_(const action_formula& operand)                                         : public action_formulas::action_formula | EI   | ActNot            | The not operator for action formulas
and_(const action_formula& left, const action_formula& right)               : public action_formulas::action_formula | EI   | ActAnd            | The and operator for action formulas
or_(const action_formula& left, const action_formula& right)                : public action_formulas::action_formula | EI   | ActOr             | The or operator for action formulas
imp(const action_formula& left, const action_formula& right)                : public action_formulas::action_formula | EI   | ActImp            | The implication operator for action formulas
forall(const data::variable_list& variables, const action_formula& body)    : public action_formulas::action_formula | EI   | ActForall         | The universal quantification operator for action formulas
exists(const data::variable_list& variables, const action_formula& body)    : public action_formulas::action_formula | EI   | ActExists         | The existential quantification operator for action formulas
at(const action_formula& operand, const data::data_expression& time_stamp)  : public action_formulas::action_formula | EI   | ActAt             | The at operator for action formulas
multi_action(const process::action_list& actions)                           : public action_formulas::action_formula | EI   | ActMultAct        | The multi action for action formulas
'''

MODAL_FORMULA_CLASSES = r'''
state_formula_specification(const data::data_specification& data, const process::action_label_list& action_labels, const state_formula& formula) | SMW | StateFrmSpec    | A state formula specification
'''

# N.B. This one is problematic due to the optional time in deadlock/multi_action.
LPS_CLASSES = r'''
deadlock(const data::data_expression& time)                                                                                                                                                                                                                   | CMS  | None              | A deadlock
multi_action(const process::action_list& actions, const data::data_expression& time)                                                                                                                                             : public atermpp::aterm      | CIUs | TimedMultAct      | A timed multi-action
deadlock_summand(const data::variable_list& summation_variables, const data::data_expression& condition, const lps::deadlock& deadlock)                                                                                                                       | CMS  | None              | A deadlock summand
action_summand(const data::variable_list& summation_variables, const data::data_expression& condition, const lps::multi_action& multi_action, const data::assignment_list& assignments)                                                                       | CMS  | None              | An action summand
process_initializer(const data::data_expression_list& expressions)                                                                                                                                                               : public atermpp::aterm      | CIUS | LinearProcessInit | A process initializer
linear_process(const data::variable_list& process_parameters, const deadlock_summand_vector& deadlock_summands, const action_summand_vector& action_summands)                                                                                                 | MSW  | LinearProcess     | A linear process
specification(const data::data_specification& data, const process::action_label_list& action_labels, const std::set<data::variable>& global_variables,const linear_process& process, const process_initializer& initial_process)                              | MSW  | LinProcSpec       | A linear process specification
stochastic_distribution(const data::variable_list& variables, const data::data_expression& distribution)                                                                                                                         : public atermpp::aterm      | CIU  | Distribution | A stochastic distribution
stochastic_action_summand(const data::variable_list& summation_variables, const data::data_expression& condition, const lps::multi_action& multi_action, const data::assignment_list& assignments, const stochastic_distribution& distribution) : public lps::action_summand | CMS  | None              | A stochastic action summand
stochastic_linear_process(const data::variable_list& process_parameters, const deadlock_summand_vector& deadlock_summands, const stochastic_action_summand_vector& action_summands) : public linear_process      | MSW  | LinearProcess     | A stochastic linear process
stochastic_specification(const data::data_specification& data, const process::action_label_list& action_labels, const std::set<data::variable>& global_variables, const stochastic_linear_process& process, const process_initializer& initial_process) : public lps::specification       | MSW  | LinProcSpec       | A stochastic linear process specification
stochastic_process_initializer(const data::data_expression_list& expressions, const stochastic_distribution& distribution)                                                                                                  : public lps::process_initializer | CIS  | LinearProcessInit | A stochastic process initializer
'''

PROCESS_CLASSES = r'''
action_label(const core::identifier_string& name, const data::sort_expression_list& sorts)                                                                                       : public atermpp::aterm | CI   | ActId              | An action label
process_specification(const data::data_specification& data, const process::action_label_list& action_labels, const std::set<data::variable>& global_variables, const std::vector<process::process_equation>& equations, const process_expression& init)           | SMW | ProcSpec    | A process specification
process_identifier(const core::identifier_string& name, const data::variable_list& variables)                                                                                    : public atermpp::aterm | CIUs | ProcVarId          | A process identifier
process_equation(const process_identifier& identifier, const data::variable_list& formal_parameters, const process_expression& expression)                                       : public atermpp::aterm | CI   | ProcEqn            | A process equation
rename_expression(core::identifier_string& source, core::identifier_string& target)                                                                                              : public atermpp::aterm | CI   | RenameExpr         | A rename expression
communication_expression(const action_name_multiset& action_name, const core::identifier_string& name)                                                                           : public atermpp::aterm | CI   | CommExpr           | A communication expression
action_name_multiset(const core::identifier_string_list& names)                                                                                                                  : public atermpp::aterm | CI   | MultActName        | A multiset of action names
untyped_multi_action(const data::untyped_data_parameter_list& actions)                                                                                                           : public atermpp::aterm | CI   | UntypedMultiAction | An untyped multi action or data application
'''

PROCESS_EXPRESSION_CLASSES = r'''
process_expression()                                                                                                                    : public atermpp::aterm              | XIC | ProcExpr                 | A process expression
action(const action_label& label, const data::data_expression_list& arguments)                                                          : public process::process_expression | EIC | Action                   | An action
process_instance(const process_identifier& identifier, const data::data_expression_list& actual_parameters)                             : public process::process_expression | EI  | Process                  | A process
process_instance_assignment(const process_identifier& identifier, const data::assignment_list& assignments)                             : public process::process_expression | EI  | ProcessAssignment        | A process assignment
delta()                                                                                                                                 : public process::process_expression | EI  | Delta                    | The value delta
tau()                                                                                                                                   : public process::process_expression | EI  | Tau                      | The value tau
sum(const data::variable_list& variables, const process_expression& operand)                                                            : public process::process_expression | EI  | Sum                      | The sum operator
block(const core::identifier_string_list& block_set, const process_expression& operand)                                                 : public process::process_expression | EI  | Block                    | The block operator
hide(const core::identifier_string_list& hide_set, const process_expression& operand)                                                   : public process::process_expression | EI  | Hide                     | The hide operator
rename(const rename_expression_list& rename_set, const process_expression& operand)                                                     : public process::process_expression | EI  | Rename                   | The rename operator
comm(const communication_expression_list& comm_set, const process_expression& operand)                                                  : public process::process_expression | EI  | Comm                     | The communication operator
allow(const action_name_multiset_list& allow_set, const process_expression& operand)                                                    : public process::process_expression | EI  | Allow                    | The allow operator
sync(const process_expression& left, const process_expression& right)                                                                   : public process::process_expression | EI  | Sync                     | The synchronization operator
at(const process_expression& operand, const data::data_expression& time_stamp)                                                          : public process::process_expression | EI  | AtTime                   | The at operator
seq(const process_expression& left, const process_expression& right)                                                                    : public process::process_expression | EI  | Seq                      | The sequential composition
if_then(const data::data_expression& condition, const process_expression& then_case)                                                    : public process::process_expression | EI  | IfThen                   | The if-then operator
if_then_else(const data::data_expression& condition, const process_expression& then_case, const process_expression& else_case)          : public process::process_expression | EI  | IfThenElse               | The if-then-else operator
bounded_init(const process_expression& left, const process_expression& right)                                                           : public process::process_expression | EI  | BInit                    | The bounded initialization
merge(const process_expression& left, const process_expression& right)                                                                  : public process::process_expression | EI  | Merge                    | The merge operator
left_merge(const process_expression& left, const process_expression& right)                                                             : public process::process_expression | EI  | LMerge                   | The left merge operator
choice(const process_expression& left, const process_expression& right)                                                                 : public process::process_expression | EI  | Choice                   | The choice operator
stochastic_operator(const data::variable_list& variables, const data::data_expression& distribution, const process_expression& operand) : public process::process_expression | EI  | StochasticOperator       | The distribution operator
untyped_process_assignment(const core::identifier_string& name, const data::untyped_identifier_assignment_list& assignments)            : public process::process_expression | EI  | UntypedProcessAssignment | An untyped process assginment
'''

PBES_CLASSES = r'''
fixpoint_symbol()                                                                                                                                                                                                            : public atermpp::aterm | XCU | FixPoint    | A fixpoint symbol
propositional_variable(const core::identifier_string& name, const data::variable_list& parameters)                                                                                                                           : public atermpp::aterm | CIU | PropVarDecl | A propositional variable declaration
pbes_equation(const fixpoint_symbol& symbol, const propositional_variable& variable, const pbes_expression& formula)                                                                                                                                 | SMC | PBEqn       | A PBES equation
pbes(const data::data_specification& data, const std::set<data::variable>& global_variables, const std::vector<pbes_system::pbes_equation>& equations, const propositional_variable_instantiation& initial_state)                                    | SMW | PBES        | A PBES
'''

PBES_EXPRESSION_CLASSES = r'''
pbes_expression()                                                                                                       : public atermpp::aterm               | XCI   | PBExpr            | A pbes expression
propositional_variable_instantiation(const core::identifier_string& name, const data::data_expression_list& parameters) : public pbes_system::pbes_expression | ECIUs | PropVarInst       | A propositional variable instantiation
not_(const pbes_expression& operand)                                                                                    : public pbes_system::pbes_expression | EI    | PBESNot           | The not operator for pbes expressions
and_(const pbes_expression& left, const pbes_expression& right)                                                         : public pbes_system::pbes_expression | EI    | PBESAnd           | The and operator for pbes expressions
or_(const pbes_expression& left, const pbes_expression& right)                                                          : public pbes_system::pbes_expression | EI    | PBESOr            | The or operator for pbes expressions
imp(const pbes_expression& left, const pbes_expression& right)                                                          : public pbes_system::pbes_expression | EI    | PBESImp           | The implication operator for pbes expressions
forall(const data::variable_list& variables, const pbes_expression& body)                                               : public pbes_system::pbes_expression | EI    | PBESForall        | The universal quantification operator for pbes expressions
exists(const data::variable_list& variables, const pbes_expression& body)                                               : public pbes_system::pbes_expression | EI    | PBESExists        | The existential quantification operator for pbes expressions
'''

PRES_CLASSES = r'''
pres_equation(const fixpoint_symbol& symbol, const propositional_variable& variable, const pres_expression& formula)                                                                                                                                      | SMC | PREqn       | A PRES equation
pres(const data::data_specification& data, const std::set<data::variable>& global_variables, const std::vector<pres_system::pres_equation>& equations, const propositional_variable_instantiation& initial_state)                                         | SMW | PRES        | A PRES
'''

PRES_EXPRESSION_CLASSES = r'''
pres_expression()                                                                                                       : public atermpp::aterm               | XCI   | PRExpr                  | A pres expression
propositional_variable_instantiation(const core::identifier_string& name, const data::data_expression_list& parameters) : public pres_system::pres_expression | ECIUs | PropVarInst             | A propositional variable instantiation
minus(const pres_expression& operand)                                                                                   : public pres_system::pres_expression | EI    | PRESMinus               | The not operator for pres expressions
and_(const pres_expression& left, const pres_expression& right)                                                         : public pres_system::pres_expression | EI    | PRESAnd                 | The and operator for pres expressions
or_(const pres_expression& left, const pres_expression& right)                                                          : public pres_system::pres_expression | EI    | PRESOr                  | The or operator for pres expressions
imp(const pres_expression& left, const pres_expression& right)                                                          : public pres_system::pres_expression | EI    | PRESImp                 | The implication operator for pres expressions
plus(const pres_expression& left, const pres_expression& right)                                                         : public pres_system::pres_expression | EI    | PRESPlus                | The addition operator for pres expressions
const_multiply(const data::data_expression& left, const pres_expression& right)                                         : public pres_system::pres_expression | EI    | PRESConstantMultiply    | The multiplication with a positive constant with the constant at the left.
const_multiply_alt(const pres_expression& left, const data::data_expression& right)                                     : public pres_system::pres_expression | EI    | PRESConstantMultiplyAlt | The multiplication with a positive constant with the constant at the right.
infimum(const data::variable_list& variables, const pres_expression& body)                                              : public pres_system::pres_expression | EI    | PRESInfimum             | The infimum over a data type for pres expressions
supremum(const data::variable_list& variables, const pres_expression& body)                                             : public pres_system::pres_expression | EI    | PRESSupremum            | The supremeum over a data type for pres expressions
sum(const data::variable_list& variables, const pres_expression& body)                                                  : public pres_system::pres_expression | EI    | PRESSum                 | The generic sum operator for pres expressions
eqinf(const pres_expression& operand)                                                                                   : public pres_system::pres_expression | EI    | PRESEqInf               | The indicator whether the argument is infinite
eqninf(const pres_expression& operand)                                                                                  : public pres_system::pres_expression | EI    | PRESEqNInf              | The indicator whether the argument is -infinite
condsm(const pres_expression& arg1, const pres_expression& arg2, const pres_expression& arg3)                           : public pres_system::pres_expression | EI    | PRESCondSm              | Conditional operator with condition smaller than 0
condeq(const pres_expression& arg1, const pres_expression& arg2, const pres_expression& arg3)                           : public pres_system::pres_expression | EI    | PRESCondEq              | Conditional operator with condition smaller than or equal to 0
'''

ADDITIONAL_EXPRESSION_CLASS_DEPENDENCIES = {
  'state_formulas::state_formula'     : [ 'data::data_expression', 'data::untyped_data_parameter' ],
  'action_formulas::action_formula'   : [ 'data::data_expression', 'data::untyped_data_parameter', 'process::untyped_multi_action' ],
  'regular_formulas::regular_formula' : [ 'action_formulas::action_formula', 'data::data_expression' ],
  'pbes_system::pbes_expression'      : [ 'data::data_expression', 'data::variable', 'data::untyped_data_parameter' ],
  'pres_system::pres_expression'      : [ 'data::data_expression', 'data::variable', 'data::untyped_data_parameter' ],
  'process::process_expression'       : [ 'data::untyped_data_parameter' ],
}

# removes 'const' and '&' from a type
def extract_type(text):
    text = re.sub(r'\bconst\b', '', text)
    text = re.sub(r'\s*&$', '', text)
    text = text.strip()
    return text

#indents the text with the given prefix
def indent_text(text, indent):
    lines = []
    for line in text.split('\n'):
        lines.append(indent + line)
    return '\n'.join(lines)

# Represents a function parameter like the following;
#
# const core::identifier_string& name
#
# self.modifiers: ['const', '&', '*'] etc.
# self.aterm: The corresponding aterm (default: None)
class Parameter:
    def __init__(self, name, type, modifiers, namespace):
        self.name_ = name
        self.type_ = type
        self.modifiers_ = modifiers
        self.aterm_ = None
        self.namespace_ = namespace

    # Returns the namespace of the parameter
    #
    def namespace(self):
        return self.namespace_

    # Returns true if the type is a list
    #
    # TODO: make this test accurate using the aterm grammar
    def is_list(self):
        return self.type(False).endswith('_list')

    # Returns the type of the parameter
    #
    # 'const core::identifier_string&'
    def type(self, include_modifiers = True, include_namespace = False, remove_templates = False):
        type1 = self.type_
        if remove_templates:
            pos1 = type1.find('<')
            pos2 = type1.rfind('>')
            if pos1 != -1 and pos2 > pos1:
                type1 = type1[pos1+1:pos2]
        if include_namespace and extract_namespace(type1) == None:
            type1 = '%s::%s' % (self.namespace(), type1)
        if include_modifiers:
            if 'const' in self.modifiers_:
                prefix = 'const '
            else:
                prefix = ''
            if '*' in self.modifiers_:
                postfix = '*'
            elif '&' in self.modifiers_:
                postfix = '&'
            else:
                postfix = ''
            return '%s%s%s' % (prefix, type1, postfix)
        return type1

    # Returns the name of the parameter
    #
    def name(self):
        return self.name_

    # Returns the corresponding aterm of the parameter, or None if no such term exists
    #
    def aterm(self):
        return self.aterm_

    def __repr__(self):
        return '%s %s' % (self.type(), self.name())

# Represents a function declaration like the following;
#
# variable(const core::identifier_string& name, const data::data_expression_list& arguments)
class FunctionDeclaration:
    def __init__(self, text, namespace):
        self.text = text.strip()
        self.namespace_ = namespace

        # compute argument text
        text = self.text
        text = re.sub('.*\(', '', text)
        text = re.sub('\).*', '', text)
        self.argument_text_ = text

        # compute parameters
        parameters = []
        words = [word.strip() for word in text.split(',')]
        for word in words:
            if word == '':
                continue
            modifiers = []
            if re.search(r'\*', word) != None:
                modifiers.append('*')
            if re.search(r'&', word) != None:
                modifiers.append('&')
            if re.search(r'\bconst\b', word) != None:
                modifiers.append('const')
            word = re.sub(r'\*|&|(\bconst\b)', '', word).strip()
            w = word.split()
            type = w[0]
            name = w[1]
            parameters.append(Parameter(name, type, modifiers, namespace))
        self.parameters_ = parameters

        # compute template arguments
        self.template_parameters_ = []
        m = re.search('<(.*)>', self.name())
        if m != None:
            text = m.group(1)
            text = re.sub(r',', '', text)
            self.template_parameters_ = text.rsplit(r'\s')

    def __repr__(self):
        return self.text

    # Returns the namespace of the function (or the class containing this member function)
    #
    def namespace(self):
        return self.namespace_

    # returns the name of the function including a namespace qualification, if available
    #
    # 'name'
    def qualified_name(self):
        return re.sub('\(.*', '', self.text)

    # returns the namespace qualifier of the function
    #
    def qualifier(self):
        if self.qualified_name().find('::') == -1:
            return ''
        return re.sub('::.*$', '::', self.qualified_name())

    # returns the name of the function without namespace qualification
    #
    def name(self):
        return re.sub(r'^.*::', '', self.qualified_name())

    # returns true if it is a template function
    #
    def is_template(self):
        return self.name().find('<') != -1

    # returns the template arguments of the function
    #
    def template_parameters(self):
        return self.template_parameters_

    # returns the argument text of the function
    #
    # 'const core::identifier_string& name, const data::data_expression_list& arguments'
    def argument_text(self):
        return self.argument_text_

    # returns the parameters of the function as a sequence of VariableDeclarations
    def parameters(self):
        return self.parameters_

# Represents a class member function
class MemberFunction:
    def __init__(self, classname, modifiers, return_type, name, arg):
        self.classname = classname
        self.modifiers = modifiers
        self.return_type = return_type
        self.name = name
        self.arg  = arg

    def expand_text(self, text):
        text = re.sub('<CLASSNAME>', self.classname, text)
        text = re.sub('<RETURN_TYPE>', self.return_type, text)
        text = re.sub('<NAME>', self.name, text)
        text = re.sub('<ARG>',  self.arg , text)
        return text

    def inline_definition(self):
        if 'i' in self.modifiers:
           text = '''    <RETURN_TYPE> <NAME>() const
    {
      return atermpp::down_cast<atermpp::aterm_int>(static_cast<const atermpp::aterm&>(*this)).value();
    }'''
        else:
           text = '''    const <RETURN_TYPE>& <NAME>() const
    {
      return atermpp::down_cast<<RETURN_TYPE>>((*this)[<ARG>]);
    }'''
        return self.expand_text(text)

    def declaration(self):
        text = '''    const <RETURN_TYPE>& <NAME>() const;'''
        return self.expand_text(text)

    def definition(self, inline = False):
        text = '''    <INLINE> const<RETURN_TYPE>& <CLASSNAME>::<NAME>() const
    {
      return atermpp::down_cast<<RETURN_TYPE>>((*this)[<ARG>]);
    }'''
        if inline:
            text = re.sub('<INLINE>',  'inline\n    ', text)
        else:
            text = re.sub('<INLINE>',  '', text)
        return self.expand_text(text)

# Represents a class constructor
class Constructor:
    def __init__(self, classname, modifiers, arguments, superclass, namespace, aterm, parameters, template_parameters):
        self.classname           = classname
        self.modifiers           = modifiers
        self.arguments           = arguments
        self.superclass          = superclass
        self.namespace           = namespace
        self.aterm               = aterm
        self.aterm_namespace     = 'core'
        self.parameters          = parameters
        self.template_parameters = template_parameters

    def expand_text(self, text):
        text = re.sub('<CLASSNAME>'          , self.classname          , text)
        text = re.sub('<ARGUMENTS>'          , self.arguments          , text)
        text = re.sub('<SUPERCLASS>'         , self.superclass         , text)
        text = re.sub('<ATERM_NAMESPACE>'    , self.aterm_namespace    , text)
        text = re.sub('<ATERM>'              , self.aterm              , text)
        text = re.sub('<PARAMETERS>'         , self.parameters         , text)
        text = re.sub('<TEMPLATE_PARAMETERS>', self.template_parameters, text)
        if len(self.parameters)>0 and self.parameters.find(',')<0:     # There is only one parameter. Add explicit. 
            text = re.sub('<EXPLICIT>','explicit ',text)
        else:
           text = re.sub('<EXPLICIT>','',text)
        return text

    def inline_definition(self):
        if self.superclass == 'atermpp::aterm':
            text = r'''    /// \\\\brief Constructor Z12.
    <EXPLICIT><CLASSNAME>(<ARGUMENTS>)
      : atermpp::aterm(core::detail::function_symbol_<ATERM>(), <PARAMETERS>)
    {}'''
        elif 'i' in self.modifiers:
            text = r'''    /// \\\\brief Constructor Z13.
    <CLASSNAME>(<ARGUMENTS>)
      : <SUPERCLASS>(atermpp::aterm_int(<PARAMETERS>))
    {}'''
        else:
            text = r'''    /// \\\\brief Constructor Z14.
    <EXPLICIT><CLASSNAME>(<ARGUMENTS>)
      : <SUPERCLASS>(atermpp::aterm(core::detail::function_symbol_<ATERM>(), <PARAMETERS>))
    {}'''
        return self.expand_text(text)

    def declaration(self):
        text = r'''    /// \\\\brief Constructor Z15.
    <CLASSNAME>(<ARGUMENTS>);'''
        return self.expand_text(text)

    def definition(self, inline = False):
        if self.superclass == 'atermpp::aterm':
            text = r'''    /// \\\\brief Constructor Z16.
    <EXPLICIT><INLINE><CLASSNAME>::<CLASSNAME>(<ARGUMENTS>)
      : atermpp::aterm(core::detail::function_symbol_<ATERM>(), <PARAMETERS>)
    {}'''
        elif 'i' in self.modifiers:
            text = r'''    /// \\\\brief Constructor Z17.
    <INLINE><CLASSNAME>::<CLASSNAME>(<ARGUMENTS>)
      : <SUPERCLASS>(atermpp::aterm_int(<PARAMETERS>))
    {}'''
        else:
            text = r'''    /// \\\\brief Constructor Z18.
    <EXPLICIT><INLINE><CLASSNAME>::<CLASSNAME>(<ARGUMENTS>)
      : <SUPERCLASS>(atermpp::aterm(core::detail::function_symbol_<ATERM>(), <PARAMETERS>))
    {}'''
        if inline:
            text = re.sub('<INLINE>',  'inline\n    ', text)
        else:
            text = re.sub('<INLINE>',  '', text)
        return self.expand_text(text)

# Represents a default class constructor
class DefaultConstructor(Constructor):
    def __init__(self, classname, modifiers, arguments, superclass, namespace, aterm, parameters, template_parameters):
        self.classname           = classname
        self.modifiers           = modifiers
        self.arguments           = arguments
        self.superclass          = superclass
        self.namespace           = namespace
        self.aterm               = aterm
        self.aterm_namespace     = 'core'
        self.parameters          = parameters
        self.template_parameters = template_parameters

    def inline_definition(self):
        if 'i' in self.modifiers:
          text = r'''    /// \\\\brief Default constructor X2.
    <CLASSNAME>()
      : <SUPERCLASS>(atermpp::aterm_int(std::size_t(0)))
    {}'''
        else:
          text = r'''    /// \\\\brief Default constructor X3.
    <CLASSNAME>()
      : <SUPERCLASS>(<ATERM_NAMESPACE>::detail::default_values::<ATERM>)
    {}'''
        return self.expand_text(text)

    def declaration(self):
        text = r'''    /// \\\\brief Default constructor X4.
    <CLASSNAME>();'''
        return self.expand_text(text)

    def definition(self, inline = False):
        text = r'''    /// \\\\brief Default constructor X5.
    <INLINE><CLASSNAME>::<CLASSNAME>()
      : <SUPERCLASS>(<ATERM_NAMESPACE>::detail::default_values::<ATERM>)
    {}'''
        if inline:
            text = re.sub('<INLINE>',  'inline\n    ', text)
        else:
            text = re.sub('<INLINE>',  '', text)
        return self.expand_text(text)

# Represents an overloaded class constructor
class OverloadedConstructor(Constructor):
    def __init__(self, classname, modifiers, arguments, superclass, namespace, aterm, parameters, template_parameters):
        self.classname           = classname
        self.modifiers           = modifiers
        self.arguments           = arguments
        self.superclass          = superclass
        self.namespace           = namespace
        self.aterm               = aterm
        self.aterm_namespace     = 'core'
        self.parameters          = parameters
        self.template_parameters = template_parameters

    def inline_definition(self):
        if self.superclass == 'atermpp::aterm':
            text = r'''    /// \\\\brief Constructor Z1.
    <TEMPLATE_PARAMETERS><CLASSNAME>(<ARGUMENTS>)
      : atermpp::aterm(core::detail::function_symbol_<ATERM>(), <PARAMETERS>)
    {}'''
        else:
            text = r'''    /// \\\\brief Constructor Z2.
    <TEMPLATE_PARAMETERS><CLASSNAME>(<ARGUMENTS>)
      : <SUPERCLASS>(atermpp::aterm(core::detail::function_symbol_<ATERM>(), <PARAMETERS>))
    {}'''
        return self.expand_text(text)

    def declaration(self):
        text = r'''    /// \\\\brief Constructor Z3.
    <TEMPLATE_PARAMETERS><CLASSNAME>(<ARGUMENTS>);'''
        return self.expand_text(text)

    def definition(self, inline = False):
        if 'i' in self.modifiers:
             text = r'''    /// \\\\brief Default constructor X1.
      <INLINE><CLASSNAME>::<CLASSNAME>()
        : <SUPERCLASS>(atermpp::aterm_int(std::size_t(0)))
      {}'''
        elif self.superclass == 'atermpp::aterm':
            text = r'''    /// \\\\brief Constructor Z4.
    <TEMPLATE_PARAMETERS><INLINE><CLASSNAME>::<CLASSNAME>(<ARGUMENTS>)
      : atermpp::aterm(core::detail::function_symbol_<ATERM>(), <PARAMETERS>)
    {}'''
        else:
            text = r'''    /// \\\\brief Constructor Z5.
    <TEMPLATE_PARAMETERS><INLINE><CLASSNAME>::<CLASSNAME>(<ARGUMENTS>)
      : <SUPERCLASS>(atermpp::aterm(core::detail::function_symbol_<ATERM>(), <PARAMETERS>))
    {}'''
        if inline and len(self.template_parameters) == 0:
            text = re.sub('<INLINE>', 'inline', text)
        else:
            text = re.sub('<INLINE>', '', text)
        return self.expand_text(text)

# Represents a class constructor taking an additional class as argument
class AdditionalConstructor(Constructor):
    def __init__(self, classname, additional_classname, superclass, additional_class):
        self.classname            = classname
        self.additional_classname = additional_classname
        self.superclass           = superclass

    def inline_definition(self):
        text = r'''    /// \\\\brief Constructor Z6.
    <CLASSNAME>(const <ADDITIONAL_CLASSNAME>& x)
      : <SUPERCLASS>(x)
    {}'''
        return self.expand_text(text)

    def declaration(self):
        text = r'''    /// \\\\brief Constructor Z7.
    <CLASSNAME>(const <ADDITIONAL_CLASSNAME>& x);'''
        return self.expand_text(text)

    def expand_text(self, text):
        text = re.sub('<CLASSNAME>', self.classname, text)
        text = re.sub('<ADDITIONAL_CLASSNAME>', self.additional_classname, text)
        text = re.sub('<SUPERCLASS>', self.superclass, text)
        return text

    def definition(self, inline = False):
        text = r'''    /// \\\\brief Constructor Z8.
    /// \\\\param term A term
    <INLINE><CLASSNAME>::<CLASSNAME>(const <ADDITIONAL_CLASSNAME>& x)
      : <SUPERCLASS>(x)
    {}'''
        if inline:
            text = re.sub('<INLINE>',  'inline\n    ', text)
        else:
            text = re.sub('<INLINE>',  '', text)
        return self.expand_text(text)

# Represents a class constructor taking an aterm as argument
class ATermConstructor(Constructor):
    def __init__(self, classname, modifiers, arguments, superclass, namespace, aterm, parameters, template_parameters):
        self.classname           = classname
        self.modifiers           = modifiers
        self.arguments           = arguments
        self.superclass          = superclass
        self.namespace           = namespace
        self.aterm               = aterm
        self.aterm_namespace     = 'core'
        self.parameters          = parameters
        self.template_parameters = template_parameters

    def inline_definition(self):
        if 'i' in self.modifiers:
           text = r'''    /// \\\\brief Constructor based on an aterm.
    /// \\param term A term
    explicit <CLASSNAME>(const atermpp::aterm& term)
      : <SUPERCLASS>(term)
    {
      assert(this->type_is_int());
    }'''
        else:
          text = r'''    /// \\\\brief Constructor Z9.
    /// \\\\param term A term
    explicit <CLASSNAME>(const atermpp::aterm& term)
      : <SUPERCLASS>(term)
    {
      assert(<ATERM_NAMESPACE>::detail::check_term_<ATERM>(*this));
    }'''
        return self.expand_text(text)

    def declaration(self):
        text = r'''    /// \\\\brief Constructor Z10.
    /// \\\\param term A term
    explicit <CLASSNAME>(const atermpp::aterm& term);'''
        return self.expand_text(text)

    def definition(self, inline = False):
        if 'i' in self.modifiers:
          text = r'''    /// \\\\brief Constructor based on an aterm.
    /// \\param term A term
    explicit <INLINE><CLASSNAME>::<CLASSNAME>(const atermpp::aterm& term)
      : <SUPERCLASS>(term)
    {
      assert(this->type_is_int());
    }'''
        else:
          text = r'''    /// \\\\brief Constructor Z11.
    /// \\\\param term A term
    explicit <INLINE><CLASSNAME>::<CLASSNAME>(const atermpp::aterm& term)
      : <SUPERCLASS>(term)
    {
      assert(<ATERM_NAMESPACE>::detail::check_term_<ATERM>(*this));
    }'''
        if inline:
            text = re.sub('<INLINE>',  'inline\n    ', text)
        else:
            text = re.sub('<INLINE>',  '', text)
        return self.expand_text(text)

# Represents a class definition
#
# it can be initialized with a string like this:
#
# ActTrue   | true_()  | The value true for action formulas
#
# self.aterm:       the name of the corresponding aterm
# self.constructor: the constructor of the class
# self.description: a string description of the class
class Class:
    def __init__(self, aterm, constructor, description, superclass = None, namespace = None, modifiers = ''):
        self.aterm = aterm
        self.description = description
        name = re.sub('\(.*', '', constructor)
        arguments = re.sub('.*\(', '(', constructor)
        self.classname_ = re.sub('\[[^]]*\]', '', name)
        self.constructor = FunctionDeclaration(self.classname_ + arguments, namespace)
        self.superclass_ = remove_namespace(superclass)
        self.superclass_namespace_ = extract_namespace(superclass)
        self.namespace_ = namespace
        self.modifiers_ = modifiers
        self.expression_classes_ = []

    # Returns the derived expression classes of this class
    #
    def expression_classes(self):
        return self.expression_classes_

    # Returns the 'modifiers' of the class
    #
    # C : add container typedefs (list + vector)
    def modifiers(self):
        return self.modifiers_

    # Returns the namespace of the class
    #
    def namespace(self):
        return self.namespace_

    # Returns the namespace of the superclass
    #
    def superclass_namespace(self):
        return self.superclass_namespace_

    # Returns the name of the class
    #
    def classname(self, include_namespace = False):
        result = self.classname_
        if include_namespace:
            if result.find('::') == -1:
                result = '%s::%s' % (self.namespace(), result)
        return result

    # Returns the name of the class without namespace qualification
    #
    def name(self):
        return self.constructor.name()

    # Returns the superclass of the class, or None if no superclass exists
    #
    def superclass(self, include_namespace = False):
        if self.namespace() != self.superclass_namespace():
            include_namespace = True
        if include_namespace and self.superclass_namespace() != None:
            return '%s::%s' % (self.superclass_namespace(), self.superclass_)
        else:
            return self.superclass_

    # Returns the name of the class including a namespace qualification, if available
    #
    def qualified_name(self):
        return self.constructor.qualified_name()

    # Returns the namespace qualifier of the class (or '' if not available)
    def qualifier(self):
        return self.constructor.qualifier()

    # Returns true if the attributes of the class can be modified 'in place'
    def modifiable(self):
        return self.aterm == None or 'M' in self.modifiers()

    # Returns the constructors of the class
    def constructors(self):
        add_constructor_overloads = 'O' in self.modifiers()
        add_string_overload_constructor = add_constructor_overloads,
        add_container_overload_constructor = add_constructor_overloads
        classname = self.classname()
        modifiers = self.modifiers()
        arguments = self.constructor.argument_text()
        superclass = self.superclass()
        superclass_namespace = self.superclass_namespace()
        namespace = self.namespace()
        aterm = self.aterm

        # for constructor
        parameters = []
        arguments = []

        # for overloaded constructor
        parameters1 = []
        arguments1 = []
        arguments2 = []
        template_parameters = []

        for i, p in enumerate(self.constructor.parameters()):
            parameters.append(p.name())
            arguments.append('%s %s' % (p.type(not 'i' in self.modifiers()), p.name()))
            if p.type(False) == 'core::identifier_string' and add_string_overload_constructor:
                parameters1.append('core::identifier_string(%s)' % p.name())
                arguments1.append('const std::string& %s' % p.name())
            elif p.is_list() and add_container_overload_constructor:
                if len(template_parameters) > 0:
                    template_parameter = 'Container%d' % (i+1)
                else:
                    template_parameter = 'Container'
                template_parameters.append(template_parameter)
                arguments1.append('const %s& %s' % (template_parameter, p.name()))
                arguments2.append('typename atermpp::enable_if_container<%s, %s>::type* = nullptr' % (template_parameter, p.type(False)[:-5]))
                parameters1.append('%s(%s.begin(), %s.end())' % (p.type(False), p.name(),p.name()))
            else:
                parameters1.append(p.name())
                arguments1.append('%s %s' % (p.type(), p.name()))

        parameters  = ', '.join(parameters)
        arguments   = ', '.join(arguments)

        if len(template_parameters) > 0:
            template_parameters = 'template <typename %s>\n    ' % ', typename'.join(template_parameters)
        else:
            template_parameters = ''
        parameters1 = ', '.join(parameters1)
        arguments1  = ', '.join(arguments1 + arguments2)

        if superclass == None:
            superclass = 'atermpp::aterm'

        constructors = []
        constructors.append(DefaultConstructor(classname, modifiers, arguments, superclass, namespace, aterm, parameters, template_parameters))
        constructors.append(ATermConstructor(classname, modifiers, arguments, superclass, namespace, aterm, parameters, template_parameters))
        if len(self.constructor.parameters()) > 0:
            constructors.append(Constructor(classname, modifiers, arguments, superclass, namespace, aterm, parameters, template_parameters))
        if len(self.constructor.parameters()) > 0 and (add_string_overload_constructor or add_container_overload_constructor) and (parameters != parameters1):
            constructors.append(OverloadedConstructor(classname, modifiers, arguments1, superclass, namespace, aterm, parameters1, template_parameters))
        if self.classname(True) in ADDITIONAL_EXPRESSION_CLASS_DEPENDENCIES:
            for additional_classname in ADDITIONAL_EXPRESSION_CLASS_DEPENDENCIES[self.classname(True)]:
                constructors.append(AdditionalConstructor(classname, additional_classname, superclass, additional_classname))
        return constructors

    # Returns a specialization of the swap function for the std namespace
    # TODO: deal with template classes
    def swap_overload(self):
        text = r'''/// \\brief swap overload
inline void swap(<CLASSNAME>& t1, <CLASSNAME>& t2)
{
  t1.swap(t2);
}
'''
        return self.expand_text(text, parameters = '', constructors = '', member_functions = '', namespace = self.namespace())

    # Returns the member functions of the super class.
    # If the super class is atermpp::aterm or None, the empty list is returned
    def superclass_member_functions(self, all_classes):
        result = []
        if self.superclass() != None:
            name = self.superclass(include_namespace = True)
            if not name in ['atermpp::aterm']:
                c = all_classes[name]
                result = c.member_functions(all_classes, False)
                result = [f.name for f in result]
        return result

    # Returns the member functions of the class
    def member_functions(self, all_classes, exclude_superclass_member_functions = True):
        if exclude_superclass_member_functions:
            skip = self.superclass_member_functions(all_classes)
        else:
            skip = []
        result = []
        index = 1
        for p in self.constructor.parameters():
            arg, n = str(p), index
            index = index + 1
            p = arg.rpartition(' ')
            return_type = extract_type(p[0].strip())
            name = p[2].strip()
            arg = str(n - 1)
            if not name in skip:
                result.append(MemberFunction(self.classname(), self.modifiers(), return_type, name, arg))
        return result

    def move_semantics_text(self):
        text = r'''

    /// Move semantics
    <CLASSNAME>(const <CLASSNAME>&) noexcept = default;
    <CLASSNAME>(<CLASSNAME>&&) noexcept = default;
    <CLASSNAME>& operator=(const <CLASSNAME>&) noexcept = default;
    <CLASSNAME>& operator=(<CLASSNAME>&&) noexcept = default;'''
        return re.sub('<CLASSNAME>', self.classname(), text)

    def expand_text(self, text, parameters = '', constructors = '', member_functions = '', namespace = None, move_semantics = ''):
        superclass = self.superclass()
        if superclass == None:
            superclass_declaration = ': public atermpp::aterm'
        else:
            superclass_declaration = ': public %s' % superclass
        if 'U' in self.modifiers():
            user_section = '\n//--- start user section %s ---//\n//--- end user section %s ---//' % (self.classname(), self.classname())
        else:
            user_section = ''
        if namespace == None:
            namespace = self.namespace()
        text = re.sub('<SUPERCLASS_DECLARATION>', superclass_declaration, text)
        text = re.sub('<DESCRIPTION>'     , self.description, text)
        text = re.sub('<CLASSNAME>'       , self.classname(), text)
        text = re.sub('<ATERM>'           , self.aterm, text)
        text = re.sub('<CONSTRUCTORS>'    , constructors, text)
        text = re.sub('<PARAMETERS>'      , parameters, text)
        text = re.sub('<NAMESPACE>'       , namespace, text)
        text = re.sub('<USER_SECTION>'    , user_section, text)
        if superclass != None:
            text = re.sub('<SUPERCLASS>'  , superclass, text)
        text = re.sub('<MEMBER_FUNCTIONS>', member_functions, text)
        text = re.sub('<MOVE_SEMANTICS>'  , move_semantics, text)
        return text

    # Generate the make_...  functions that allow the construction of a class member in situ. 
    def make_function(self):
        if 'i' in self.modifiers():
          text = r'''/// \\brief Make_<CLASSNAME> constructs a new term into a given address.
/// \\ \param t The reference into which the new <CLASSNAME> is constructed. 
template <class... ARGUMENTS>
inline void make_<CLASSNAME>(atermpp::aterm& t, size_t n)
{
  atermpp::make_aterm_int(reinterpret_cast<atermpp::aterm_int&>(t), n);
}'''
        else:
          text = r'''/// \\brief Make_<CLASSNAME> constructs a new term into a given address.
/// \\ \param t The reference into which the new <CLASSNAME> is constructed. 
template <class... ARGUMENTS>
inline void make_<CLASSNAME>(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl<HASINDEX>(t, core::detail::function_symbol_<ATERM>(), args...);
}'''
        text = re.sub('<ATERM>', self.aterm, text)
        text = re.sub('<ARGUMENTS>', ', '.join([p.type() + ' ' + p.name() for p in self.constructor.parameters()]), text)
        text = re.sub('<PARAMETERS>', ', '.join([p.name() for p in self.constructor.parameters()]), text)
        if 'N' in self.modifiers():
          if len(self.constructor.parameters()) == 1:
              text = re.sub('<HASINDEX>', '_with_index<<CLASSNAME>,<PARAMETER_SORTS>>', text)
          else:
              text = re.sub('<HASINDEX>', '_with_index<<CLASSNAME>,std::pair<<PARAMETER_SORTS>>>', text)
        else:
          text = re.sub('<HASINDEX>', '', text)
        text = re.sub('<CLASSNAME>', self.classname(), text)
        text = re.sub('<PARAMETER_SORTS>',', '.join([re.sub('const ','',re.sub('&','',p.type())) for p in self.constructor.parameters()]), text)
        return text

    # Returns typedefs for term lists and term vectors.
    def container_typedefs(self):
        text = r'''/// \\brief list of <CLASSNAME>s
typedef atermpp::term_list<<CLASSNAME>> <CLASSNAME>_list;

/// \\brief vector of <CLASSNAME>s
typedef std::vector<<CLASSNAME>>    <CLASSNAME>_vector;'''
        text = re.sub('<CLASSNAME>', self.classname(), text)
        return text

    # Returns an is_<classname> function
    def is_function(self, all_classes):
        text = r'''%s/// \\brief Test for a %s expression
/// \\param x A term
/// \\return True if \\a x is a %s expression
inline
bool is_%s(const %s& x)
{
  return %s;
}'''
        name = self.name()
        if name[-1] == '_':
            name = name[:-1]
        if 'X' in self.modifiers():
            prototypes = '// prototypes\n' + '\n'.join(['inline bool is_%s(const atermpp::aterm& x);' % re.sub('_$', '', c.classname()) for c in self.derived_classes(all_classes) if c.namespace() == self.namespace()]) + '\n\n'
            return_value = ' ||\n         '.join(['%s::is_%s(x)' % (c.namespace(), re.sub('_$', '', c.classname())) for c in self.derived_classes(all_classes)])
        else:
            prototypes = ''
            return_value = 'x.function() == core::detail::function_symbols::%s' % self.aterm
        text = text % (prototypes, name, name, name, 'atermpp::aterm', return_value)
        return text

    def is_function_name(self, include_namespace = True):
        name = self.name()
        if name[-1] == '_':
            name = name[:-1]
        result = 'is_' + name
        if include_namespace:
            result = self.namespace() + '::' + result
        return result

    # Returns the class definition
    def class_inline_definition(self, all_classes):
        text = ''

        # generate the class definition
        if not 'S' in self.modifiers():
            ptext = ', '.join([p.name() for p in self.constructor.parameters()])
            if 's' in self.modifiers():
                ctext = ''
            else:
                ctext = '\n\n'.join([x.inline_definition() for x in self.constructors()])
            mtext = ''.join(['\n\n' + x.inline_definition() for x in self.member_functions(all_classes)])

            text = r'''/// \\brief <DESCRIPTION>
class <CLASSNAME><SUPERCLASS_DECLARATION>
{
  public:
<CONSTRUCTORS><MOVE_SEMANTICS><MEMBER_FUNCTIONS><USER_SECTION>
};'''
            text = self.expand_text(text = text,
                                    parameters = ptext,
                                    constructors = ctext,
                                    member_functions = mtext,
                                    namespace = self.namespace(),
                                    move_semantics = self.move_semantics_text()
                                   )

            if 'X' in self.modifiers():
                text = re.sub('check_term', 'check_rule', text)

        # generate additional functions
        # if not ('s' in self.modifiers() or 'S' in self.modifiers()):
        if not 'S' in self.modifiers() and self.constructor.parameters():
            text = text + '\n\n' + self.make_function()
        if 'C' in self.modifiers():
            text = text + '\n\n' + self.container_typedefs()
        if 'I' in self.modifiers():
            text = text + '\n\n' + self.is_function(all_classes)
        text = text + '\n\n' + self.ostream_overload()
        if not 'W' in self.modifiers():
            text = text + '\n\n' + self.swap_overload()

        return text + '\n'

    # Returns an overload of operator<<(ostream&)
    def ostream_overload(self):
        text = r'''// prototype declaration
std::string pp(const <CLASSNAME>& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const <CLASSNAME>& x)
{
  return out << <NAMESPACE>::pp(x);
}'''
        text = self.expand_text(text)
        return text

    # Returns the class declaration
    def class_declaration(self, all_classes, namespace = 'core'):
        ptext = ', '.join([p.name() for p in self.constructor.parameters()])
        ctext = '\n\n'.join([x.declaration() for x in self.constructors(namespace)])
        mtext = ''.join(['\n\n' + x.declaration() for x in self.member_functions(all_classes)])

        text = r'''/// \\brief <DESCRIPTION>
class <CLASSNAME><SUPERCLASS_DECLARATION>
{
  public:
<CONSTRUCTORS><MEMBER_FUNCTIONS>
};'''
        text = self.expand_text(text, ptext, ctext, mtext, namespace)

        if 'C' in self.modifiers():
            text = text + '\n\n' + self.container_typedefs()
        if 'I' in self.modifiers():
            text = text + '\n\n' + self.is_function(all_classes)
        text = text + '\n\n' + self.ostream_overload()
        if not 'W' in self.modifiers():
            text = text + '\n\n' + self.swap_overload()

        return text + '\n'

    def derived_classes(self, all_classes):
        classes = [all_classes[name] for name in self.expression_classes()]
        classes.sort(key = lambda x: x.index)
        classes.sort(key = lambda x: 'X' in x.modifiers(), reverse=True)
        return classes

    def traverser_function(self, all_classes, dependencies):
        text = r'''void apply(const <CLASS_NAME>& x)
{
  static_cast<Derived&>(*this).enter(x);<VISIT_TEXT>
  static_cast<Derived&>(*this).leave(x);
}
'''
        classname = self.classname(True)
        visit_text = ''
        dependent = False
        if 'X' in self.modifiers():
            classes = self.derived_classes(all_classes)
            updates = []
            for c in classes:
                is_function = c.is_function_name(True)
                if 'M' in c.modifiers():
                    cast = '%s(atermpp::down_cast<atermpp::aterm>(x))' % c.classname(True)
                else:
                    cast = 'atermpp::down_cast<%s>(x)' % c.classname(True)
                updates.append('''if (%s(x))
{
  static_cast<Derived&>(*this).apply(%s);
}''' % (is_function, cast))
            if len(updates) == 0:
                visit_text = '// skip'
            else:
                visit_text = '\nelse '.join(updates)
        else:
            if 'E' in self.modifiers():
                return_type = self.superclass(include_namespace = True)
            else:
                return_type = classname

            updates = []
            f = self.constructor
            for p in f.parameters():
                ptype = p.type(include_modifiers = False, include_namespace = True)
                if is_dependent_type(dependencies, ptype):
                    dependent = True
                    # special case for arguments of a data application
                    if self.classname(True) == 'data::application' and p.name() == 'arguments':
                        update = 'for (const data_expression& t: x) { static_cast<Derived&>(*this).apply(t); }'
                    else:
                        update = 'static_cast<Derived&>(*this).apply(x.%s());' % p.name()
                    # special case for stochastic distribution
                    if self.classname(True) == 'lps::stochastic_distribution' and p.name() == 'distribution':
                        update = 'if (x.is_defined()) { ' + update + ' }'
                    updates.append(update)
            if dependent:
                visit_text = '\n'.join(updates)
            else:
                visit_text = '// skip'

        # fix the layout
        if visit_text != '':
            visit_text = '\n' + indent_text(visit_text, '  ')

        text = re.sub('<CLASS_NAME>', classname, text)
        text = re.sub('<VISIT_TEXT>', visit_text, text)
        if self.constructor.is_template():
            text = 'template <typename ' + ', typename '.join(f.template_parameters()) + '>\n' + text
        return text

    #def builder_return_type(self, all_classes, modifiability_map):
    #    classname = self.classname(True)
    #    # N.B. the order of the statements below is important!
    #    if is_modifiable_type(classname, modifiability_map):
    #        result = 'void'
    #    elif 'E' in self.modifiers():
    #        result = self.superclass(include_namespace = True)
    #        c = all_classes[result]
    #        if 'E' in c.modifiers():
    #            result = c.builder_return_type(all_classes, modifiability_map)
    #    else:
    #        result = classname
    #    return result

    def builder_function(self, all_classes, dependencies, modifiability_map):
        text = r'''<TEMPLATE>void <METHOD>(<RESULT><CONST><CLASS_NAME>& x)
{ 
  <ASSERT>static_cast<Derived&>(*this).enter(x);<VISIT_TEXT>
  static_cast<Derived&>(*this).leave(x);<RETURN_STATEMENT>
}
'''
        classname = self.classname(True)
        visit_text = ''
        dependent = False

        method = 'apply'
        if is_modifiable_type(classname, modifiability_map):
            template = ''
            return_type = 'void'
            const = ''
            method = 'update'
            result = ''
            assertion = ''
        else:
            const = 'const '
            template = 'template <class T>'
            # We currently return the same class as we get passed as parameter.
            # The following statement was not present in previous versions of this
            # file.
            return_type = classname

            # It turns out that the builder_return_type function that computes the return
            # type has been bypassed. Unfortunately the new solution does not handle the
            # artificial layer of data abstractions correctly. This is a quick fix to deal
            # with it.
            if 'data::abstraction' in return_type:
                return_type = 'data::data_expression'
            result = 'T& result, '
            assertion = 'assert(&result!=&x);'

        if is_modifiable_type(classname, modifiability_map):
            return_statement = ''

            updates = []
            f = self.constructor
            for p in f.parameters():
                ptype = p.type(include_modifiers = False, include_namespace = True, remove_templates = True)
                pclass = None
                if ptype in all_classes:
                    pclass = all_classes[ptype]
                qtype = p.type(include_modifiers = False, include_namespace = True)
                if is_dependent_type(dependencies, ptype):
                    dependent = True
                    if is_modifiable_type(qtype, modifiability_map):
                        updates.append('static_cast<Derived&>(*this).update(x.%s());' % p.name())
                    else:
                        local_type = re.sub('const ','',re.sub('&','',p.type()))
                        if classname == 'lps::stochastic_specification' and local_type == 'process_initializer':
                            local_type = 'stochastic_process_initializer'   # Unclear why this needs to be done. Appears to be a bug. 
                        local_variable = 'result_%s' % p.name()
                        updates.append('%s %s;\nstatic_cast<Derived&>(*this).apply(%s, x.%s());\nx.%s() = %s;' \
                                         % (local_type, local_variable, local_variable, p.name(), p.name(), local_variable))
                else:
                    continue
            if dependent:
                visit_text = '\n'.join(updates)
            else:
                visit_text = '// skip'
        else:
            if 'X' in self.modifiers():
                classes = self.derived_classes(all_classes)
                updates = []
                for c in classes:
                    is_function = c.is_function_name(True)
                    if 'M' in c.modifiers():
                        cast = '%s(atermpp::down_cast<atermpp::aterm>(x))' % c.classname(True)
                    else:
                        cast = 'atermpp::down_cast<%s>(x)' % c.classname(True)
                    updates.append('''if (%s(x))
{
  static_cast<Derived&>(*this).apply(result, %s);
}''' % (is_function, cast))
                if len(updates) == 0:
                    visit_text = '// skip'
                    return_statement = 'result = x;'
                else:
                    visit_text = '\nelse '.join(updates)
                    return_statement = ''
            else:
                updates = []
                f = self.constructor
                for p in f.parameters():
                    ptype = p.type(include_modifiers = False, include_namespace = True)
                    pclass = None
                    if ptype in all_classes:
                        pclass = all_classes[ptype]
                    if is_dependent_type(dependencies, ptype):
                        dependent = True
                        updates.append('[&](%s result){ static_cast<Derived&>(*this).apply(result, x.%s()); }' \
                             % (re.sub('const ','',p.type()), p.name()))
                    else:
                        updates.append('x.%s()' % p.name())
                if dependent:
                    # special case for arguments of a data application
                    make_class_function = re.sub('::', '::make_', classname)
                    if self.classname(True) == 'data::application' and p.name() == 'arguments':
                        visit_text = '''data::make_application(result,
   x.head(),
   x.begin(),
   x.end(),
   [&](data_expression& result, const data::data_expression& t){ static_cast<Derived&>(*this).apply(result,t);} );''' 
                    # special case for stochastic distribution
                    elif return_type == 'lps::stochastic_distribution':
                        visit_text = 'result = x; if (x.is_defined()) { %s(result, %s); }' % (make_class_function, ', '.join(updates))
                    else:
                        visit_text = '%s(result, %s);' % (make_class_function, ', '.join(updates))
                    return_statement = ''
                else:
                    visit_text = '// skip'
                    return_statement = 'result = x;'

        # fix the layout
        if template != '':
            template = template + '\n'
        if return_statement != '':
            return_statement = '\n  ' + return_statement
        if assertion != '':
            assertion = return_statement +'\n  '
        if visit_text != '':
            visit_text = '\n' + indent_text(visit_text, '  ')

        text = re.sub('<TEMPLATE>', template, text)
        text = re.sub('<RETURN_TYPE>', return_type, text)
        text = re.sub('<RESULT>', result, text)
        text = re.sub('<CONST>', const, text)
        text = re.sub('<VISIT_TEXT>', visit_text, text)
        text = re.sub('<CLASS_NAME>', classname, text)
        text = re.sub('<RETURN_STATEMENT>', return_statement, text)
        text = re.sub('<METHOD>', method, text)
        text = re.sub('<ASSERT>', assertion, text)
        if self.constructor.is_template():
            text = 'template <typename ' + ', typename '.join(f.template_parameters()) + '>\n' + text
        return text

##########################################################################################################################

def extract_namespace(text):
    if text == None:
        return None
    text = re.sub('\(.*', '', text)
    pos = text.find(':')
    if pos >= 0:
        return text[:pos]
    else:
        return None

def remove_namespace(text):
    if text == None:
        return None
    pos = text.find('::')
    if pos >= 0:
        return text[pos+2:]
    else:
        return text

# parses lines that contain entries separated by '|'
# empty lines are removed
#
# example input:
#
# ActTrue   | true_()  | The value true for action formulas
# ActFalse  | false_() | The value false for action formulas
#
# Each line is split w.r.t. the '|' character; the words of the line
# are put in a tuple, and the sequence of tuples is returned
def parse_classes(text, namespace = None):
    result = []
    lines = text.rsplit('\n')
    superclass = None
    classes = {}
    for line in lines:
        if line.startswith('%'):
            continue
        words = [word.strip() for word in line.split('|')]
        if len(words) < 4:
            continue
        constructor, modifiers, aterm, description = words

        m = re.search(r'\:\s*public\s*(.*)$', constructor)
        if m != None:
            superclass = m.group(1)
            constructor = re.sub(r'\:\s*public\s*.*', '', constructor).strip()

        result.append(Class(aterm, constructor, description, superclass, namespace, modifiers))
    return result

def print_dependencies(dependencies, message):
    print(message)
    for type in sorted(dependencies):
        print((type, dependencies[type]))

def is_dependent_type(dependencies, type):
    if type in dependencies:
        return dependencies[type]
    m = re.search('<(.+)>', type)
    if m != None:
        return dependencies[m.group(1)]
    return False

def update_dependency(classname, all_classes, dependencies, value = True):
    dependencies[classname] = value
    if 'C' in all_classes[classname].modifiers():
        dependencies[classname + '_list'] = value
        dependencies[classname + '_vector'] = value

def update_dependencies(all_classes, dependencies):
    changed = False
    for classname in all_classes:
        c = all_classes[classname]
        if dependencies[classname] == True:

            # this is needed because of poor choices in the internal format
            #if 'X' in c.modifiers():
            #   for expr in all_classes[classname].expression_classes():
            #       update_dependency(expr, all_classes, dependencies, value = True)
            continue

        # check expression class dependencies
        if 'X' in c.modifiers():
            for expr in c.expression_classes():
                if dependencies[expr] == True:
                    update_dependency(classname, all_classes, dependencies)
                    changed = True
                    break

        # check parameter dependencies
        for p in c.constructor.parameters():
            type = p.type(include_modifiers = False, include_namespace = True, remove_templates = True)
            if is_dependent_type(dependencies, type):
                update_dependency(classname, all_classes, dependencies)
                changed = True

    return changed

def find_dependencies(all_classes, type):
    dependencies = {} # maps class names to True/False

    # initially set all dependencies to False
    for classname in all_classes:
        update_dependency(classname, all_classes, dependencies, value = False)

    # initial dependency: the expression class depends on itself, and all derivatives too
    update_dependency(type, all_classes, dependencies, value = True)

    #for expr in all_classes[type].expression_classes():
    #    update_dependency(expr, all_classes, dependencies, value = True)
    #    print(expr, dependencies[expr])

    while update_dependencies(all_classes, dependencies):
        pass
    return dependencies

# Computes a mapping m, such that m[<type>] returns true if <type> is a type that can be modified in place.
def make_modifiability_map(all_classes):
    result = {}
    for classname in all_classes:
        c = all_classes[classname]
        value = (c.aterm == None) or ('M' in c.modifiers())
        result[classname] = value
        if 'C' in c.modifiers():
            result[classname + '_list'] = False
            result[classname + '_vector'] = True
    return result

def is_modifiable_type(type, modifiability_map):
    if type in modifiability_map:
        return modifiability_map[type]
    elif type.startswith('std::vector<'):
        return True
    elif type.startswith('std::set<'):
        return True
    else:
        raise Exception('is_modifiable_type(' + type + ') is unknown')

def parse_class_map(class_map):
    result = {}
    index = 0 # give each class an index, used for sorting
    for namespace in class_map:
        class_text = class_map[namespace]
        classes = parse_classes(class_text, namespace)
        for c in classes:
            c.index = index
            index = index + 1
            classname = c.classname(include_namespace = True)
            result[classname] = c

    # set expression class dependencies
    for classname in list(result.keys()):
        c = result[classname]
        if 'E' in c.modifiers():
            superclass = c.superclass(include_namespace = True)
            result[superclass].expression_classes().append(classname)
        if 'X' in c.modifiers() and classname in ADDITIONAL_EXPRESSION_CLASS_DEPENDENCIES:
            for name in ADDITIONAL_EXPRESSION_CLASS_DEPENDENCIES[classname]:
                result[classname].expression_classes().append(name)
        if 'X' in c.modifiers():
            c.expression_classes().sort(key = lambda x: result[x].index)
    return result

def parse_classnames(text, namespace):
    result = []
    lines = text.rsplit('\n')
    classes = {}
    for line in lines:
        if line.startswith('%'):
            continue
        words = [word.strip() for word in line.split('|')]
        if len(words) < 4:
            continue
        constructor, modifiers, aterm, description = words
        classname = re.sub(r'\(.*', '', constructor)
        classname = '%s::%s' % (namespace, classname)

        # duplicates are not allowed
        if not classname in result:
            result.append(classname)
    return result

def mcrl2_class_map():
    return {
          'core'             : CORE_CLASSES,
          'data'             : DATA_EXPRESSION_CLASSES + ASSIGNMENT_EXPRESSION_CLASSES + SORT_EXPRESSION_CLASSES + CONTAINER_TYPES + BINDER_TYPES + ABSTRACTION_EXPRESSION_CLASSES + STRUCTURED_SORT_ELEMENTS + DATA_CLASSES,
          'state_formulas'   : STATE_FORMULA_CLASSES + MODAL_FORMULA_CLASSES,
          'regular_formulas' : REGULAR_FORMULA_CLASSES,
          'action_formulas'  : ACTION_FORMULA_CLASSES,
          'lps'              : LPS_CLASSES,
          'process'          : PROCESS_CLASSES + PROCESS_EXPRESSION_CLASSES,
          'pbes_system'      : PBES_CLASSES + PBES_EXPRESSION_CLASSES,
          'pres_system'      : PRES_CLASSES + PRES_EXPRESSION_CLASSES,
        }
