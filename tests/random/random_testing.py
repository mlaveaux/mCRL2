#!/usr/bin/env python

# Copyright 2015 Wieger Wesselink.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

import argparse
import os
import os.path
import random
import re
import shutil
import sys
import traceback

# Makes sure that the script can find the modules when ran directly.
sys.path.append(os.path.join(os.path.dirname(__file__),'../../'))

from tests.utility.process_expression import Literal
import tests.utility.random_process_expression as random_process_expression
import tests.utility.random_state_formula_generator as random_state_formula_generator
from tests.utility.random_bes_generator import make_bes
from tests.utility.random_pbes_generator import make_pbes
from tests.utility.testing import YmlTest
from tests.utility.text_utility import write_text

MCRL2_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

def ymlfile(file):
    return f'{MCRL2_ROOT}/tests/specifications/{file}.yml'

def mcrl2file(file):
    return os.path.join(MCRL2_ROOT, file)

class RandomTest(YmlTest):
    def __init__(self, name, _ymlfile, settings):
        super().__init__(name, _ymlfile, [], settings)

    # create input files for the random test, and add the filenames to self.inputfiles
    def create_inputfiles(self, runpath = '.'):
        raise NotImplementedError

    # removes input files that are in the runpath directory
    def remove_inputfiles(self, runpath = '.'):
        for filename in self.inputfiles:
            if os.path.abspath(runpath) == os.path.abspath(os.path.dirname(filename)):
                os.remove(filename)

    def execute(self, runpath = '.'):
        self.create_inputfiles(runpath)
        super().execute()
        self.remove_inputfiles(runpath)

class ProcessTest(RandomTest):
    def __init__(self, name, _ymlfile, settings):
        super().__init__(name, _ymlfile, settings)
        self.actions = ['a', 'b', 'c', 'd']
        self.process_identifiers = ['P', 'Q', 'R']
        self.process_size = 13
        self.parallel_operator_generators = random_process_expression.default_parallel_operator_generators
        self.process_expression_generators = random_process_expression.default_process_expression_generators
        self.init = None
        self.generate_process_parameters = False

    def create_inputfiles(self, runpath = '.'):
        filename = f'{self.name}.mcrl2'
        p = random_process_expression.make_process_specification(self.parallel_operator_generators, self.process_expression_generators, self.actions, self.process_identifiers, self.process_size, self.init, self.generate_process_parameters)
        write_text(filename, str(p))
        self.inputfiles += [filename]

# generates stochastic random processes
class StochasticProcessTest(ProcessTest):
    def __init__(self, name, _ymlfile, settings):
        super().__init__(name, _ymlfile, settings)
        self.process_expression_generators = {
            random_process_expression.ActionGenerator(): 8,
            random_process_expression.DeltaGenerator(): 1,
            random_process_expression.TauGenerator(): 1,
            random_process_expression.ProcessInstanceGenerator(): 2,
            random_process_expression.SumGenerator(): 2,
            random_process_expression.IfThenGenerator(): 2,
            random_process_expression.IfThenElseGenerator(): 2,
            random_process_expression.ChoiceGenerator(): 5,
            random_process_expression.SeqGenerator(): 5,
            random_process_expression.MultiActionGenerator(): 1,
            random_process_expression.DistGenerator(): 3,
        }

# generates random process with higher probability of tau transitions
class ProcessTauTest(ProcessTest):
    def __init__(self, name, testfile, settings):
        super().__init__(name, testfile, settings)
        self.actions = ['a', 'b', 'c']
        self.init = Literal(r'hide({a}, allow({a, b, c}, P || Q || R))')
        self.process_expression_generators = {
                        random_process_expression.ActionGenerator(): 8,
                        random_process_expression.DeltaGenerator(): 1,
                        random_process_expression.TauGenerator(): 4,
                        random_process_expression.ProcessInstanceGenerator(): 1,
                        random_process_expression.SumGenerator(): 0,
                        random_process_expression.IfThenGenerator(): 0,
                        random_process_expression.IfThenElseGenerator(): 0,
                        random_process_expression.ChoiceGenerator(): 5,
                        random_process_expression.SeqGenerator(): 5,
                        random_process_expression.MultiActionGenerator(): 1,
                        random_process_expression.DistGenerator(): 0,
                             }

class AlphabetReduceTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('alphabet-reduce'), settings)
        self.actions = ['a', 'b', 'c', 'd', 'e']

class LpsSuminstTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lpssuminst'), settings)

class LpsSumelmTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lpssumelm'), settings)

class LpsParelmTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lpsparelm'), settings)
        self.generate_process_parameters = True

class LpsOnePointRuleRewriteTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lpstransform'), settings)
        self.add_command_line_options('t2', ['-alps-one-point-rule-rewriter'])

class LpsConfcheckTest(ProcessTauTest):
    def __init__(self, name, confluence_type, settings):
        self.option_map = { 'commutative' : 'C',
                            'commutative-disjoint' : 'c',
                            'disjoint' : 'd',
                            'triangular' : 'T',
                            'trivial' : 'Z'
                          }
        assert confluence_type in self.option_map
        super().__init__(name, ymlfile('lpsconfcheck'), settings)
        self.add_command_line_options('t2', ['-x' + self.option_map[confluence_type]])

class LtscompareTest(ProcessTauTest):
    def __init__(self, name, equivalence_type, settings):
        assert equivalence_type in ['bisim', 'bisim-gv', 'bisim-gjkw', 'bisim-gj', 'bisim-gj-lazy-BLC', 'branching-bisim', 'branching-bisim-gv', 'branching-bisim-gjkw', 'branching-bisim-gj', 'branching-bisim-gj-lazy-BLC', 'dpbranching-bisim', 'dpbranching-bisim-gv', 'dpbranching-bisim-gjkw', 'dpbranching-bisim-gj', 'dpbranching-bisim-gj-lazy-BLC', 'weak-bisim', 'dpweak-bisim', 'sim', 'ready-sim' , 'trace', 'weak-trace']
        super().__init__(name, ymlfile('ltscompare'), settings)
        self.add_command_line_options('t3', ['-e' + equivalence_type])
        self.add_command_line_options('t4', ['-e' + equivalence_type])

class LtsCombineTest(ProcessTest):
    def __init__(self, name, parallel, settings):
        super().__init__(name, ymlfile('ltscombine'), settings)

        if parallel:
            self.add_command_line_options('t5', ['--threads=8'])


    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)

        # Create a second mCRL2 specification to combine
        filename = f'{self.name}2.mcrl2'
        p = random_process_expression.make_process_specification(self.parallel_operator_generators, self.process_expression_generators, self.actions, self.process_identifiers, self.process_size, self.init, self.generate_process_parameters)
        write_text(filename, str(p))
        self.inputfiles += [filename]

class LtscompareCounterexampleTest(ProcessTauTest):
    def __init__(self, name, equivalence_type, hide_actions, settings):
        assert equivalence_type in ['bisim', 'branching-bisim', 'trace', 'impossible-futures']
        super().__init__(name, ymlfile('ltscompare-counter-example'), settings)

        if equivalence_type == 'impossible-futures':
            self.add_command_line_options('t4', ['-p' + equivalence_type])
            self.add_command_line_options('t5', ['-p' + equivalence_type])
        else:
            self.add_command_line_options('t4', ['-e' + equivalence_type])
            self.add_command_line_options('t5', ['-e' + equivalence_type])

        if hide_actions:
            self.add_command_line_options('t4', ['--tau=b'])
            self.add_command_line_options('t5', ['--tau=b'])
            self.add_command_line_options('t6', ['--tau=b'])
            self.add_command_line_options('t7', ['--tau=b'])

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        # This is a hack to ensure that the counter example formula always exists for the further steps.
        filename = runpath + '/l7.mcf'
        write_text(filename, "true")

class StochasticLtscompareTest(StochasticProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('stochastic-ltscompare'), settings)

class BisimulationTest(ProcessTauTest):
    def __init__(self, name, equivalence_type, settings):
        assert equivalence_type in ['bisim', 'bisim-gv', 'bisim-gjkw', 'bisim-gj', 'bisim-gj-lazy-BLC', 'branching-bisim', 'branching-bisim-gv', 'branching-bisim-gjkw', 'branching-bisim-gj', 'branching-bisim-gj-lazy-BLC', 'weak-bisim']
        super().__init__(name, ymlfile('bisimulation'), settings)
        self.add_command_line_options('t3', ['-e' + equivalence_type])
        self.add_command_line_options('t4', ['-e' + equivalence_type])
        if equivalence_type in ['branching-bisim-gv', 'branching-bisim-gjkw', 'branching-bisim-gj', 'branching-bisim-gj-lazy-BLC']:
            self.add_command_line_options('t7', ['-bbranching-bisim'])
        elif equivalence_type in ['bisim', 'bisim-gv', 'bisim-gjkw', 'bisim-gj', 'bisim-gj-lazy-BLC']:
            self.add_command_line_options('t7', ['-bstrong-bisim'])
        else:
            self.add_command_line_options('t7', ['-b' + equivalence_type])

class Lps2ltsParallelTest(ProcessTauTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lps2lts_parallel'), settings)

class Lps2ltsAlgorithmsTest(ProcessTauTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lps2lts-algorithms'), settings)
        # randomly choose an algorithm
        actions = random.choice(['a', 'a,b', 'a,b,c'])
        options = [random.choice(['--deadlock', '--divergence', '--nondeterminism', f'--action={actions}'])]
        options = [random.choice(['--deadlock', '--nondeterminism', f'--action={actions}'])]
        if 'divergence' in options[0]:
            tau_actions = random.choice(['', '', 'b', 'b,c'])
            if tau_actions:
                options.append('--tau={}'.format(tau_actions))
        self.add_command_line_options('t2', options)
        self.add_command_line_options('t3', options)


class LpsConstelmTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lpsconstelm'), settings)
        self.generate_process_parameters = True

class LpsBinaryTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lpsbinary'), settings)
        self.generate_process_parameters = True

class LpsstategraphTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lpsstategraph'), settings)
        self.generate_process_parameters = True

class Lps2pbesTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lps2pbes'), settings)

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        self.inputfiles.append(mcrl2file('examples/modal-formulas/nodeadlock.mcf'))

class Lps2presTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lps2pres'), settings)

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        self.inputfiles.append(mcrl2file('examples/modal-formulas/nodeadlock.mcf'))

class LpsrewrTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lpsrewr'), settings)
        self.add_command_line_options('t2', ['-p' + 'prune-dataspec'])

class Lts2pbesTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('lts2pbes'), settings)

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        # Use the nodeadlock property to generate the PBES.
        self.inputfiles.append(mcrl2file('examples/modal-formulas/nodeadlock.mcf'))

class LtsconvertsymbolicTest(ProcessTest):
    def __init__(self, name, arguments, settings):
        super().__init__(name, ymlfile('ltsconvertsymbolic'), settings)

        if arguments:
            self.add_command_line_options('t3', arguments)


class PbesTest(RandomTest):
    def __init__(self, name, ymlfile, settings):
        super().__init__(name, ymlfile, settings)
        self.equation_count = 4
        self.atom_count = 4
        self.propvar_count = 3
        self.use_quantifiers = True
        self.use_integers = True

    def create_inputfiles(self, runpath = '.'):
        filename = f'{self.name}.txt'
        p = make_pbes(self.equation_count, self.atom_count, self.propvar_count, self.use_quantifiers, use_integers=self.use_integers)
        write_text(filename, str(p))
        self.inputfiles += [filename]

# N.B. does not work yet due to unusable abstraction map
class PbesabsintheTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbesabsinthe'), settings)

# N.B. This test has been disabled, since the tool has been deprecated.
class PbesabstractTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbesabstract'), settings)

class PbesconstelmTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbesconstelm'), settings)

class PbesfixpointsolveTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbesfixpointsolve'), settings)

class PbesparelmTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbesparelm'), settings)

class PbespareqelmTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbespareqelm'), settings)

class PbeschainTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbeschain'), settings)

class Pbespor1Test(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbespor1'), settings)

class Pbespor2Test(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbespor2'), settings)

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        filename = f'{self.name}.mcf'
        formula = random_state_formula_generator.make_modal_formula()
        write_text(filename, str(formula))
        self.inputfiles += [filename]

class PbesrewrTest(PbesTest):
    def __init__(self, name, rewriter, settings):
        super().__init__(name, ymlfile('pbesrewr'), settings)
        self.add_command_line_options('t2', ['-p' + rewriter])

class PbestransformTest(PbesTest):
    def __init__(self, name, rewriter, settings):
        super().__init__(name, ymlfile('pbestransform'), settings)
        self.add_command_line_options('t2', ['-a' + rewriter])

class PbesinstTest(PbesTest):
    def __init__(self, name, options, settings):
        super().__init__(name, ymlfile('pbesinst'), settings)
        self.add_command_line_options('t2', options)

class PbespgsolveTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbespgsolve'), settings)

class PbesstategraphTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbesstategraph'), settings)

class PbessymbolicbisimTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbessymbolicbisim'), settings)

class PbessolvesymbolicTest(PbesTest):
    def __init__(self, name, arguments, settings):
        super().__init__(name, ymlfile('pbessolvesymbolic'), settings)

        if arguments:
            self.add_command_line_options('t3', arguments)

class PbessolvesymbolicCounterexampleTest(ProcessTest):
    def __init__(self, name, arguments, settings):
        super().__init__(name, ymlfile('pbessolvesymbolic-counter-example'), settings)

        if arguments:
            self.add_command_line_options('t3', arguments)

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        filename = f'{self.name}.mcf'
        formula = random_state_formula_generator.make_modal_formula()
        write_text(filename, str(formula))
        self.inputfiles += [filename]

class Pbes2boolTest(PbesTest):
    def __init__(self, name, parallel, settings):
        super().__init__(name, ymlfile('pbessolve'), settings)

        if parallel:
            self.add_command_line_options('t2', ["--threads=4"])
            self.add_command_line_options('t3', ["--threads=4"])
            self.add_command_line_options('t4', ["--threads=4"])
            self.add_command_line_options('t5', ["--threads=4"])
            self.add_command_line_options('t6', ["--threads=4"])
            self.add_command_line_options('t7', ["--threads=4"])
            self.add_command_line_options('t8', ["--threads=4"])
            self.add_command_line_options('t9', ["--threads=4"])

class Pres2boolTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pressolve'), settings)

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        filename = f'{self.name}.mcf'
        formula = random_state_formula_generator.make_modal_formula()
        write_text(filename, str(formula))
        self.inputfiles += [filename]

class Pbes2boolDepthFirstTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbessolve'), settings)
        self.add_command_line_options('t2', ['-zdepth-first'])
        self.add_command_line_options('t3', ['-zdepth-first'])
        self.add_command_line_options('t4', ['-zdepth-first'])
        self.add_command_line_options('t5', ['-zdepth-first'])
        self.add_command_line_options('t6', ['-zdepth-first'])
        self.add_command_line_options('t7', ['-zdepth-first'])
        self.add_command_line_options('t8', ['-zdepth-first'])

class Pbes2bool_counter_exampleTest(ProcessTest):
    def __init__(self, name, optimization, settings):
        super().__init__(name, ymlfile('pbessolve-counter-example'), settings)
        if optimization in [4, 5]:
            self.add_command_line_options('t3', [f'-l{optimization}', '--aggressive', '--prune-todo-list'])
        else:
            self.add_command_line_options('t3', [f'-l{optimization}', '--prune-todo-list'])

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        filename = f'{self.name}.mcf'
        formula = random_state_formula_generator.make_modal_formula()
        write_text(filename, str(formula))
        self.inputfiles += [filename]

class Pbes2bool_evidence_ltsTest(ProcessTest):
    def __init__(self, name, optimization, settings):
        super().__init__(name, ymlfile('pbessolve-evidence-lts'), settings)
        if optimization in [4, 5]:
            self.add_command_line_options('t3', [f'-l{optimization}', '--aggressive', '--prune-todo-list'])
            self.add_command_line_options('t4', [f'-l{optimization}', '--aggressive', '--prune-todo-list'])
        else:
            self.add_command_line_options('t3', [f'-l{optimization}', '--prune-todo-list'])
            self.add_command_line_options('t4', [f'-l{optimization}', '--prune-todo-list'])

    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        filename = f'{self.name}.mcf'
        formula = random_state_formula_generator.make_modal_formula()
        write_text(filename, str(formula))
        self.inputfiles += [filename]

class Pbes2bool_counter_example_parelmTest(ProcessTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbessolve-parelm'), settings)


    def create_inputfiles(self, runpath = '.'):
        super().create_inputfiles(runpath)
        filename = f'{self.name}.mcf'
        formula = random_state_formula_generator.make_modal_formula()
        write_text(filename, str(formula))
        self.inputfiles += [filename]

class Pbes_unify_parametersTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbes-unify-parameters'), settings)

class Pbes_srfTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('pbes-srf'), settings)

# N.B does not work due to unknown expressions (F_or)
class SymbolicExplorationTest(PbesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('symbolic_exploration'), settings)

class BesTest(RandomTest):
    def __init__(self, name, ymlfile, settings):
        super().__init__(name, ymlfile, settings)
        self.equation_count = 4
        self.term_size = 3

    def create_inputfiles(self, runpath = '.'):
        filename = '{0}.txt'.format(self.name)
        p = make_bes(self.equation_count, self.term_size)
        write_text(filename, str(p))
        self.inputfiles += [filename]

class BessolveTest(BesTest):
    def __init__(self, name, settings):
        super().__init__(name, ymlfile('bessolve'), settings)

available_tests = {
    'alphabet-reduce'                             : lambda name, settings: AlphabetReduceTest(name, settings)                                          ,
    'lpssuminst'                                  : lambda name, settings: LpsSuminstTest(name, settings)                                              ,
    'lpssumelm'                                   : lambda name, settings: LpsSumelmTest(name, settings)                                               ,
    'lpsparelm'                                   : lambda name, settings: LpsParelmTest(name, settings)                                               ,
    'lps-quantifier-one-point'                    : lambda name, settings: LpsOnePointRuleRewriteTest(name, settings)                                  ,
    'lpsconfcheck-commutative'                    : lambda name, settings: LpsConfcheckTest(name, 'commutative', settings)                             ,
    'lpsconfcheck-commutative-disjoint'           : lambda name, settings: LpsConfcheckTest(name, 'commutative-disjoint', settings)                    ,
    'lpsconfcheck-disjoint'                       : lambda name, settings: LpsConfcheckTest(name, 'disjoint', settings)                                ,
    'lpsconfcheck-triangular'                     : lambda name, settings: LpsConfcheckTest(name, 'triangular', settings)                              ,
    'lpsconfcheck-trivial'                        : lambda name, settings: LpsConfcheckTest(name, 'trivial', settings)                                 ,
    'lpsconstelm'                                 : lambda name, settings: LpsConstelmTest(name, settings)                                             ,
    'lpsbinary'                                   : lambda name, settings: LpsBinaryTest(name, settings)                                               ,
    'lps2lts-algorithms'                          : lambda name, settings: Lps2ltsAlgorithmsTest(name, settings)                                       ,
    'lps2lts-parallel'                            : lambda name, settings: Lps2ltsParallelTest(name, settings)                                         ,
    'lps2pbes'                                    : lambda name, settings: Lps2pbesTest(name, settings)                                                ,
    'lps2pres'                                    : lambda name, settings: Lps2presTest(name, settings)                                                ,
    'lpsstategraph'                               : lambda name, settings: LpsstategraphTest(name, settings)                                           ,
    'lts2pbes'                                    : lambda name, settings: Lts2pbesTest(name, settings)                                                ,
    'lpsrewr'                                     : lambda name, settings: LpsrewrTest(name, settings)                                                ,
    'ltscompare-bisim'                            : lambda name, settings: LtscompareTest(name, 'bisim', settings)                                     ,
    'ltscompare-bisim-counter-example'            : lambda name, settings: LtscompareCounterexampleTest(name, 'bisim', False, settings)                ,
    'ltscompare-bisim-counter-example-hidden'     : lambda name, settings: LtscompareCounterexampleTest(name, 'bisim', True, settings)                 ,
    'ltscompare-bisim-gv'                         : lambda name, settings: LtscompareTest(name, 'bisim-gv', settings)                                  ,
    'ltscompare-bisim-gjkw'                       : lambda name, settings: LtscompareTest(name, 'bisim-gjkw', settings)                                ,
    'ltscompare-bisim-gj'                         : lambda name, settings: LtscompareTest(name, 'bisim-gj', settings)                                  ,
    'ltscompare-bisim-gj-lazy-BLC'                : lambda name, settings: LtscompareTest(name, 'bisim-gj-lazy-BLC', settings)                         ,
    'ltscompare-branching-bisim'                  : lambda name, settings: LtscompareTest(name, 'branching-bisim', settings)                           ,
    'ltscompare-branching-bisim-counter-example'  : lambda name, settings: LtscompareCounterexampleTest(name, 'branching-bisim', False, settings)      ,
    'ltscompare-branching-bisim-counter-example-hidden' : lambda name, settings: LtscompareCounterexampleTest(name, 'branching-bisim', True, settings) ,
    'ltscompare-branching-bisim-gv'               : lambda name, settings: LtscompareTest(name, 'branching-bisim-gv', settings)                        ,
    'ltscompare-branching-bisim-gjkw'             : lambda name, settings: LtscompareTest(name, 'branching-bisim-gjkw', settings)                      ,
    'ltscompare-branching-bisim-gj'               : lambda name, settings: LtscompareTest(name, 'branching-bisim-gj', settings)                        ,
    'ltscompare-branching-bisim-gj-lazy-BLC'      : lambda name, settings: LtscompareTest(name, 'branching-bisim-gj-lazy-BLC', settings)               ,
    'ltscompare-dpbranching-bisim'                : lambda name, settings: LtscompareTest(name, 'dpbranching-bisim', settings)                         ,
    'ltscompare-dpbranching-bisim-gv'             : lambda name, settings: LtscompareTest(name, 'dpbranching-bisim-gv', settings)                      ,
    'ltscompare-dpbranching-bisim-gjkw'           : lambda name, settings: LtscompareTest(name, 'dpbranching-bisim-gjkw', settings)                    ,
    'ltscompare-dpbranching-bisim-gj'             : lambda name, settings: LtscompareTest(name, 'dpbranching-bisim-gj', settings)                      ,
    'ltscompare-dpbranching-bisim-gj-lazy-BLC'    : lambda name, settings: LtscompareTest(name, 'dpbranching-bisim-gj-lazy-BLC', settings)             ,
    'ltscompare-weak-bisim'                       : lambda name, settings: LtscompareTest(name, 'weak-bisim', settings)                                ,
    'ltscompare-dpweak-bisim'                     : lambda name, settings: LtscompareTest(name, 'dpweak-bisim', settings)                              ,
    'ltscompare-sim'                              : lambda name, settings: LtscompareTest(name, 'sim', settings)                                       ,
    'ltscompare-ready-sim'                        : lambda name, settings: LtscompareTest(name, 'ready-sim', settings)                                 ,
    'ltscompare-trace'                            : lambda name, settings: LtscompareTest(name, 'trace', settings)                                     ,
    'ltscompare-trace-counter-example'            : lambda name, settings: LtscompareCounterexampleTest(name, 'trace', False, settings)                ,
    'ltscompare-trace-counter-example-hidden'     : lambda name, settings: LtscompareCounterexampleTest(name, 'trace', True, settings)                 ,
    'ltscompare-impossible-futures-counter-example' : lambda name, settings: LtscompareCounterexampleTest(name, 'impossible-futures', True, settings)  ,
    'ltscompare-weak-trace'                       : lambda name, settings: LtscompareTest(name, 'weak-trace', settings)                                ,
    'bisimulation-bisim'                          : lambda name, settings: BisimulationTest(name, 'bisim', settings)                                   ,
    'bisimulation-bisim-gv'                       : lambda name, settings: BisimulationTest(name, 'bisim-gv', settings)                                ,
    'bisimulation-bisim-gjkw'                     : lambda name, settings: BisimulationTest(name, 'bisim-gjkw', settings)                              ,
    'bisimulation-bisim-gj'                       : lambda name, settings: BisimulationTest(name, 'bisim-gj', settings)                                ,
    'bisimulation-bisim-gj-lazy-BLC'              : lambda name, settings: BisimulationTest(name, 'bisim-gj-lazy-BLC', settings)                       ,
    'bisimulation-branching-bisim'                : lambda name, settings: BisimulationTest(name, 'branching-bisim', settings)                         ,
    'bisimulation-branching-bisim-gv'             : lambda name, settings: BisimulationTest(name, 'branching-bisim-gv', settings)                      ,
    'bisimulation-branching-bisim-gjkw'           : lambda name, settings: BisimulationTest(name, 'branching-bisim-gjkw', settings)                    ,
    'bisimulation-branching-bisim-gj'             : lambda name, settings: BisimulationTest(name, 'branching-bisim-gj', settings)                      ,
    'bisimulation-branching-bisim-gj-lazy-BLC'    : lambda name, settings: BisimulationTest(name, 'branching-bisim-gj-lazy-BLC', settings)             ,
    'bisimulation-weak-bisim'                     : lambda name, settings: BisimulationTest(name, 'weak-bisim', settings)                              ,
    'pbesconstelm'                                : lambda name, settings: PbesconstelmTest(name, settings)                                            ,
    'pbesrewr-simplify'                           : lambda name, settings: PbesrewrTest(name, 'simplify', settings)                                    ,
    'pbesrewr-pfnf'                               : lambda name, settings: PbesrewrTest(name, 'pfnf', settings)                                        ,
    'pbesrewr-quantifier-all'                     : lambda name, settings: PbesrewrTest(name, 'quantifier-all', settings)                              ,
    'pbesrewr-quantifier-finite'                  : lambda name, settings: PbesrewrTest(name, 'quantifier-finite', settings)                           ,
    'pbesrewr-quantifier-inside'                  : lambda name, settings: PbesrewrTest(name, 'quantifier-inside', settings)                           ,
    'pbesrewr-quantifier-one-point'               : lambda name, settings: PbesrewrTest(name, 'quantifier-one-point', settings)                        ,
    'pbesrewr-data-rewriter'                      : lambda name, settings: PbestransformTest(name, 'pbes-data-rewriter', settings)                     ,
    'pbesrewr-simplify-rewriter'                  : lambda name, settings: PbestransformTest(name, 'pbes-simplify-rewriter', settings)                 ,
    'pbesrewr-simplify-data-rewriter'             : lambda name, settings: PbestransformTest(name, 'pbes-simplify-data-rewriter', settings)            ,
    'pbesrewr-simplify-quantifiers-rewriter'      : lambda name, settings: PbestransformTest(name, 'pbes-simplify-quantifiers-rewriter', settings)     ,
    'pbesrewr-simplify-quantifiers-data-rewriter' : lambda name, settings: PbestransformTest(name, 'pbes-simplify-quantifiers-data-rewriter', settings),
    'pbesinst-lazy'                               : lambda name, settings: PbesinstTest(name, ['-slazy'], settings)                                    ,
    'pbesinst-alternative_lazy'                   : lambda name, settings: PbesinstTest(name, ['-salternative-lazy'], settings)                        ,
    'pbesinst-finite'                             : lambda name, settings: PbesinstTest(name, ['-sfinite', '-f*(*:Bool)'], settings)                   ,
    'pbespgsolve'                                 : lambda name, settings: PbespgsolveTest(name, settings)                                             ,
    'pbessolve'                                   : lambda name, settings: Pbes2boolTest(name, False, settings)                                        ,
    'pbessolve-parallel'                          : lambda name, settings: Pbes2boolTest(name, True, settings)                                         ,
    'pbessolve-depth-first'                       : lambda name, settings: Pbes2boolDepthFirstTest(name, settings)                                     ,
    'pbessolve-counter-example-optimization-0'    : lambda name, settings: Pbes2bool_counter_exampleTest(name, 0, settings)                            ,
    'pbessolve-counter-example-optimization-1'    : lambda name, settings: Pbes2bool_counter_exampleTest(name, 1, settings)                            ,
    'pbessolve-counter-example-optimization-2'    : lambda name, settings: Pbes2bool_counter_exampleTest(name, 2, settings)                            ,
    'pbessolve-counter-example-optimization-3'    : lambda name, settings: Pbes2bool_counter_exampleTest(name, 3, settings)                            ,
    'pbessolve-counter-example-optimization-4'    : lambda name, settings: Pbes2bool_counter_exampleTest(name, 4, settings)                            ,
    'pbessolve-counter-example-optimization-5'    : lambda name, settings: Pbes2bool_counter_exampleTest(name, 5, settings)                            ,
    'pbessolve-counter-example-optimization-6'    : lambda name, settings: Pbes2bool_counter_exampleTest(name, 6, settings)                            ,
    'pbessolve-counter-example-optimization-7'    : lambda name, settings: Pbes2bool_counter_exampleTest(name, 7, settings)                            ,
    'pbessolve-evidence-lts-optimization-0'       : lambda name, settings: Pbes2bool_evidence_ltsTest(name, 0, settings)                               ,
    'pbessolve-evidence-lts-optimization-1'       : lambda name, settings: Pbes2bool_evidence_ltsTest(name, 1, settings)                               ,
    'pbessolve-evidence-lts-optimization-2'       : lambda name, settings: Pbes2bool_evidence_ltsTest(name, 2, settings)                               ,
    'pbessolve-evidence-lts-optimization-3'       : lambda name, settings: Pbes2bool_evidence_ltsTest(name, 3, settings)                               ,
    'pbessolve-evidence-lts-optimization-4'       : lambda name, settings: Pbes2bool_evidence_ltsTest(name, 4, settings)                               ,
    'pbessolve-evidence-lts-optimization-5'       : lambda name, settings: Pbes2bool_evidence_ltsTest(name, 5, settings)                               ,
    'pbessolve-evidence-lts-optimization-6'       : lambda name, settings: Pbes2bool_evidence_ltsTest(name, 6, settings)                               ,
    'pbessolve-evidence-lts-optimization-7'       : lambda name, settings: Pbes2bool_evidence_ltsTest(name, 7, settings)                               ,
    'pbessolve-parelm'                            : lambda name, settings: Pbes2bool_counter_example_parelmTest(name, settings)                        ,
    'pbesstategraph'                              : lambda name, settings: PbesstategraphTest(name, settings)                                          ,
    'pbes-unify-parameters'                       : lambda name, settings: Pbes_unify_parametersTest(name, settings)                                   ,
    'pbes-srf'                                    : lambda name, settings: Pbes_srfTest(name, settings)                                                ,
    'stochastic-ltscompare'                       : lambda name, settings: StochasticLtscompareTest(name, settings)                                    ,
}

available_experimental_tests = {
    # ML: Disable until the tests pass.
    'ltscombine'                                  : lambda name, settings: LtsCombineTest(name, False, settings)                                       ,
    'ltscombine-parallel'                         : lambda name, settings: LtsCombineTest(name, True, settings)                                        ,
    'pbesparelm'                                  : lambda name, settings: PbesparelmTest(name, settings)                                              ,
    'pressolve'                                   : lambda name, settings: Pres2boolTest(name, settings)                                               ,
    'pbespor2'                                    : lambda name, settings: Pbespor2Test(name, settings)                                                ,
    'bessolve'                                    : lambda name, settings: BessolveTest(name, settings)                                                ,
    'pbeschain'                                   : lambda name, settings: PbeschainTest(name, settings)                                               ,
    'pbesfixpointsolve'                           : lambda name, settings: PbesfixpointsolveTest(name, settings)                                       ,
    'pbespareqelm'                                : lambda name, settings: PbespareqelmTest(name, settings)
}

if shutil.which("z3") is not None:
    available_tests.update({'pbessymbolicbisim' : lambda name, settings: PbessymbolicbisimTest(name, settings)})

# These test do not work on Windows due to dependencies.
if os.name != 'nt':
    available_tests.update({'pbessolvesymbolic' : lambda name, settings: PbessolvesymbolicTest(name, [], settings) })
    available_tests.update({'pbessolvesymbolic-parallel' : lambda name, settings: PbessolvesymbolicTest(name, ['--threads=8'], settings) })
    available_tests.update({'pbessolvesymbolic-total' : lambda name, settings: PbessolvesymbolicTest(name, ['--total'], settings) })
    available_tests.update({'pbessolvesymbolic-chaining' : lambda name, settings: PbessolvesymbolicTest(name, ['--chaining'], settings) })
    available_tests.update({'pbessolvesymbolic-total-chaining' : lambda name, settings: PbessolvesymbolicTest(name, ['--total', '--chaining'], settings) })
    available_tests.update({'pbessolvesymbolic-saturation' : lambda name, settings: PbessolvesymbolicTest(name, ['--total', '--saturation'], settings) })
    available_tests.update({'pbessolvesymbolic-split' : lambda name, settings: PbessolvesymbolicTest(name, ['--total', '--split-conditions'], settings) })

    available_tests.update({'pbessolvesymbolic-partial-s1' : lambda name, settings: PbessolvesymbolicTest(name, ['-s1', '--aggressive'], settings) })
    available_tests.update({'pbessolvesymbolic-partial-s2' : lambda name, settings: PbessolvesymbolicTest(name, ['-s2', '--aggressive'], settings) })
    available_tests.update({'pbessolvesymbolic-partial-s3' : lambda name, settings: PbessolvesymbolicTest(name, ['-s3', '--aggressive'], settings) })
    available_tests.update({'pbessolvesymbolic-partial-s4' : lambda name, settings: PbessolvesymbolicTest(name, ['-s4', '--aggressive'], settings) })
    available_tests.update({'pbessolvesymbolic-partial-s5' : lambda name, settings: PbessolvesymbolicTest(name, ['-s5', '--aggressive'], settings) })
    available_tests.update({'pbessolvesymbolic-partial-s6' : lambda name, settings: PbessolvesymbolicTest(name, ['-s6', '--aggressive'], settings) })
    available_tests.update({'pbessolvesymbolic-partial-s7' : lambda name, settings: PbessolvesymbolicTest(name, ['-s7', '--aggressive'], settings) })

    available_tests.update({'pbessolvesymbolic-counter-example' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, [], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-s1' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['-s1'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-s2' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['-s2'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-s3' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['-s3'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-s4' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['-s4'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-s7' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['-s7'], settings) })

    available_tests.update({'pbessolvesymbolic-counter-example-saturation' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--saturation'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-saturation-s1' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--saturation', '-s1'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-saturation-s2' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--saturation', '-s2'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-saturation-s3' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--saturation', '-s3'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-saturation-s4' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--saturation', '-s4'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-saturation-s7' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--saturation', '-s7'], settings) })

    available_tests.update({'pbessolvesymbolic-counter-example-chaining' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--chaining'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-chaining-s1' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--chaining', '-s1'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-chaining-s2' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--chaining', '-s2'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-chaining-s3' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--chaining', '-s3'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-chaining-s4' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--chaining', '-s4'], settings) })
    available_tests.update({'pbessolvesymbolic-counter-example-chaining-s7' : lambda name, settings: PbessolvesymbolicCounterexampleTest(name, ['--chaining', '-s7'], settings) })
    # This test fails because of odr-violation.
    # available_tests.update({'ltsconvertsymbolic' :  lambda name, settings: LtsconvertsymbolicTest(name, [], settings)})
    # This test fails because of data races.
    # available_tests.update({'ltsconvertsymbolic-parallel' : lambda name, settings: LtsconvertsymbolicTest(name, ['--threads=8'], settings) })

def print_names(tests):
    for name in sorted(tests):
        print(name)

# Return all tests that match with pattern. In case of an exact match, only this exact match is returned.
def matching_tests(tests, pattern):
    matches = [name for name in sorted(tests) if re.search(pattern, name)]
    if pattern in matches:
        return [pattern]
    return matches

def main():
    cmdline_parser = argparse.ArgumentParser()
    cmdline_parser.add_argument('-t', '--toolpath', dest='toolpath', help='The path where the mCRL2 tools are installed')
    cmdline_parser.add_argument('-r', '--repetitions', dest='repetitions', metavar='N', default='10', help='Perform N repetitions of each test')
    cmdline_parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Display additional progress messages.')
    cmdline_parser.add_argument('-k', '--keep-files', dest='keep_files', action='store_true', help='Keep the files produced by the test')
    cmdline_parser.add_argument('-n', '--names', dest='names', action='store_true', help='Print the names of the available tests')
    cmdline_parser.add_argument('-p', '--pattern', dest='pattern', metavar='P', default='.', action='store', help='Run the tests that match with pattern P')
    cmdline_parser.add_argument('-o', '--output', dest='output', metavar='o', action='store', help='Run the tests in the given directory')
    cmdline_parser.add_argument('-m', '--timeout', dest='timeout', metavar='t', action='store', help='Set the timeout for each test to t seconds')
    cmdline_parser.add_argument('-e', '--experimental', dest='experimental', action='store_true', help='Run random tests using experimental tools.')
    cmdline_parser.add_argument('-i', '--python', dest='python', action='store', help='Sets the path to the Python interpreter that is used.')

    args = cmdline_parser.parse_args()
    tests = available_tests
    if args.experimental:
        tests.update(available_experimental_tests)

    if args.names:
        print_names(tests)
        return

    python_path = None
    if args.python:
        # Check if the given path exists
        python_path = args.python
    else:
        python_path = shutil.which("python3")

    if python_path is None:
        print("Cannot find python in the PATH environment, and --python was not supplied")
        sys.exit(-1)

    if not args.toolpath:
        print("To run tests the --toolpath argument is required")
        sys.exit(-1)

    settings = {'toolpath': os.path.abspath(args.toolpath), 'verbose': args.verbose, 'cleanup_files': not args.keep_files, 'allow-non-zero-return-values': True, 'python_path': python_path}

    if args.timeout:
        settings['timeout'] = int(args.timeout)

    repeats = range(int(args.repetitions))

    if args.output:
        if not os.path.exists(args.output):
            os.mkdir(args.output)
        os.chdir(args.output)

    test_failed = False
    for name in matching_tests(tests, args.pattern):
        try:
            for i in repeats:
                test = tests[name](f'{name}_{i}', settings)
                test.execute_in_sandbox()
        except Exception as e:
            print('An exception occurred:', e.__class__, e)
            traceback.print_exc()
            test_failed = True
            sys.exit(-1)

    if test_failed:
        sys.exit(-1)

if __name__ == '__main__':
    main()
