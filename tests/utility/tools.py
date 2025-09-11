#!/usr/bin/env python

# Copyright 2013, 2014 Mark Geelen.
# Copyright 2014, 2015 Wieger Wesselink.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

import os.path
import platform
import re

from .text_utility import read_text
from .run_process import RunProcess

def is_list_of(l, types):
    if not isinstance(l, list):
        return False
    for x in l:
        if not isinstance(x, types):
            return False
    return True

class Node:
    def __init__(self, label, _type, value):
        self.label = label
        self.type = _type
        self.value = value
        return

    def __str__(self):
        return 'Node(label = {0}, type = {1}, value = {2})'.format(self.label, self.type, self.value)

    def filename(self):
        return '{}.{}'.format(self.label, self.type)

class Tool(object):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        assert is_list_of(input_nodes, Node)
        assert is_list_of(output_nodes, Node)
        self.label = label
        self.name = name
        self.toolpath = toolpath
        self.input_nodes = input_nodes
        self.output_nodes = output_nodes
        self.args = args
        self.stderr = ""
        self.executed = False
        self.value = {}

    # If no_paths is True, then all paths in the command are excluded
    def arguments(self, working_directory = None, no_paths = False):
        if not working_directory:
            working_directory = os.getcwd()
        input_filenames = [node.filename() for node in self.input_nodes]
        output_filenames = [node.filename() for node in self.output_nodes if node.type != 'Bool']
        filenames = input_filenames + output_filenames
        if not no_paths:
            filenames = [os.path.join(working_directory, filename) for filename in filenames]
        return filenames

    def assign_outputs(self):
        for node in self.output_nodes:
            if node.type == 'text':
                text = read_text(node.filename())
                node.value = text
            else:
                node.value = 'executed'

    # value[key] is an integer
    def parse_number(self, text, key, regex):
        m = re.search(regex, text)
        if m is not None:
            self.value[key] = int(m.group(1))

    # value[key] is an integer
    def parse_numbers(self, text, key1, key2, regex):
        m = re.search(regex, text)
        if m is not None:
            self.value[key1] = int(m.group(1))
            self.value[key2] = int(m.group(2))

    # value[key] is a set of strings
    # All occurrences of regex are processed
    def parse_action(self, text, key, regex):
        for m in re.finditer(regex, text):
            if not key in self.value:
                self.value[key] = set([])
            self.value[key].add(m.group(1))

    # value[key] is a boolean
    # multiple regular expressions are checked
    def parse_boolean_regexes(self, text, key, regexes):
        result = False
        for regex in regexes:
            if re.search(regex, text, re.DOTALL) is not None:
                result = True
        self.value[key] = result

    # value[key] is a boolean
    def parse_boolean(self, text, key, regex, negated_regex = None):
        if negated_regex:
            m = re.search(negated_regex, text, re.DOTALL)
            if m is not None:
                self.value[key] = False
        if regex:
            m = re.search(regex, text, re.DOTALL)
            if m is not None:
                self.value[key] = True

    def parse_output(self):
        text = self.stdout + self.stderr
        self.parse_number(text, 'summand-count'               , r'Number of summands                  : (\d+)')
        self.parse_number(text, 'tau-summand-count'           , r'Number of tau-summands              : (\d+)')
        self.parse_number(text, 'global-variable-count'       , r'Number of declared global variables : (\d+)')
        self.parse_number(text, 'process-parameter-count'     , r'Number of process parameters        : (\d+)')
        self.parse_number(text, 'action-label-count'          , r'Number of declared action labels    : (\d+)')
        self.parse_number(text, 'used-action-label-count'     , r'Number of used actions              : (\d+)')
        self.parse_number(text, 'used-multi-action-count'     , r'Number of used multi-actions        : (\d+)')
        self.parse_number(text, 'state-count'                 , r'Number of states: (\d+)')
        self.parse_number(text, 'state-label-count'           , r'Number of state labels: (\d+)')
        self.parse_number(text, 'action-label-count'          , r'Number of action labels: (\d+)')
        self.parse_number(text, 'transition-count'            , r'Number of transitions: (\d+)')
        self.parse_number(text, 'equation-count'              , r'Number of equations: (\d+)')
        self.parse_number(text, 'mu-count'                    , r"Number of mu's:      (\d+)")
        self.parse_number(text, 'nu-count'                    , r"Number of nu's:      (\d+)")
        self.parse_number(text, 'block-nesting-depth'         , r'Block nesting depth: (\d+)')
        self.parse_number(text, 'vertex-count'                , r'Number of vertices in the structure graph: (\d+)')
        self.parse_numbers(text, 'confluent-tau-summand-count', 'tau-summand-count', r'(\d+) of (\d+) tau summands were found to be confluent')
        self.parse_boolean(text, 'has-state-labels'           , 'Has state labels.', 'Does not have state labels.')
        self.parse_boolean(text, 'has-action-labels'          , 'Has action labels.')
        self.parse_boolean(text, 'is-deterministic'           , 'LTS is deterministic.', 'LTS is not deterministic.')
        self.parse_boolean(text, 'is-closed'                  , 'is closed', 'is not closed')
        self.parse_boolean(text, 'is-well-formed'             , 'well formed', 'not well formed')
        self.parse_boolean(text, 'is-well-typed'              , 'is well typed', 'is not well typed')
        self.parse_boolean(text, 'result'                     , r'LTSs are strongly bisimilar', 'LTSs are not strongly bisimilar')
        self.parse_boolean(text, 'result'                     , r'LTSs are branching bisimilar', 'LTSs are not branching bisimilar')
        self.parse_boolean(text, 'result'                     , r'LTSs are equal \(branching bisimilarity using the almost-O\(m log n\) Groote/Wijs algorithm\)', r'LTSs are not equal \(branching bisimilarity using the almost-O\(m log n\) Groote/Wijs algorithm\)')
        self.parse_boolean(text, 'result'                     , r'LTSs are equal \(branching bisimilarity using the O\(m log n\) Groote/Keiren/Jansen/Wijs algorithm\)', r'LTSs are not equal \(branching bisimilarity using the O\(m log n\) Groote/Keiren/Jansen/Wijs algorithm\)')
        self.parse_boolean(text, 'result'                     , r'LTSs are divergence preserving branching bisimilar', 'LTSs are not divergence preserving branching bisimilar')
        self.parse_boolean(text, 'result'                     , r'LTSs are weak bisimilar', 'LTSs are not weak bisimilar')
        self.parse_boolean(text, 'result'                     , r'LTSs are divergence preserving weak bisimilar', 'LTSs are not divergence preserving weak bisimilar')
        self.parse_boolean(text, 'result'                     , r'LTSs are strongly simulation equivalent', 'LTSs are not strongly simulation equivalent')
        self.parse_boolean(text, 'result'                     , r'LTSs are strongly trace equivalent', 'LTSs are not strongly trace equivalent')
        self.parse_boolean(text, 'result'                     , r'LTSs are weak trace equivalent', 'LTSs are not weak trace equivalent')
        self.parse_boolean(text, 'result'                     , r'LTSs are equal', 'LTSs are not equal')
        self.parse_boolean(text, 'result'                     , r'is included in', 'is not included in')
        self.parse_action(text, 'actions'                     , r"Detected action '(\w+)'")
        self.parse_action(text, 'actions'                     , r"Action '(\w+)' found")
        self.parse_boolean_regexes(text, 'has-deadlock'       , [r'deadlock-detect: deadlock found', r'Deadlock found'])
        self.parse_boolean_regexes(text, 'has-divergence'     , [r'divergence-detect: divergence found', r'Divergent state found'])
        self.parse_boolean(text, 'has-nondeterminism'         , r'Nondeterministic state found')

    # If no_paths is True, then all paths in the command are excluded
    def command(self, working_directory = None, no_paths = False):
        args = self.arguments(working_directory, no_paths)
        name = self.name
        if not no_paths:
            name = os.path.join(self.toolpath, name)
        return ' '.join([name] + args + self.args)


    def execute(self, timeout, memlimit, verbose):
        args = []
        if self.name.endswith('.py'):
            # If it is a python script then the toolpath is the Python interpreter.
            args = [self.name]
            name = self.toolpath
        else:
            name = os.path.join(self.toolpath, self.name)
        if verbose:
            print('Executing ' + ' '.join([name] + args + self.args))
            
        args.extend(self.arguments())

        process = RunProcess(name, args + self.args, memlimit, timeout)
   
        self.executed = True
        self.stdout = process.stdout
        self.stderr = process.stderr
        self.assign_outputs()
        self.parse_output()
        return process.returncode

    def __str__(self):
        import io
        out = io.StringIO()
        out.write('label    = ' + str(self.label)    + '\n')
        out.write('name     = ' + str(self.name)     + '\n')
        out.write('input    = [{0}]\n'.format(', '.join([str(x) for x in self.input_nodes])))
        out.write('output   = [{0}]\n'.format(', '.join([str(x) for x in self.output_nodes])))
        out.write('args     = ' + str(self.args)     + '\n')
        out.write('stderr   = ' + str(self.stderr)    + '\n')
        out.write('executed = ' + str(self.executed) + '\n')
        return out.getvalue()

class SolveTool(Tool):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        super(SolveTool, self).__init__(label, name, toolpath, input_nodes, output_nodes, args)

    def assign_outputs(self):
        text = self.stdout.strip() + self.stderr.strip()
        if text.endswith('true'):
            value = True
        elif text.endswith('false'):
            value = False
        else:
            value = None
        self.value['solution'] = value

class Lps2PbesTool(Tool):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        assert len(input_nodes) == 2
        assert len(output_nodes) == 1
        super(Lps2PbesTool, self).__init__(label, name, toolpath, input_nodes, output_nodes, args)

    def arguments(self, working_directory = None, no_paths = False):
        args = super(Lps2PbesTool, self).arguments(working_directory, no_paths)
        args.insert(1, '-f')
        return args

class Lts2PbesTool(Tool):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        assert len(input_nodes) == 2
        assert len(output_nodes) == 1
        super(Lts2PbesTool, self).__init__(label, name, toolpath, input_nodes, output_nodes, args)

    def arguments(self, working_directory = None, no_paths = False):
        args = super(Lts2PbesTool, self).arguments(working_directory, no_paths)
        args.insert(1, '-f')
        return args

class Lts2LpsTool(Tool):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        assert len(input_nodes) == 2
        assert len(output_nodes) == 1
        super(Lts2LpsTool, self).__init__(label, name, toolpath, input_nodes, output_nodes, args)

    def arguments(self, working_directory = None, no_paths = False):
        args = super(Lts2LpsTool, self).arguments(working_directory, no_paths)
        args.insert(1, '-l')
        return args

class Lps2LtsTool(Tool):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        super(Lps2LtsTool, self).__init__(label, name, toolpath, input_nodes, output_nodes, args)

    def assign_outputs(self):
        self.value['has-deadlock'] = None
        self.value['has-nondeterminism'] = None
        self.value['has-divergence'] = None
        self.value['actions'] = set([])
        super(Lps2LtsTool, self).assign_outputs()

class PbesSolveTool(Tool):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        super(PbesSolveTool, self).__init__(label, name, toolpath, input_nodes, output_nodes, args)

    def arguments(self, working_directory = None, no_paths = False):
        args = super(PbesSolveTool, self).arguments(working_directory, no_paths)

        # counter example generation
        if len(self.input_nodes) > 1:
            args[1] = '--file={}'.format(args[1])
        if len(self.output_nodes) > 0:
            args[2] = '--evidence-file={}'.format(args[2])
        return args

    def assign_outputs(self):
        text = self.stdout.strip() + self.stderr.strip()
        # N.B. In verbose mode, the solution may appear at the start.
        if text.startswith('true') or text.endswith('true'):
            value = True
        elif text.startswith('false') or text.endswith('false'):
            value = False
        else:
            value = None
        self.value['solution'] = value

        # mark the evidence file as executed
        if len(self.output_nodes) == 1:
            self.output_nodes[0].value = 'executed'
  
class PresSolveTool(Tool):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        super(PresSolveTool, self).__init__(label, name, toolpath, input_nodes, output_nodes, args)

    def assign_outputs(self):
        text = self.stdout.strip() + self.stderr.strip()
        # N.B. In verbose mode, the solution may appear at the start.
        if text.startswith('false') or text.endswith('false') or text.startswith('-inf') or text.endswith('-inf'):
            value = False
        elif text.startswith('true') or text.endswith('true') or text.startswith('inf') or text.endswith('inf'):
            value = True
        else:
            value = self.stdout.strip()
        self.value['solution'] = value

class LtscompareTool(Tool):
    def __init__(self, label, name, toolpath, input_nodes, output_nodes, args):
        super(LtscompareTool, self).__init__(label, name, toolpath, input_nodes, output_nodes, args)

    def arguments(self, working_directory = None, no_paths = False):
        args = super(LtscompareTool, self).arguments(working_directory, no_paths)

        # counter example generation
        if len(self.output_nodes) > 0:
            args[2] = '--counter-example-file={}'.format(args[2])
        return args

class ToolFactory(object):
    def create_tool(self, label, name, toolpath, input_nodes, output_nodes, args, python_path):
        if name in ['lps2pbes', 'lps2pres']:
            return Lps2PbesTool(label, name, toolpath, input_nodes, output_nodes, args)
        elif name == 'lts2pbes':
            return Lts2PbesTool(label, name, toolpath, input_nodes, output_nodes, args)
        elif name in ['generatelts', 'lps2lts']:
            return Lps2LtsTool(label, name, toolpath, input_nodes, output_nodes, args)
        elif name == 'lts2lps':
            return Lts2LpsTool(label, name, toolpath, input_nodes, output_nodes, args)
        elif name in ['pbespgsolve', 'bessolve']:
            return SolveTool(label, name, toolpath, input_nodes, output_nodes, args)
        elif name in ['pbessolve', 'pbessolve', 'pbessolvesymbolic', 'pbessymbolicbisim']:
            return PbesSolveTool(label, name, toolpath, input_nodes, output_nodes, args)
        elif name in ['pressolve']:
            return PresSolveTool(label, name, toolpath, input_nodes, output_nodes, args)
        elif name == 'ltscompare':
            return LtscompareTool(label, name, toolpath, input_nodes, output_nodes, args)
        elif name == 'ltscombine_naive.py':
            return Tool(label, os.path.join(os.path.abspath(os.path.dirname(__file__)), '../scripts/', name), python_path, input_nodes, output_nodes, [] + args)
    
        return Tool(label, name, toolpath, input_nodes, output_nodes, args)
