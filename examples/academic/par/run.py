#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'par.mcrl2', 'par.lps'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'par.lps', 'par.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'par.nodeadlock.pbes'], check=True)

# The following creates state space, and stores trace to deadlock.
subprocess.run(['lps2lts', '-v', '-Dt', 'par.lps', 'par.aut'], check=True)
