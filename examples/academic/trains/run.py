#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'trains.mcrl2', 'trains.lps'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'trains.lps', 'trains.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'trains.nodeadlock.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'infinitely_often_enabled_then_infinitely_often_taken_enter.mcf', 'trains.lps', 'trains.infinitely_often_enabled_then_infinitely_often_taken_enter.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'trains.infinitely_often_enabled_then_infinitely_often_taken_enter.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'mutual_exclusion.mcf', 'trains.lps', 'trains.mutual_exclusion.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'trains.mutual_exclusion.pbes'], check=True)
