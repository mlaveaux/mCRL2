#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'alma.mcrl2', 'alma.lps'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'alma.lps', 'alma.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'alma.nodeadlock.pbes'], check=True)

subprocess.run(['lps2lts', '-v', 'alma.lps', 'alma.aut'], check=True)

