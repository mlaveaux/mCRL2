#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'brp.mcrl2', 'brp.lps'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'brp.lps', 'brp.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'brp.nodeadlock.pbes'], check=True)

subprocess.run(['lps2lts', '-v', 'brp.lps', 'brp.aut'], check=True)

