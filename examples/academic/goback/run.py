#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', 'goback.mcrl2', 'goback.lps'], check=True)

print('Checking property nodeadlock.mcf for goback.mcrl2')
subprocess.run(['lps2pbes', '-f', 'nodeadlock.mcf', 'goback.lps', 'goback.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', 'goback.nodeadlock.pbes'], check=True)