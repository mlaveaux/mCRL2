#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-vnf', 'parallel_counting.mcrl2', 'temp.lps'], check=True)
subprocess.run(['lps2pbes', '-v', '-fparallel_counting.mcf', 'temp.lps', 'temp.pbes'], check=True)
subprocess.run(['pbessolve', '-v', '-s2', 'temp.pbes'], check=True)
