#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', '1394-fin.mcrl2', '1394-fin.lps'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', '1394-fin.lps', '1394-fin.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', '-v', '1394-fin.nodeadlock.pbes'], check=True)

subprocess.run(['lps2lts', '-v', '--cached', '1394-fin.lps', '1394-fin.aut'], check=True)

