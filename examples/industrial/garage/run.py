#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

run = subprocess.run(['mcrl22lps', '-v', 'garage-ver.mcrl2'], stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpssuminst', '-v', '-f'], input=run.stdout, stdout=subprocess.PIPE, check=True)
subprocess.run(['lpssumelm', '-v', '-', 'garage-ver1.lps'], input=run.stdout, check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'garage-ver1.lps', 'garage-ver.nodeadlock.pbes'], check=True)

if '-rjittyc' in argv:
    subprocess.run(['lps2lts', '-v', '-rjittyc', 'garage-ver1.lps'], check=True)
    subprocess.run(['pbessolve', '-v', '-rjittyc', 'garage-ver.nodeadlock.pbes'], check=True)
