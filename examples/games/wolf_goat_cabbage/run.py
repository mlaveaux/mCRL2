#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'wolf_goat_cabbage.mcrl2', 'wolf_goat_cabbage.lps'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'wolf_goat_cabbage.lps', 'wolf_goat_cabbage.nodeadlock.pbes'], check=True)

if '-rjittyc' in argv:
    subprocess.run(['lps2lts', '-vrjittyc', 'wolf_goat_cabbage.lps', 'wolf_goat_cabbage.lts'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'wolf_goat_cabbage.nodeadlock.pbes'], check=True)

