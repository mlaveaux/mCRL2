#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'clobber.mcrl2', 'clobber.lps'], check=True)

if "-rjittyc" in argv:
    subprocess.run(['lps2lts', '-vrjittyc', 'clobber.lps', 'clobber.aut'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'clobber.lps', 'clobber.nodeadlock.pbes'], check=True)

if "-rjittyc" in argv:
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'clobber.nodeadlock.pbes'], check=True)