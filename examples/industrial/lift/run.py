#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-vn', 'lift3-final.mcrl2', 'lift3-final.lps'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'lift3-final.lps', 'lift3-final.nodeadlock.pbes'], check=True)

subprocess.run(['mcrl22lps', '-vn', 'lift3-init.mcrl2', 'lift3-init.lps'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'lift3-init.lps', 'lift3-init.nodeadlock.pbes'], check=True)

if '-rjittyc' in argv:
    subprocess.run(['lps2lts', '-vrjittyc', 'lift3-final.lps', 'lift3-final.aut'], check=True)
    subprocess.run(['lps2lts', '-vrjittyc', 'lift3-init.lps', 'lift3-init.aut'], check=True)

    subprocess.run(['ltscompare', '-v', '-ebisim', 'lift3-init.aut', 'lift3-final.aut'], check=True)
    subprocess.run(['lpsbisim2pbes', '-v', '-bstrong-bisim', 'lift3-init.lps', 'lift3-final.lps', 'lift3-bisim.pbes'], check=True)

    subprocess.run(['pbessolve', '-vrjittyc', 'lift3-final.nodeadlock.pbes'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', 'lift3-init.nodeadlock.pbes'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', 'lift3-bisim.pbes'], check=True)