#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'chatbox.mcrl2', 'chatbox.lps'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'chatbox.lps', 'chatbox.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'chatbox.nodeadlock.pbes'], check=True)

subprocess.run(['lps2lts', '-v', 'chatbox.lps', 'chatbox.aut'], check=True)

