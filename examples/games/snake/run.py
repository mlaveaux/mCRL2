#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-vD', 'snake.mcrl2', 'snake.lps'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'snake.lps', 'snake.nodeadlock.pbes'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'white_can_win.mcf', 'snake.lps', 'snake.white_can_win.pbes'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'black_can_win.mcf', 'snake.lps', 'snake.black_can_win.pbes'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'eventually_white_or_black_wins.mcf', 'snake.lps', 'snake.eventually_white_or_black_wins.pbes'], check=True)

if '-rjittyc' in argv:
    subprocess.run(['lps2lts', '-vrjittyc', 'snake.lps', 'snake.aut'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'snake.nodeadlock.pbes'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'snake.white_can_win.pbes'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'snake.black_can_win.pbes'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'snake.eventually_white_or_black_wins.pbes'], check=True)
