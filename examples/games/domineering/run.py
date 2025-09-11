#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-vD', 'domineering.mcrl2', 'domineering.lps'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'domineering.lps', 'domineering.nodeadlock.pbes'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'player1_can_win.mcf', 'domineering.lps', 'domineering.player1_can_win.pbes'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'player2_can_win.mcf', 'domineering.lps', 'domineering.player2_can_win.pbes'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'eventually_player1_or_player2_wins.mcf', 'domineering.lps', 'domineering.eventually_player1_or_player2_wins.pbes'], check=True)

if "-rjittyc" in argv:
    subprocess.run(['lps2lts', '-vrjittyc', 'domineering.lps', 'domineering.lts'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'domineering.nodeadlock.pbes'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'domineering.player1_can_win.pbes'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'domineering.player2_can_win.pbes'], check=True)
    subprocess.run(['pbessolve', '-vrjittyc', '-s1', 'domineering.eventually_player1_or_player2_wins.pbes'], check=True)
