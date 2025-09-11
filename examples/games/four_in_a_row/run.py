#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

run = subprocess.run(['mcrl22lps', '-vDn', 'four_in_a_row.mcrl2'], stdout=subprocess.PIPE, check=True)

# We use sum instantiation and parunfold to split the list of lists for the board into individual simple parameters.
run =  subprocess.run(['lpssuminst', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsparunfold', '-sBoard', '-n5'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsparunfold', '-sRow', '-n5'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsrewr', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
subprocess.run(['lpsconstelm', '-v', '-', 'temp.lps'], input=run.stdout, check=True)

if '-rjittyc' in argv:
    # The state space of the original 7x6 four in a row is so big that explicit
    # state space exploration will not terminate. Reducing the number of columns,
    # e.g. to four instead of default 7 does lead to a result. In that case the
    # there are 16M states. Solving the formulas red_wins.mcf requires 7M bes
    # variables (with result false). Five columns lead to appr. 100Mstates. See
    # four_in_a_row_symbolic for an approach using the symbolic tools that are able
    # to solve the full game on sufficient hardware (approx 1TB of memory).
    subprocess.run(['lps2lts', '-v', '-rjittyc', 'temp.lps'], check=True)

# Also pbessolve does not tend to terminate on the full game, as it is too big.
# os.system('lps2pbes -vfred_wins.mcf temp2.lps temp.pbes')
# os.system('pbesrewr -pquantifier-all temp.pbes | pbesconstelm -v | pbesrewr -v > temp1.pbes')
# os.system('pbessolve -rjittyc -v -s2 -zdepth-first temp1.pbes ')#



