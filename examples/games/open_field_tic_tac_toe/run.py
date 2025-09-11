#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

run = subprocess.run(['mcrl22lps', '-v', 'open_field_tictactoe.mcrl2'], stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpssuminst', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsparunfold', '-sBoard', '-n6'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsparunfold', '-sRow', '-n8'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsrewr', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
subprocess.run(['lpsconstelm', '-v', '-t', '-', 'temp.lps'], input=run.stdout, check=True)


run = subprocess.run(['lps2pbes', '-fyellow_has_a_winning_strategy.mcf', 'temp.lps'], stdout=subprocess.PIPE, check=True)
run = subprocess.run(['pbesconstelm', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['pbesparelm', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
subprocess.run(['pbesrewr', '-pquantifier-all', '-', 'temp_yellow.pbes'], input=run.stdout, check=True)

run = subprocess.run(['lps2pbes', '-fred_has_a_winning_strategy.mcf', 'temp.lps'], stdout=subprocess.PIPE, check=True)
run = subprocess.run(['pbesconstelm', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['pbesparelm', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
subprocess.run(['pbesrewr', '-pquantifier-all', '-', 'temp_red.pbes'], input=run.stdout, check=True)

if "-rjittyc" in argv:
    subprocess.run(['lps2lts', '-v', '-rjittyc', 'temp.lps', 'tictactoe.aut'], check=True)
    subprocess.run(['pbessolve', '-v', '-s2', '-rjittyc', '-esome', 'temp_yellow.pbes'], check=True)
    subprocess.run(['pbessolve', '-v', '-s2', '-rjittyc', '-esome', 'temp_red.pbes'], check=True)