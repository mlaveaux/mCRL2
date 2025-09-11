#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-vo', 'othello.mcrl2', 'othello.lps'], check=True)
run = subprocess.run(['lpssuminst', '-v', '-rjitty', 'othello.lps'], stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsparunfold', '-n6', '-v', '-sBoard'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsrewr', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsconstelm', '-c', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsrewr', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsparunfold', '-n6', '-v', '-sRow'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsrewr', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsconstelm', '-c', '-v'], input=run.stdout, stdout=subprocess.PIPE,check=True)
subprocess.run(['lpsrewr', '-v', '-', 'temp.lps'], input=run.stdout, check=True)

if '-rjittyc' in argv:
    subprocess.run(['lps2lts', '-v', '-rjittyc', 'temp.lps', 'othello.lts'], check=True)

print('---------------------------------------------- ')
print(' Red can always win, i.e. red has a winning strategy.')
print('---------------------------------------------- ')
subprocess.run(['lps2pbes', '-v', '-fred_wins_always.mcf', 'temp.lps', 'othello.pbes'], check=True)
if '-rjittyc' in argv:
    subprocess.run(['pbessolve', '-v', '-rjittyc', 'othello.pbes'], check=True)

print('----------------------------------------------')
print(' White can always win, i.e. white has a winning strategy.')
print('----------------------------------------------')
subprocess.run(['lps2pbes', '-v', '-fwhite_wins_always.mcf', 'temp.lps', 'othello.pbes'], check=True)
if '-rjittyc' in argv:
    subprocess.run(['pbessolve', '-v', '-rjittyc', 'othello.pbes'], check=True)

print('----------------------------------------------')
print(' Nodeadlock.')
print('----------------------------------------------')
subprocess.run(['lps2pbes', '-v', '-fnodeadlock.mcf', 'temp.lps', 'othello.pbes'], check=True)
if '-rjittyc' in argv:
    subprocess.run(['pbessolve', '-v', '-rjittyc', 'othello.pbes'], check=True)

print('----------------------------------------------')
print(' There exists a draw.')
print('----------------------------------------------')
subprocess.run(['lps2pbes', '-v', '-fexists_draw.mcf', 'temp.lps', 'othello.pbes'], check=True)
if '-rjittyc' in argv:
    subprocess.run(['pbessolve', '-v', '-rjittyc', 'othello.pbes'], check=True)

print('----------------------------------------------')
print(' White can win.')
print('----------------------------------------------')
subprocess.run(['lps2pbes', '-v', '-fwhite_can_win.mcf', 'temp.lps', 'othello.pbes'], check=True)
if '-rjittyc' in argv:
    subprocess.run(['pbessolve', '-v', '-rjittyc', 'othello.pbes'], check=True)

print('----------------------------------------------')
print(' Red can win.')
print('----------------------------------------------')
subprocess.run(['lps2pbes', '-v', '-fred_can_win.mcf', 'temp.lps', 'othello.pbes'], check=True)
if '-rjittyc' in argv:
    subprocess.run(['pbessolve', '-v', '-rjittyc', 'othello.pbes'], check=True)

