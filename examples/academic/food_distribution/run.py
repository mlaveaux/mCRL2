#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

run = subprocess.run(['mcrl22lps', 'food_package.mcrl2'], stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpssuminst'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lps2pbes', '-f', 'sustained_delivery.mcf'], input=run.stdout, stdout=subprocess.PIPE, check=True)
subprocess.run(['pbesconstelm', '-e', '-', 'sustained_delivery.pbesconstelm.pbes'], input=run.stdout, stdout=subprocess.PIPE, check=True)

if '-rjittyc' in argv:
    # Note that the generated bes is huge.
    subprocess.run(['pbessolve', '-v', '-rjittyc', '-zdepth-first', '-s3', 'sustained_delivery.pbesconstelm.pbes'], check=True)
else:
    print('Example requires -rjittyc to run all steps')