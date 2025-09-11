#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', 'mpsu.mcrl2', 'mpsu.lps'], check=True)

for property in [
    'nodeadlock.mcf',
    'mpsu1.mcf',
    'mpsu2.mcf',
    'mpsu3.mcf',
    'mpsu4.mcf',
    'mpsu5.mcf',
    'mpsu6.mcf'
]:
    print(f'Checking property {property} for "mpsu.lps')
    name, _ = os.path.splitext(property)

    subprocess.run(['lps2pbes', '-f', property, 'mpsu.lps', f'mpsu.{name}.pbes'], check=True)
    subprocess.run(['pbessolve', f'mpsu.{name}.pbes'], check=True)