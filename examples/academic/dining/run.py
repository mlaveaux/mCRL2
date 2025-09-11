#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

for name in [
    'dining3',
    'dining3_cs',
    'dining3_cs_seq',
    'dining3_ns',
    'dining3_ns_seq',
    'dining3_schedule',
    'dining3_schedule_seq'
]:
    subprocess.run(['mcrl22lps', f'{name}.mcrl2', f'{name}.lps'], check=True)

    print(f'Checking property nodeadlock.mcf for specification {name}.mcrl2')
    subprocess.run(['lps2pbes', '-f', 'nodeadlock.mcf', f'{name}.lps', f'{name}.nodeadlock.pbes'], check=True)
    subprocess.run(['pbessolve', f'{name}.nodeadlock.pbes'], check=True)
  
    print(f'Checking property nostuffing.mcf for specification {name}.mcrl2')
    subprocess.run(['lps2pbes', '-f', 'nostuffing.mcf', f'{name}.lps', f'{name}.nostuffing.pbes'], check=True)
    subprocess.run(['pbessolve', f'{name}.nostuffing.pbes'], check=True)

    print(f'Checking property nostarvation.mcf for specification {name}.mcrl2')
    subprocess.run(['lps2pbes', '-f', 'nostarvation.mcf', f'{name}.lps', f'{name}.nostarvation.pbes'], check=True)
    subprocess.run(['pbessolve', f'{name}.nostarvation.pbes'], check=True)
