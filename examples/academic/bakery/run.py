#!/usr/bin/env', 'python3

import subprocess
import os

from shutil import which

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

mcrl22lps = subprocess.run(['mcrl22lps', '-f', '-D', '-n', '-w', 'bakery.mcrl2'], stdout=subprocess.PIPE, check=True)
lpsconstelm = subprocess.run(['lpsconstelm'], input=mcrl22lps.stdout, stdout=subprocess.PIPE, check=True)
subprocess.run(['lpsparelm', '-', 'bakery.lps'], input=lpsconstelm.stdout, check=True)

# Two different abstractions used (an aggressive and less aggressive variation)
# the following four properties can be proved/disproved using pbesabsinthe:
#  1. get_at_least_number_circulating
#  2. always_can_get_number
#  3. request_must_eventually_enter
#  4. request_can_eventually_enter

# the remaining two properties cannot (currently) be proved/disproved using pbesabsinthe.
for formula in [
    'nodeadlock',
    'request_can_eventually_enter',
    'request_must_eventually_enter',
    'mutual_exclusion',
    'always_can_get_number',
    'get_at_least_number_circulating',
]:
    print(f'Solving property {formula}')
    run = subprocess.run(['lps2pbes', '-f', f'{formula}.mcf', 'bakery.lps'], stdout=subprocess.PIPE, check=True)
    run = subprocess.run(['pbesrewr', '-psimplify'], input=run.stdout, stdout=subprocess.PIPE, check=True)
    run = subprocess.run(['pbesconstelm'], input=run.stdout, stdout=subprocess.PIPE, check=True)
    subprocess.run(['pbesparelm', '-', f'bakery.{formula}.pbes'], input=run.stdout, stdout=subprocess.PIPE, check=True)
    
    # Does not terminate since the statespace is infinite
    #subprocess.run(['pbessolve', '-v', '-s3', '-rjitty', f'bakery.{formula}.pbes'], check=True)

    # Instead we can abstract the data types to obtain over and underapproximations. These use the partial solvers of pbessolve to terminate early.
    pbesabsinthe = which('pbesabsinthe')
    if pbesabsinthe is not None:
        print('Answer by under-approximation: ')
        subprocess.run([pbesabsinthe, '-sunder', '-a', 'abstraction.txt', f'bakery.{formula}.pbes', f'bakery.{formula}.absinthe.pbes'], check=True)
        subprocess.run(['pbessolve', '-s3', '-rjitty', f'bakery.{formula}.absinthe.pbes'], check=True)

        print('Answer by over-approximation: ')
        subprocess.run([pbesabsinthe, '-sover', '-a', 'abstraction.txt', f'bakery.{formula}.pbes', f'bakery.{formula}.absinthe.pbes'], check=True)
        subprocess.run(['pbessolve', '-s3', '-rjitty', f'bakery.{formula}.absinthe.pbes'], check=True)

        # Alternatively, we can abstract more aggressive to ensure that the parity game is finite.
        print('Answer by aggressive under-approximation: ')
        subprocess.run([pbesabsinthe, '-sunder', '-a', 'aggressive_abstraction.txt', f'bakery.{formula}.pbes', f'bakery.{formula}.absinthe.pbes'], check=True)
        subprocess.run(['pbessolve', '-rjitty', f'bakery.{formula}.absinthe.pbes'], check=True)

        print('Answer by aggressive over-approximation: ')
        subprocess.run([pbesabsinthe, '-sover', '-a', 'aggressive_abstraction.txt', f'bakery.{formula}.pbes', f'bakery.{formula}.absinthe.pbes'], check=True)
        subprocess.run(['pbessolve', '-rjitty', f'bakery.{formula}.absinthe.pbes'], check=True)