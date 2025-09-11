#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def run(mcffile, lpsfile):
    print(f'verifying property {os.path.basename(mcffile)}')

    run = subprocess.run(['lps2pbes', '-f', mcffile, lpsfile], stdout=subprocess.PIPE, check=True)
    subprocess.run(['pbessolve'], input=run.stdout, check=True)

lpsfile = '11073.lps'
subprocess.run(['mcrl22lps', '11073.mcrl2', lpsfile], check=True)

run('data_can_be_communicated.mcf', lpsfile)
run('infinite_data_communication_is_possible.mcf', lpsfile)
run('nodeadlock.mcf', lpsfile)
run('no_inconsistent_operating_states.mcf', lpsfile)
run('no_successful_transmission_in_inconsistent_operating_states.mcf', lpsfile)
