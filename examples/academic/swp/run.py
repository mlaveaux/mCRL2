#!/usr/bin/env python3

import subprocess
import shutil
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'swp_lists.mcrl2', 'swp_lists.lps'], check=True)

for prop in [
    'nodeadlock.mcf',
    'infinitely_often_enabled_then_infinitely_often_taken.mcf', 
    'infinitely_often_lost.mcf',
    'infinitely_often_receive_d1.mcf',
    'infinitely_often_receive_for_all_d.mcf',
    'read_then_eventually_send.mcf',
    'read_then_eventually_send_if_fair.mcf',
    'no_generation_of_messages.mcf',
    'no_duplication_of_messages.mcf'
]:
    path, _ = os.path.splitext(prop)
    name = os.path.basename(path)

    print('verifying property {base} for swp_lists.lps')
    subprocess.run(['lps2pbes', '-v', '-f', prop, 'swp_lists.lps', f'swp_lists.{name}.pbes'], check=True)
    subprocess.run(['pbessolve', '-v', f'swp_lists.{name}.pbes'], check=True)

# SWP with Tanenbaum's bug
subprocess.run(['mcrl22lps', '-v', 'swp_with_tanenbaums_bug.mcrl2', 'swp_with_tanenbaums_bug.lps'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'swp_with_tanenbaums_bug.lps', 'swp_with_tanenbaums_bug.nodeadlock.pbes'], check=True)

lpsreach = shutil.which('lpsreach')
if lpsreach and '--long' in argv:
    subprocess.run([lpsreach, '-v', '--cached', '--chaining', '--saturation', 'swp_with_tanenbaums_bug.lps'], check=True)

pbessolvesymbolic = shutil.which('pbessolvesymbolic')
if pbessolvesymbolic and '--long' in argv:
    subprocess.run([pbessolvesymbolic, '-v', '-c', '--cached', '--chaining', '--saturation', 'swp_with_tanenbaums_bug.nodeadlock.pbes'], check=True)

# It is possible to generate the state space, while looking for occurrence of
# error action, because the error action can be found without exploring the
# whole state space. The command below stores trace, and stops generating after
# finding a single trace. Removal of --cached leads to less memory requirements,
# at the expense of a longer running time. 
subprocess.run(['lps2lts', '-v', '--cached', '-aerror', '-t1', 'swp_with_tanenbaums_bug.lps'], check=True)
