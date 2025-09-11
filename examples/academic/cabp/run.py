#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', 'cabp.mcrl2', 'cabp.lps'], check=True)

for prop in [
    'nodeadlock.mcf',
    'infinitely_often_enabled_then_infinitely_often_taken.mcf',
    'infinitely_often_lost.mcf',
    'infinitely_often_receive_for_all_d.mcf',
    'read_then_eventually_send.mcf',
    'infinitely_often_receive_d1.mcf',
    'read_then_eventually_send_if_fair.mcf',
    'no_generation_of_messages.mcf',
    'no_duplication_of_messages.mcf'
]:
    path, _ = os.path.splitext(prop)
    name = os.path.basename(path)

    print(f'verifying property {name}.mcf for cabp.lps')
    subprocess.run(['lps2pbes', '-f', prop, 'cabp.lps', f'cabp.{name}.pbes'], check=True)
    subprocess.run(['pbessolve', f'cabp.{name}.pbes'], check=True)