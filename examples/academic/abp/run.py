#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'abp.mcrl2', 'abp.lps'], check=True)
subprocess.run(['lps2lts', '-v', 'abp.lps', 'abp.aut'], check=True)

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

    print(f'verifying property {name}.mcf for abp.lps')
    subprocess.run(['lps2pbes', '-f', prop, 'abp.lps', f'abp.{name}.pbes'], check=True)
    subprocess.run(['pbessolve', f'abp.{name}.pbes'], check=True)