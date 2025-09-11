#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'onebit.mcrl2', 'onebit.lps'], check=True)
subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'onebit.lps', 'onebit.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'onebit.nodeadlock.pbes'], check=True)

# Solving the following using pbessolve is slow as the formulas have nested fixed points!
subprocess.run(['lps2pbes', '-v', '-f', 'infinitely_often_enabled_then_infinitely_often_taken.mcf', 'onebit.lps', 'onebit.infinitely_often_enabled_then_infinitely_often_taken.pbes'], check=True)
subprocess.run(['pbespgsolve', '-srecursive', '-v', 'onebit.infinitely_often_enabled_then_infinitely_often_taken.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'infinitely_often_lost.mcf', 'onebit.lps', 'onebit.infinitely_often_lost.pbes'], check=True)
subprocess.run(['pbespgsolve', '-srecursive', '-v', 'onebit.infinitely_often_lost.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'infinitely_often_receive_dat1.mcf', 'onebit.lps', 'onebit.infinitely_often_receive_dat1.pbes'], check=True)
subprocess.run(['pbespgsolve', '-srecursive', '-v', 'onebit.infinitely_often_receive_dat1.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'infinitely_often_receive_for_all_d.mcf', 'onebit.lps', 'onebit.infinitely_often_receive_for_all_d.pbes'], check=True)
subprocess.run(['pbespgsolve', '-srecursive', '-v', 'onebit.infinitely_often_receive_for_all_d.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'read_then_eventually_send.mcf', 'onebit.lps', 'onebit.read_then_eventually_send.pbes'], check=True)
subprocess.run(['pbespgsolve', '-srecursive', '-v', 'onebit.read_then_eventually_send.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'read_then_eventually_send_if_fair.mcf', 'onebit.lps', 'onebit.read_then_eventually_send_if_fair.pbes'], check=True)
subprocess.run(['pbespgsolve', '-srecursive', '-v', 'onebit.read_then_eventually_send_if_fair.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'no_generation_of_messages.mcf', 'onebit.lps', 'onebit.no_generation_of_messages.pbes'], check=True)
subprocess.run(['pbespgsolve', '-srecursive', '-v', 'onebit.no_generation_of_messages.pbes'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'no_duplication_of_messages.mcf', 'onebit.lps', 'onebit.no_duplication_of_messages.pbes'], check=True)
subprocess.run(['pbespgsolve', '-srecursive', '-v', 'onebit.no_duplication_of_messages.pbes'], check=True)
