#!/usr/bin/env python3

import subprocess
import os

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', '-v', 'allow.mcrl2', 'allow.lps'], check=True)

subprocess.run(['lps2pbes', '-v', '-f', 'nodeadlock.mcf', 'allow.lps', 'allow.nodeadlock.pbes'], check=True)
subprocess.run(['pbessolve', '-v', 'allow.nodeadlock.pbes'], check=True)
# The above shows that there is a deadlock in the specification, let's investigate.

# The following command creates the state space, and stores trace to the deadlock
subprocess.run(['lps2lts', '-v', '-Dt', 'allow.lps', 'allow.aut'], check=True)
print('Print the trace and find out what\'s wrong: ')
subprocess.run(['tracepp', 'allow.lps_dlk_0.trc'], check=True)

# The same is also immediately clear if we visualise the state space, uncomment the following line.
# subprocess.run(['ltsgraph', 'allow.aut'], check=True)