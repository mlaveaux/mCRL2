#!/usr/bin/env python3

import subprocess
import os

from sys import argv

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

run = subprocess.run(['mcrl22lps', '-vnb', 'ertms-hl3.mcrl2'], stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpssumelm', '-vc'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsconstelm', '-v'], input=run.stdout, stdout=subprocess.PIPE, check=True)
run = subprocess.run(['lpsparelm', '-v', '-', 'ertms-hl3.lps'], input=run.stdout, stdout=subprocess.PIPE, check=True)

if '-rjittyc' in argv:
    subprocess.run(['lps2lts', '-v', '-rjittyc', '--cached', '--timings=lps2lts_times.txt', 'ertms-hl3.lps', 'ertms-hl3.mod.aut'], check=True)

    subprocess.run(['ltsconvert', '--timings=ltsconvert_times.txt', '--tau=break,split_train,enter,leave,extend_EoA,move,connect,disconnect', '-v', '-edpbranching-bisim', 'ertms-hl3.mod.aut', 'ertms-hl3.mod.min.aut'], check=True)
    subprocess.run(['ltsconvert', '-vl', 'ertms-hl3.lps', 'ertms-hl3.mod.min.aut', 'ertms-hl3.mod.min.lts'], check=True)
    
    print('Verifying strong determinacy')
    subprocess.run(['lts2pbes', '-vf', 'strong_determinacy.mcf', 'ertms-hl3.mod.min.lts', '-l', 'ertms-hl3.lps', 'ertms-hl3.mod.min.strong_determinacy.pbes'], check=True)
    subprocess.run(['pbessolve', '-v', '-s2', '--timings=pbessolve_strong_determinacy_times.txt', 'ertms-hl3.mod.min.strong_determinacy.pbes'], check=True)

    print('Verifying termination')
    subprocess.run(['lts2pbes', '-vf', 'termination.mcf', 'ertms-hl3.mod.min.lts', '-l', 'ertms-hl3.lps', 'ertms-hl3.mod.min.termination.pbes'], check=True)
    subprocess.run(['pbessolve', '-v', '-s2', '--timings=pbessolve_termination_times.txt', 'ertms-hl3.mod.min.termination.pbes'], check=True)

    print('Verifying deterministic stabilisation')
    subprocess.run(['lts2pbes', '-vf', 'deterministic_stabilisation.mcf', 'ertms-hl3.mod.min.lts', '-l', 'ertms-hl3.lps', 'ertms-hl3.mod.min.deterministic_stabilisation.pbes'], check=True)
    subprocess.run(['pbessolve', '-v', '-s2', '--timings=pbessolve_times.txt', 'ertms-hl3.mod.min.deterministic_stabilisation.pbes'], check=True)

    # This requires the SU to be explored first.
    #print('checking whether IU included in SU')
    #subprocess.run(['ltscompare', '-vpweak-trace-ac', '--tau=change,continue,ptd_continue,ptd_stable', 'ertms-hl3.mod.aut', '../SU/ertms-hl3.mod.aut'], check=True)
    #print('checking whether SU included in IU')
    #subprocess.run(['ltscompare', '-vpweak-trace-ac', '--tau=change,continue,ptd_continue,ptd_stable', '../SU/ertms-hl3.mod.aut', 'ertms-hl3.mod.aut'], check=True)

