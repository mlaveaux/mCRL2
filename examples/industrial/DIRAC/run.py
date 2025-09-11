#!/usr/bin/env python3

import subprocess
import os

from shutil import which

# Change working dir to the script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run(['mcrl22lps', 'SMS.mcrl2', 'SMS.lps'], check=True)
subprocess.run(['lps2lts', '-v', 'SMS.lps'], check=True)

subprocess.run(['lps2pbes', '-f', 'properties_SMS/eventuallyDeleted.mcf', 'SMS.lps', 'eventuallyDeleted.pbes'], check=True)
subprocess.run(['lps2pbes', '-f', 'properties_SMS/noTransitFromDeleted.mcf', 'SMS.lps', 'noTransitFromDeleted.pbes'], check=True)
subprocess.run(['pbessolve', 'eventuallyDeleted.pbes'], check=True)
subprocess.run(['pbessolve', 'noTransitFromDeleted.pbes'], check=True)

subprocess.run(['mcrl22lps', '-n', 'WMS.mcrl2', 'WMS.lps'], check=True)

lpsreach = which('lpsreach')
if lpsreach is not None:
    subprocess.run([lpsreach, '-v', '--chaining', '--groups=simple', '--saturation', 'WMS.lps'], check=True)


subprocess.run(['lps2pbes', '-f', 'properties_WMS/jobFailedToDone.mcf', 'WMS.lps', 'jobFailedToDone.pbes'], check=True)
subprocess.run(['lps2pbes', '-f', 'properties_WMS/noZombieJobs.mcf', 'WMS.lps', 'noZombieJobs.pbes'], check=True)

pbessolvesymbolic = which('pbessolvesymbolic')
if pbessolvesymbolic is not None:
    subprocess.run([pbessolvesymbolic, '-v', '--chaining', '--groups=simple', '--saturation', 'jobFailedToDone.pbes'], check=True)
    subprocess.run([pbessolvesymbolic, '-v', '--chaining', '--groups=simple', '--saturation', 'noZombieJobs.pbes'], check=True)

