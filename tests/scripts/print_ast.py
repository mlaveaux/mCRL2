#!/usr/bin/env python3

# Copyright 2024-2025 Maurice Laveaux
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

import argparse
from enum import Enum
import os
from pathlib import Path
import sys
import shutil
import concurrent.futures

# Makes sure that the script can find the modules when ran directly.
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))
from tests.utility.run_process import RunProcess, TimeExceededError

class Variant(Enum):
    MCRL2 = 1
    MCF = 2

def main():
    """Execute all the run.py scripts in the example directory."""

    # Parse some configuration options
    parser = argparse.ArgumentParser(
        prog="print_ast.py",
        description="Prints all the ASTs of the the example specifications",
        epilog="",
    )

    parser.add_argument("-m", "--max_workers", action="store", default=1, type=int)
    parser.add_argument("-t", "--toolpath", action="store", type=str, required=True)

    args = parser.parse_args()

    os.environ["PATH"] = args.toolpath + os.pathsep + os.environ["PATH"]

    mcrl22lps_path = shutil.which("mcrl22lps")
    lps2pbes_path = shutil.which("lps2pbes")
    
    def run_mcrl2_file(mcrl22lps_path, path, index):
        print(f"[{index}] Printing mcrl2 {path}", flush=True)

        process = RunProcess(
            mcrl22lps_path,
            ["--print-ast", path],
            2000,
            1000,
        )

        # Write the AST output to a file
        output_path = Path(path).with_suffix('.txt')
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(process.stdout)

        return process.stdout

    def run_mcf_file(lps2pbes, path, index):
        print(f"[{index}] Printing mcf {path}", flush=True)

        process = RunProcess(
            lps2pbes,
            ["--print-ast", f"--formula={path}", "abp.lps", "-"],
            2000,
            1000,
        )

        # Write the AST output to a file
        output_path = Path(path).with_suffix('.txt')
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(process.stdout)

        return process.stdout

    # Start the linearisation process
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=args.max_workers
    ) as executor:
        # Collect all .mcrl2 files
        run_scripts = []
        for root, _, files in os.walk(
            os.path.join(os.path.realpath(os.path.dirname(__file__)), "../../examples/")
        ):
            for file in files:
                if file.endswith(".mcrl2"):
                    run_scripts.append((Variant.MCRL2, os.path.join(root, file)))
                elif file.endswith(".mcf"):
                    run_scripts.append((Variant.MCF, os.path.join(root, file)))

        # Submit the jobs
        futures = {}
        for index, (variant, path) in enumerate(run_scripts):
            if variant == Variant.MCRL2:
                futures[(executor.submit(run_mcrl2_file, mcrl22lps_path, path, index))] = path
            elif variant == Variant.MCF:
                futures[(executor.submit(run_mcf_file, lps2pbes_path, path, index))] = path

        # Wait for jobs to finish and show a progress bar.
        finished = 0
        failed = []
        for future in concurrent.futures.as_completed(futures):
            path = futures[future]

            if future.exception():
                if isinstance(future.exception(), TimeExceededError):
                    # Ignore timeouts, they are not errors.
                    print(
                        f"[{finished}/{len(run_scripts)}] Example {path} timed out",
                        flush=True,
                    )
                else:
                    print(
                        f"[{finished}/{len(run_scripts)}] Example {path} failed with {future.exception()}",
                        flush=True,
                    )
                    failed.append(path)
            elif future.cancelled():
                print(
                    f"[{finished}/{len(run_scripts)}] Example {path} was cancelled",
                    flush=True,
                )
            else:
                print(
                    f"[{finished}/{len(run_scripts)}] Finished example {path}",
                    flush=True,
                )

            finished += 1

        print("Failed tests...")
        for name in failed:
            print(f"{name} failed")

        if failed:
            sys.exit(-1)


if __name__ == "__main__":
    main()
