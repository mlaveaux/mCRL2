#!/usr/bin/env python

# Copyright 2015-2025 Wieger Wesselink, Maurice Laveaux
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

import argparse
import os
import sys
import concurrent.futures
from typeguard import typechecked

from typing import List, Dict, Any, Optional

# Makes sure that the script can find the modules when ran directly.
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))

from tests.utility.run_process import RunProcess

MCRL2_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


class LpsFile:
    pass


class LtsFile:
    pass


class ToolTest:
    @typechecked
    def __init__(self, toolname: str, file_input: str, options: List[str]):
        self.toolname: str = toolname
        self.file_input: str = file_input
        self.options: List[str] = options


available_tests: List[ToolTest] = [
    ToolTest(
        "mcrl22lps",
        "[mcrl2file]",
        [
            "-a",
            "-b",
            "-c",
            "-e",
            "-f",
            "-g",
            "-lregular",
            "-lregular2",
            "-lstack",
            "-m",
            "-n",
            "--no-constelm",
            "-o",
            "-rjitty",
            "-rjittyp",
            "--timings",
            "-w",
            "-z",
        ],
    ),
    ToolTest("lpsinfo", "[lpsfile]", ["-v"]),
    ToolTest("lpspp", "[lpsfile]", ["-fdefault", "-finternal"]),
    ToolTest(
        "lps2lts",
        "[lpsfile]",
        [
            "-ctau",
            "-D",
            "--error-trace",
            "-rjitty",
            "-rjittyp",
            "-rjittyc",
            "-sd",
            "-sb",
            "-sh --todo-max=100",
            "--verbose --suppress",
            "-u",
        ],
    ),
    ToolTest(
        "lpsconstelm",
        "[lpsfile]",
        ["-c", "-f", "-s", "-t", "-rjitty", "-rjittyp", "-rjittyc"],
    ),
    ToolTest("lpsrewr", "[lpsfile]", ["-rjitty", "-rjittyp", "-rjittyc"]),
    ToolTest("lpsparelm", "[lpsfile]", ["-v"]),
    ToolTest("lpssumelm", "[lpsfile]", ["-c"]),
    ToolTest(
        "lpsactionrename",
        "[lpsfile]",
        list(
            map(
                lambda x: "-f[rename] " + x,
                ["-m", "-o", "-t", "-rjitty", "-rjittyp", "-rjittyc"],
            )
        ),
    ),
    ToolTest("lpsuntime", "[lpsfile]", ["-v", "-i", "-e"]),
    ToolTest(
        "lpsinvelm",
        "[lpsfile]",
        list(
            map(
                lambda x: "-i[trueinv] " + x,
                [
                    "-c",
                    "-e",
                    "-l",
                    "-n",
                    "-pinvelm_test",
                    "-rjitty",
                    "-t10",
                    "-rjittyc",
                    "-zcvc",
                ],
            )
        ),
    ),
    ToolTest(
        "lpsconfcheck",
        "[lpsfile]",
        [
            "-a",
            "-c",
            "-g",
            "-i[trueinv]",
            "-i[falseinv]",
            "-p[dotfile] -i[trueinv]",
            "-m",
            "-n",
            "-o",
            "-rjitty",
            "-rjittyc",
            "-s10 -i[trueinv]",
            "-t10 -i[trueinv]",
        ],
    ),
    ToolTest(
        "lpsbinary",
        "[lpsfile]",
        ["-v", "-rjitty", "-rjittyp", "-rjittyc"],
    ),
    ToolTest(
        "lpsparunfold",
        "[lpsfile]",
        list(
            map(
                lambda x: "-sNat " + x,
                ["-a", "-p", "-n10", "-rjitty", "-rjittyp", "-rjittyc"],
            )
        )
        + ["-i0"],
    ),
    ToolTest(
        "lpssuminst",
        "[lpsfile]",
        list(map(lambda x: "--finite " + x, ["-rjitty", "-rjittyp", "-rjittyc"])),
    ),
    ToolTest(
        "ltsinfo",
        "[ltsfile]",
        ["-v"],
    ),
    ToolTest(
        "lts2lps",
        "[ltsfile]",
        ["-m [mcrl2file]"],
    ),
    ToolTest(
        "ltsconvert",
        "[ltsfile]",
        list(
            map(
                lambda x: "-oaut " + x,
                [
                    "-D",
                    "-ebisim",
                    "-ebranching-bisim",
                    "-edpbranching-bisim",
                    "-esim",
                    "-etau-star",
                    "-etrace",
                    "-eweak-trace",
                    "-n",
                    "--no-reach",
                    "--tau=[actions]",
                ],
            )
        ),
    ),
    ToolTest(
        "ltscompare",
        "[ltsfile] [ltsfile]",
        [
            "-c -ebranching-bisim",
            "-c -edpbranching-bisim",
            "-c -esim",
            "-c -etrace",
            "-c -eweak-trace",
            "-psim",
            "-pweak-trace",
            "-ebisim --tau=[actions]",
        ],
    ),
    ToolTest("lps2pbes", "[lpsfile]", ["-f[mcffile]", "-t -f[mcffile]"]),
    ToolTest("lts2pbes", "[ltsfile]", ["-f[mcffile]"]),
    ToolTest("pbessolve", "[pbesfile]", ["-v"]),
    ToolTest("pbespp", "[pbesfile]", ["-finternal", "-fdefault"]),
    ToolTest(
        "pbesconstelm",
        "[pbesfile]",
        [
            "-c",
            "-e",
            "-psimplify",
            "-pquantifier-all",
            "-pquantifier-finite",
            "-ppfnf",
            "-rjitty",
            "-rjittyp",
            "-rjittyc",
        ],
    ),
    ToolTest("pbesparelm", "[pbesfile]", ["-v"]),
    ToolTest(
        "pbesrewr",
        "[pbesfile]",
        [
            "-psimplify",
            "-pquantifier-all",
            "-pquantifier-finite",
            "-ppfnf",
            "-rjitty",
            "-rjittyp",
            "-rjittyc",
        ],
    ),
    ToolTest("txt2lps", "[txtlpsfile]", ["-v"]),
    ToolTest(
        "pbes2bes",
        "[pbesfile]",
        [
            "-opbes",
            "-otext",
            "-obes",
            "-opgsolver",
            "-rjitty",
            "-rjittyp",
            "-rjittyc",
            "-s0",
            "-s1",
            "-s2",
            "-s3",
            "-u",
            "-zb",
            "-zd",
            "--erase=none",
            "--erase=some",
            "--erase=all",
        ],
    ),
    ToolTest("pbesinfo", "[pbesfile]", ["-v", "-f"]),
    ToolTest("besinfo", "[pbesfile]", ["-v", "-f"]),
    ToolTest("bespp", "[pbesfile]", ["-fdefault"]),
    # TODO: "-sspm" causes a segfault.
    ToolTest("bessolve", "[pbesfile]", ["-sgauss"]),
    ToolTest("besconvert", "[besfile]", ["-ebisim", "-estuttering"]),
    ToolTest(
        "pbesinst",
        "[pbesfile]",
        [
            "-v",
            "-opbes",
            "-obes",
            "-otext",
            "-opgsolver",
            "-slazy",
            '-sfinite -f"*(*:Bool)"',
            "-rjitty",
            "-rjittyp",
            "-rjittyc",
        ],
    ),
    ToolTest(
        "pbespareqelm",
        "[pbesfile]",
        [
            "-v",
            "-I",
            "-psimplify",
            "-pquantifier-all",
            "-pquantifier-finite",
            "-ppfnf",
            "-rjitty",
            "-rjittyp",
            "-rjittyc",
        ],
    ),
    ToolTest(
        "pbespgsolve",
        "[pbesfile]",
        [
            "-v",
            "-c",
            "-C",
            "-L",
            "-e",
            "-sspm",
            "-saltspm",
            "-srecursive",
            "-rjitty",
            "-rjittyp",
            "-rjittyc",
        ],
    ),
    # TODO: This test takes too much time
    # ToolTest(
    #     "lpsbisim2pbes",
    #     "[lpsfile]",
    #     [
    #         "-bstrong-bisim",
    #         "-bweak-bisim",
    #         "-bbranching-bisim",
    #         "-bbranching-sim",
    #         "-bstrong-bisim -n",
    #     ],
    # ),
    ToolTest("txt2bes", "[txtpbesfile]", ["-v"]),
    ToolTest(
        "lpsrealelm", "[lpsfile]", ["-v", "--max=10", "-rjitty", "-rjittyp", "-rjittyc"]
    ),
]


@typechecked
def run_test(toolpath: str, arguments: List[str]) -> Any:
    return RunProcess(toolpath, arguments, 1000, 600)

@typechecked
def main(tests: List[ToolTest]) -> Optional[int]:
    # Parse the input arguments.
    cmdline_parser = argparse.ArgumentParser()
    cmdline_parser.add_argument(
        "-t",
        "--toolpath",
        dest="toolpath",
        help="The path where the mCRL2 tools are installed",
    )
    cmdline_parser.add_argument(
        "-o",
        "--output",
        dest="output",
        metavar="o",
        action="store",
        help="Run the tests in the given directory, where the outputs are produced",
    )
    cmdline_parser.add_argument(
        "-m", "--max_workers", action="store", default=1, type=int
    )
    cmdline_parser.add_argument("--names", action="store_true")
    cmdline_parser.add_argument("--verbose", default=False, action="store_true")
    cmdline_parser.add_argument(
        "--experimental",
        action="store_true",
        help="Enable testing of experimental tests",
    )
    cmdline_parser.add_argument("--pattern", action="store", default="", type=str)
    cmdline_parser.add_argument(
        "--jittyc", action="store_true", help="Enables tests using the -rjittyc flag"
    )
    cmdline_parser.add_argument(
        "-cvc3",
        action="store_true",
        help="Enables tests using the -zcvc flag, requires CVC3",
    )

    args: argparse.Namespace = cmdline_parser.parse_args()

    output: str
    if args.output:
        output = os.path.join(os.getcwd(), args.output)
    else:
        output = "."

    if args.names:
        # Generate a name per line for a test
        names: str
        names = ""
        for test in tests:
            names += test.toolname + "\n"
        print(names)
        return

    for example, rename_spec, mcrl2_file in [
        (
            "abp",
            "act renamed;\n var x:D;\n rename r1(x) => renamed;",
            os.path.join(MCRL2_ROOT, "examples/academic/abp/abp.mcrl2"),
        ),
        # ("dice", "act renamed;\n var x:Bool;\n rename b1(x) => renamed;", os.path.join(MCRL2_ROOT, "examples/probabilistic/coins_simulate_dice/dice.mcrl2"))
    ]:
        print(f"Write intermediate files into {output}", flush=True)
        name: str
        name, _ = os.path.splitext(os.path.basename(mcrl2_file))
        lps_file: str = os.path.join(output, f"{name}.lps")
        lts_file: str = os.path.join(output, f"{name}.lts")

        # Generate textual versions from the lps and pbes
        txtlps_file: str = os.path.join(output, f"{name}.txtlps")
        txtpbes_file: str = os.path.join(output, f"{name}.txtpbes")

        # Create input files for lpsconfcheck and lpsinvelm
        true_file: str = os.path.join(output, "true.inv")
        false_file: str = os.path.join(output, "false.inv")

        with open(true_file, "w", encoding="utf-8") as f:
            f.write("true")
        with open(false_file, "w", encoding="utf-8") as f:
            f.write("false")

        # Create input files for lpsactionrename
        rename_file: str = os.path.join(output, f"{example}.rename")
        with open(rename_file, "w", encoding="utf-8") as f:
            f.write(rename_spec)

        # Create a generic nodeadlock formula
        mcf_file: str = os.path.join(output, "nodeadlock.mcf")
        with open(mcf_file, "w", encoding="utf-8") as f:
            f.write("[true*]<true>true")
        pbes_file: str = os.path.join(output, f"{example}.nodeadlock.pbes")
        bes_file: str = os.path.join(output, f"{example}.nodeadlock.bes")

        # Generate the input files for all tests beforehand
        run_test(os.path.join(args.toolpath, "mcrl22lps"), ["-D", mcrl2_file, lps_file])
        run_test(os.path.join(args.toolpath, "lps2lts"), [lps_file, lts_file])
        run_test(
            os.path.join(args.toolpath, "lps2pbes"),
            [lps_file, "-f", mcf_file, pbes_file],
        )
        run_test(os.path.join(args.toolpath, "lpspp"), [lps_file, txtlps_file])

        run_test(os.path.join(args.toolpath, "pbespp"), [pbes_file, txtpbes_file])
        run_test(os.path.join(args.toolpath, "pbes2bes"), [pbes_file, bes_file])

        try:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=args.max_workers
            ) as executor:
                # Start the test cases
                futures: Dict[Any, str] = {}
                for test in tests:
                    for argument in test.options:
                        toolpath: str = os.path.join(args.toolpath, test.toolname)

                        # Replace placeholders by actual filenames
                        final_arguments: List[str] = []
                        for arg in argument.split(" ") + test.file_input.split(" "):
                            arg = arg.replace("[mcrl2file]", mcrl2_file)
                            arg = arg.replace("[lpsfile]", lps_file)
                            arg = arg.replace("[ltsfile]", lts_file)
                            arg = arg.replace("[trueinv]", true_file)
                            arg = arg.replace("[falseinv]", false_file)
                            arg = arg.replace("[rename]", rename_file)
                            arg = arg.replace("[mcffile]", mcf_file)
                            arg = arg.replace("[pbesfile]", pbes_file)
                            arg = arg.replace("[besfile]", bes_file)
                            arg = arg.replace("[txtlpsfile]", txtlps_file)
                            arg = arg.replace("[txtpbesfile]", txtpbes_file)
                            arg = arg.replace("[actions]", "r1")
                            final_arguments.append(arg)

                        if not args.jittyc and "-rjittyc" in argument:
                            print(f"Skipped jittyc test {test.toolname} {' '.join(final_arguments)}")
                            continue
                        if not args.cvc3 and "-zcvc" in argument:
                            print(f"Skipped cvc test {test.toolname} {' '.join(final_arguments)}")
                            continue
                        if args.pattern not in test.toolname:
                            print(
                                f"Skipped tool {test.toolname} with does not fit pattern {args.pattern}"
                            )
                            continue
                        if not args.experimental and test.toolname in [
                            "besconvert",
                            "pbespareqelm",
                            "lpsrealelm",
                        ]:
                            print(f"Skipped experimental tool {test.toolname}")
                            continue

                        if args.verbose:
                            print(f"Executing {test.toolname} {final_arguments}", flush=True)
                        futures[
                            executor.submit(run_test, toolpath, final_arguments)
                        ] = f"{test.toolname} {argument}"

                # Wait for the results and keep track of failed tests.
                failed: List[str] = []
                for future in concurrent.futures.as_completed(futures):
                    name: str = futures[future]

                    try:
                        future.result()
                        if args.verbose:
                            print(f"{name} {(100 - len(name)) * '.'} ok", flush=True)
                    except Exception as e:
                        print(
                            f"{name} {(100 - len(name)) * '.'} failed \n  {e}",
                            flush=True,
                        )
                        failed.append(name)

                print("Failed tests...")
                for name in failed:
                    print(f"{name} failed")

                if failed:
                    return os._exit(1)
        except KeyboardInterrupt:
            print("Killed")


if __name__ == "__main__":
    main(available_tests)
