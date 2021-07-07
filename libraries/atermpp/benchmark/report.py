#!/bin/python

import sys
import re

tool_regex = re.compile(r"\+ (benchmark.*) .* (.*)$")

# For the output of time
time_regex = re.compile(r"real	([0-9]*)m([0-9]*)\.([0-9]*)s$")
time_internal_regex = re.compile(r"time: (.*)$")

def main():
    if len(sys.argv) <= 1:
        print("Usage: report.py <filename> (where the filename contains the output of the run.sh scripts)")
        return -1

    name = ""
    k = 0

    table = {}

    # Open the given file and iterate over all lines.
    with open(sys.argv[1]) as file:

        for line in file:
            # The name of the benchmark
            info = tool_regex.search(line)
            if info is not None:
                name = info.group(1)
                k = info.group(2)
                #print(name + " " + k)

            info = time_regex.search(line)
            if info is not None:
                minutes = int(info.group(1))
                remaining = float("." + info.group(3)) 

                seconds = int(info.group(2)) + remaining

                #print("Required {}m{}s ({} seconds) time".format(minutes, round(seconds, 2), round(minutes * 60 + seconds, 2)))

                if not name in table:
                    table[name] = {}

                if not k in table[name]:
                    table[name][k] = {'count': 0, 'total': 0}

                #table[name][k]['total'] += minutes * 60 + seconds                
                #table[name][k]['count'] += 1

            info = time_internal_regex.search(line)
            if info is not None:
                seconds = float(info.group(1))

                if not name in table:
                    table[name] = {}

                if not k in table[name]:
                    table[name][k] = {'count': 0, 'total': 0}

                table[name][k]['total'] += seconds
                table[name][k]['count'] += 1

    print("")
    print("plot entries:\n")
    for name in table:
        print(name)
        for k in table[name]:
            print("({}, {})".format(k, round(table[name][k]['total'] / table[name][k]['count'], 2)))
        print("\n")


if __name__ == "__main__":
    main()
