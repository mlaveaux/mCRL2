#!/bin/bash

benchmark() {
    for n in 1 2 3 4 5:
    do
        time $1 $2 1
    done
}

(
    set -x
    
    benchmark benchmark_target_atermpp_shared_creation 2 &&
    benchmark benchmark_target_atermpp_shared_lookup 2 &&
    benchmark benchmark_target_atermpp_shared_inspect 2 &&

    benchmark benchmark_target_atermpp_unique_creation 2 &&
    benchmark benchmark_target_atermpp_unique_lookup 2 &&
    benchmark benchmark_target_atermpp_unique_inspect 2
) 2>&1 | tee $1 &&

./report.py $1
