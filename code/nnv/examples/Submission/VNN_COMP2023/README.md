# Local Testing for VNN-COMP Scripts

## File Summaries
`benchmark.sh`        -- a script for local testing of the benchmark testing sequence that will be launched for the competition.
`prepare_instance.sh` -- a script for starting a local matlab engine within which `nnv` can be ran on an instance.
`run_instance.sh`     -- a script for running `nnv` on an instance.

## `benchmark.sh`
For local testing only. Launches the `prepare_instance.sh` and `run_instance.sh` scripts for each row of a provided .csv file--which includes the parameters of the test instance.

## `prepare_instance.sh`
Provided for the competition. Handles error-handling and set up of MATLAB engine prior to running `nnv` on the test instance.

## `run_instance.sh`
Provided for the competition. Executes `nnv` on the test instance.


# Preparing for 2024 competition

TODO:

- Choose optimal reachability options for each benchmark.

- Change to the new ONNX import function (will this work from command line so no need to convert models prior to competition?)

- There were some possible errors with counterexamples based on last year report, take a look at this (some already fixed).

- Test the python integration with newest MATLAB version. Can we simplify the code and reduce computation time using the matlab engine for python? (This did not work properly on 2023a)


# Instructions to test locally (2023 benchmarks)

1) Ensure NNV is installed and the corresponding toolboxes

2) Go to https://github.com/ChristopherBrix/vnncomp2023_benchmarks, clone the repo and run (in your desired folder) ./setup.sh

3) Change line 4 in 'test_some_instances.m' to the folder where you downloaded the benchmarks.


3) 