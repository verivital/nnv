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


