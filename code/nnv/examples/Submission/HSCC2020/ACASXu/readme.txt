# These codes are generating the comparison between NNV tool, Reluplex, Marabou, Marabou_DNC, and ReluVal.

1. The data file "other_p12" contains results of testing the first 9 ACASXu networks with properties p1 and p2 on Reluplex, Marabou, Marabou_DNC, and ReluVal.
2. The data file "other_p34" contains results of testing the 45 ACASXu networks with properties p3 and p4 on Reluplex, Marabou, Marabou_DNC, and ReluVal.
3. run "verify_p3_p4" to test the NNV tool for properties p3 and p4, and generate a comparison table with the aforementioned approaches (recommended).   The NNV test results are stored into the variable "nnv". Running time: around 3 hours if we use 4 cores. 

4. run "verify_p1_p2" to test the NNV tool for properties p1 and p2, and generate a comparison table with the aforementioned approaches (not recommended due to large computational time and potential out-of-memory error). The NNV test results are stored into the variable "nnv". Running time: around 45 hours.
