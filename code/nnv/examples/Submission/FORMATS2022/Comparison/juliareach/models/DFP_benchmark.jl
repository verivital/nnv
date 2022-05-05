using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

# How to keep track of time
# eltime3 = @elapsed sol3 = solve(prob, T=Tf, alg=alg, params=param3);

#SUITE = BenchmarkGroup()
#model = "DFP"
#cases = ["1"]
#SUITE[model] = BenchmarkGroup()

include("CTRNN_DFP.jl")
include("convert_sets.jl")

# ------------------------
# Case 1 
# ------------------------

#prob = CTRNN_DFP()
#alg = TMJets(abstol=1e-10, orderT=8, orderQ=1);

#time_1 = @elapsed sol_1 = solve(prob, T=10.0, alg=alg);
#solz_1 = overapproximate(sol_1, Zonotope);


# ------------------------
# Case 2 (long - this works better)
# ------------------------

prob = CTRNN_DFP()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=3500)
#alg = TMJets(abstol=1e-10, orderT=orderT, orderQ=1);

time_2 = @elapsed sol_2 = solve(prob, T=10.0, alg=alg);
solz_2 = overapproximate(sol_2, Zonotope);

# ------------------------
# Case 3 (mid)
# ------------------------

prob = CTRNN_DFP()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=1000)
#alg = TMJets(abstol=1e-10, orderT=orderT, orderQ=1);

time_mid = @elapsed sol_mid = solve(prob, T=2.5, alg=alg);
solz_mid = overapproximate(sol_mid, Zonotope);

# ------------------------
# Case 2 (short)
# ------------------------

prob = CTRNN_DFP()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=400)
#alg = TMJets(abstol=1e-10, orderT=orderT, orderQ=1);

time_short = @elapsed sol_short = solve(prob, T=0.5, alg=alg);
solz_short = overapproximate(sol_short, Zonotope);


# convert_sets(solz_1, "results/dfp1.mat")
convert_sets(solz_2, "results/dfp2.mat")
convert_sets(solz_mid, "results/dfp_mid.mat")
convert_sets(solz_short, "results/dfp_short.mat")

matwrite("results/dfp_time.mat", Dict("time_short" => time_short,"time_mid" => time_mid, "time_2" => time_2); compress = false)
