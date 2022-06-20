using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median


include("CTRNN_FPA.jl")
include("convert_sets.jl")


# ------------------------
# Case 2 (long - this works better)
# ------------------------

prob = CTRNN_FPA()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=3500)
#alg = TMJets(abstol=1e-10, orderT=orderT, orderQ=1);

time_2 = @elapsed sol_2 = solve(prob, T=10.0, alg=alg);
solz_2 = overapproximate(sol_2, Zonotope);

# ------------------------
# Case 3 (mid)
# ------------------------

prob = CTRNN_FPA()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=1000)
#alg = TMJets(abstol=1e-10, orderT=orderT, orderQ=1);

time_mid = @elapsed sol_mid = solve(prob, T=2.5, alg=alg);
solz_mid = overapproximate(sol_mid, Zonotope);

# ------------------------
# Case 2 (short)
# ------------------------

prob = CTRNN_FPA()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=400)
#alg = TMJets(abstol=1e-10, orderT=orderT, orderQ=1);

time_short = @elapsed sol_short = solve(prob, T=0.5, alg=alg);
solz_short = overapproximate(sol_short, Zonotope);


# convert_sets(solz_1, "results/dfp1.mat")
convert_sets(solz_2, "results/fpa2.mat")
convert_sets(solz_mid, "results/fpa_mid.mat")
convert_sets(solz_short, "results/fpa_short.mat")

matwrite("results/fpa_time.mat", Dict("time_short" => time_short,"time_mid" => time_mid, "time_2" => time_2); compress = false)
