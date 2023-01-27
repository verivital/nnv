using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

# How to keep track of time
# eltime3 = @elapsed sol3 = solve(prob, T=Tf, alg=alg, params=param3);

#SUITE = BenchmarkGroup()
#model = "Cartpole"
#cases = ["1"]
#SUITE[model] = BenchmarkGroup()

include("CTRNN_Cartpole.jl")
include("convert_sets.jl")

# ------------------------
# Case 1 (other algo works better, comment this one out)
# ------------------------

#prob = CTRNN_Cartpole()
#alg = TMJets(abstol=1e-10, orderT=8, orderQ=1);

#time_1 = @elapsed sol_1 = solve(prob, T=2.0, alg=alg);
#solz_1 = overapproximate(sol_1, Zonotope);


# ------------------------
# Case 2 (long - good one)
# ------------------------

prob = CTRNN_Cartpole()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=3500)

time_2 = @elapsed sol_2 = solve(prob, T=2.0, alg=alg);
solz_2 = overapproximate(sol_2, Zonotope);

# ------------------------
# Case 3 (mid)
# ------------------------

prob = CTRNN_Cartpole()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=1700)

time_mid = @elapsed sol_mid = solve(prob, T=1.0, alg=alg);
solz_mid = overapproximate(sol_mid, Zonotope);

# ------------------------
# Case 4 (mid)
# ------------------------

prob = CTRNN_Cartpole()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=400)

time_short = @elapsed sol_short = solve(prob, T=0.1, alg=alg);
solz_short = overapproximate(sol_short, Zonotope);


#convert_sets(solz_1, "results/cartpole1.mat")
convert_sets(solz_2, "results/cartpole2.mat") # long
convert_sets(solz_mid, "results/cartpole_mid.mat") # mid
convert_sets(solz_short, "results/cartpole_short.mat") # short

matwrite("results/Cartpole_time.mat", Dict("time_mid" => time_mid, "time_short" => time_short, "time_2" => time_2); compress = false)


