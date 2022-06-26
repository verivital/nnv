using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median



include("spiral2D.jl")
include("convert_sets.jl")

# -----------------------------------
# -----------------------------------

# alg = GLGM06(δ=1e-7, max_order=1, static=true, dim=2, ngens=4)
# alg=LGG09(δ=6e-4, vars = [x1,x2], sparse=true, cache=false)
# alg = BOX(δ=1e-7, dim = 2)
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=3500)
prob = spiralL1()

time_4 = @elapsed sol_4 = solve(prob, T=10.0, alg=alg);
#time_4 = @elapsed sol_4 = solve(prob, T=10.0)
solz_4 = overapproximate(sol_4, Zonotope);

# -----------------------------------
# -----------------------------------

# alg = GLGM06(δ=1e-7, max_order=1, static=true, dim=2, ngens=4)
# alg=LGG09(δ=6e-4, vars = [x1,x2], sparse=true, cache=false)
# alg = BOX(δ=1e-7, dim = 2)
# prob = spiralL1()

# time_4 = @elapsed sol_4 = solve(prob, T=10.0, alg=alg);
# #time_4 = @elapsed sol_4 = solve(prob, T=10.0)
# solz_4 = overapproximate(sol_4, Zonotope);

# # -----------------------------------
# # -----------------------------------

prob = spiralL2()

time_5 = @elapsed sol_5 = solve(prob, T=10.0, alg=alg);
# time_5 = @elapsed sol_5 = solve(prob, T=10.0)
solz_5 = overapproximate(sol_5, Zonotope);

# # -----------------------------------
# # -----------------------------------

prob = spiralL3()

time_6 = @elapsed sol_6 = solve(prob, T=10.0, alg=alg);
# time_6 = @elapsed sol_6 = solve(prob, T=10.0)
solz_6 = overapproximate(sol_6, Zonotope);

# -----------------------------------
# -----------------------------------


convert_sets(solz_4, "results/spiral4.mat") 
convert_sets(solz_5, "results/spiral5.mat") 
convert_sets(solz_6, "results/spiral6.mat")

matwrite("results/spiral_time.mat", Dict("time_4" => time_4, "time_5" => time_5, "time_6" => time_6); compress = false)

