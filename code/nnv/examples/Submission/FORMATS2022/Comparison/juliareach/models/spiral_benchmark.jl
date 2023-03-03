using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median



include("spiral2D.jl")
include("convert_sets.jl")

# -----------------------------------
# -----------------------------------

prob = spiralNL1()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=3500)

time_1 = @elapsed sol_1 = solve(prob, T=10.0, alg=alg);
solz_1 = overapproximate(sol_1, Zonotope);

# -----------------------------------
# -----------------------------------

prob = spiralNL2()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=3500)

time_2 = @elapsed sol_2 = solve(prob, T=10.0, alg=alg);
solz_2 = overapproximate(sol_2, Zonotope);

# -----------------------------------
# -----------------------------------

prob = spiralNL3()
orderT = 5
alg = TMJets21a(abstol=1e-10, orderT=orderT, orderQ=1, maxsteps=3500)

time_3 = @elapsed sol_3 = solve(prob, T=10.0, alg=alg);
solz_3 = overapproximate(sol_3, Zonotope);

# -----------------------------------
# -----------------------------------

# alg = GLGM06(δ=1e-7, max_order=1, static=true, dim=2, ngens=4)
# alg=LGG09(δ=6e-4, vars = [x1,x2], sparse=true, cache=false)
# alg = BOX(δ=1e-7, dim = 2)
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


convert_sets(solz_1, "results/spiral1.mat") 
convert_sets(solz_2, "results/spiral2.mat") 
convert_sets(solz_3, "results/spiral3.mat")
convert_sets(solz_4, "results/spiral4.mat") 
convert_sets(solz_5, "results/spiral5.mat") 
convert_sets(solz_6, "results/spiral6.mat")

matwrite("results/spiral_time.mat", Dict("time_1" => time_1, "time_2" => time_2, "time_3" => time_3, "time_4" => time_4, "time_5" => time_5, "time_6" => time_6); compress = false)


# Trying to run the linear models without taylorize
# prob = spiralL_1()

# # time_1 = @elapsed sol_1 = solve(prob, T=10.0, alg=alg);
# time_1 = @elapsed sol_1 = solve(prob, T=10.0)
# solz_1 = box_approximation(sol_1);

# # # -----------------------------------
# # # -----------------------------------

# prob = spiralL_2()

# # # time_5 = @elapsed sol_5 = solve(prob, T=10.0, alg=alg);
# time_2 = @elapsed sol_2 = solve(prob, T=10.0)
# solz_2 = box_approximation(sol_2);

# # # -----------------------------------
# # # -----------------------------------

# prob = spiralL_3()

# # # time_6 = @elapsed sol_6 = solve(prob, T=10.0, alg=alg);
# time_3 = @elapsed sol_3 = solve(prob, T=10.0)
# solz_3 = box_approximation(sol_3);

# # -----------------------------------
# # -----------------------------------


# convert_boxs(solz_1, "results/spiralL1.mat") 
# convert_boxs(solz_2, "results/spiralL2.mat") 
# convert_boxs(solz_3, "results/spiralL3.mat")


# matwrite("results/spiralL_time.mat", Dict("time_1" => time_1, "time_2" => time_2, "time_3" => time_3); compress = false)
