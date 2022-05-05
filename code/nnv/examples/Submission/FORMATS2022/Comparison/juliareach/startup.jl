# deactivate plot GUI, which is not available in Docker
ENV["GKSwstype"] = "100"

# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.add("MAT")

const TARGET_FOLDER = "results"

function main()
    if !isdir(TARGET_FOLDER)
        mkdir(TARGET_FOLDER)
    end

    println("Running HSCC_comparison benchmarks...")

    # Spiral 2D
    println("###\nRunning Spiral 2D benchmark\n###")
    include("models/spiral_benchmark.jl")

    # CTRNN Damped Forced Pendulum
    println("###\nRunning CTRNN DFP benchmark\n###")
    include("models/DFP_benchmark.jl")

    # CTRNN Cartpole
    println("###\nRunning CTRNN Cartpole benchmark\n###")
    include("models/Cartpole_benchmark.jl")

    println("Finished running benchmarks.")

    return 
end

main()


