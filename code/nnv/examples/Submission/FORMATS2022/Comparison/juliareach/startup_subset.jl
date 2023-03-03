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

    println("Running comparison benchmarks...")

    # Spiral 2D
    println("###\nRunning Spiral 2D benchmark\n###")
    include("models/spiral_benchmark_linear.jl")

    # CTRNN Damped Forced Pendulum
    println("###\nRunning CTRNN FPA benchmark\n###")
    include("models/FPA_benchmark.jl")

    println("Finished running benchmarks.")

    return 
end

main()


