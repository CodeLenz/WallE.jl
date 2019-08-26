using WallE
using Test



# Run tests
 

println("##########\n  Test 1\n##########")

@time include("test_unconstrained.jl")

println("\n\n##########\n  Test 2\n##########")

@time include("test_constrained.jl")