using WallE
using Test



# Run tests
 

println("#"^40, "\n"," "^30,"  Test 1\n","#"^40)

@time include("test_unconstrained.jl")

println("#"^40, "\n"," "^30,"  Test 2\n","#"^40)

@time include("test_constrained.jl")
