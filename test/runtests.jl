using WallE
using Test



# Run tests
 

println("#"^80, "\n"," "^40,"  Test 1\n","#"^80)

@time include("test_unconstrained.jl")

#println("#"^80, "\n"," "^40,"  Test 2\n","#"^80)

#@time include("test_constrained.jl")
