# WallE
Bounding Box Optimizer for large problems where the optimal solution lies on the boundary. The algorithm is a modified Steepest Descent projecting infeasible variables to the boundary of the feasible design space (defined by the side constraints ci and cs). Unless disabled by the user (GC=false), a modified Conjugate Gradient is tried at each iteration to improve convergence. Line search is performed by a modified (projected) Armijo backtracking and strong conditions can be enabled (LS_STRONG=true). 

To add this package to julia 
```julia
]add https://github.com/CodeLenz/WallE.jl.git
```
To cite this repository:
[![DOI](https://zenodo.org/badge/190200352.svg)](https://zenodo.org/badge/latestdoi/190200352)

Example

```julia
    using WallE

    function f(x) 
        100*(x[2]-x[1]^2)^2+(x[1]-1)^2
    end

       
    function df(x)
        df1 = 2.0*(x[1]-1)-400*x[1]*(x[2]-x[1]^2)
        df2 = 200.0*(x[2]-x[1]^2)
        return [df1 ; df2]
    end

    # Initial point
    x0 = [0.0 ; 3.0]

    # Side constraints
    ci = [-Inf ; 0.5]
    cs = [0.8 ; Inf] 

    # Call optimizer
    options = WallE.Init()
    options["NITER"] = 10_000
    output = WallE.Solve(f,df,x0,ci,cs,options)

    # Recovering solution
    x_opt = output["RESULT"]
    flag_converged = output["CONVERGED"]
    opt_norm = output["NORM"]

```

Default input options are

```julia
   "NITER"=>1000
   "TOL_NORM"=>1E-6
   "SHOW"=>true
   "ARMIJO_C"=>0.1
   "ARMIJO_TAU"=>0.5
   "LS_ALPHA_INI"=>100.0
   "LS_ALPHA_MIN"=>1E-12
   "LS_SIGMA"=>0.9
   "LS_STRONG"=>false
   "GC"=>true

```
where NITER is the number of iterations, TOL_NORM is the (relative) tolerance of the norm with respect to the objective function, SHOW enables a summary at the end of the optimization, ARMIJO_C is the constant associated to the expected decrease of the objective function (first  Wolfe condition), LS_ALPHA_INI is the initial step in Armijo's Backtracking line search, LS_ALPHA_MIN is the minimum allowable step, LS_SIGMA is the parameter associated to the expected decrese in curvature (second Wolfe condition) that is used only if LS_STRONG is true. GC enables the (experimental) constrained conjugate gradient. If it cannot be used in some iteration, the program automatically switch to steepest descent.

The output options are

```julia
    "RESULT"
    "FINI"
    "FOPT"
    "CONVERGED"
    "NORM"
    "COUNTER_ITER"

```
where RESULT is the vector of optimal design variables, FINI is the initial value of the objective function, FOPT is the optimal value of the objective function and CONVERGED is the flag indicating if the optimal solution satisfies first order optimality conditions. NORM is the norm of free positions (not blocked) and COUNTER_ITER
is the effective number of iterations.

