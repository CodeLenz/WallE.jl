# WallE
Bounding Box Optimizer for large problems where the optimal solution lies on the boundary. The algorithm is a modified Steepest Descent projecting infeasible variables to the boundary of the feasible design space (defined by the side constraints ci and cs). Unless disabled by the user (GC=false), a modified Conjugate Gradient is tried at each iteration to improve convergence. Line search is performed by a modified (projected) Armijo backtracking.

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
    "BETA"=>"HS_proj_sec"
    "GATE"=>"non_conservative"
    "LS_TYPE"=>"Armijo"

```
where NITER is the number of iterations, TOL_NORM is the (relative) tolerance of the norm with respect to the objective function, SHOW enables progress and a summary at the end of the optimization, ARMIJO_C is the constant associated to the expected decrease of the objective function (first Wolfe condition), LS_ALPHA_INI is the initial step in Armijo's Backtracking line search, LS_ALPHA_MIN is the minimum allowable step, and LS_SIGMA is used only by the optional strong line-search check. GC enables the constrained conjugate gradient. If it cannot be used in some iteration, the program automatically switches to steepest descent.

Warning: `LS_STRONG=true` should not be used in production runs. The current implementation is not the standard strong Wolfe condition. After the projected Armijo decrease is satisfied, it checks

```julia
dot(df(xn), Δx / norm(Δx)) >= LS_SIGMA * dot(df(x0), Δx / norm(Δx))
```

where `Δx = xn - x0` is the effective projected displacement. This is a projected curvature check based on the accepted displacement, not the classical strong Wolfe condition `abs(phi'(alpha)) <= c2 * abs(phi'(0))` along an unprojected line-search path. It is kept only for compatibility and experimentation; the default and recommended setting is `LS_STRONG=false`.

BETA selects the conjugate-gradient beta formula:

```julia
    "SD"           # steepest descent fallback
    "FR"           # Fletcher-Reeves
    "PRP"          # Polak-Ribiere-Polyak
    "HS_zero"      # Hestenes-Stiefel restricted to free variables
    "HS_proj_e"    # Hestenes-Stiefel with the original projection correction
    "HS_proj_sec"  # Hestenes-Stiefel with a secant projection correction
```

GATE selects when a conjugate-gradient direction is accepted. `"non_conservative"` is the default and tries GC after the first iteration, accepting it only when the projected descent check succeeds. `"conservative"` keeps the original rule and only tries GC when the active set is stable.

LS_TYPE selects the line search. `"Armijo"` is the default projected Armijo search. `"Wall"` enables the legacy experimental search and disables GC automatically.

The output options are

```julia
    "RESULT"
    "FINI"
    "FOPT"
    "CONVERGED"
    "NORM"
    "COUNTER_ITER"
    "NITER"
    "GC_FRACTION"
    "N_GC_USED"

```
where RESULT is the vector of optimal design variables, FINI is the initial value of the objective function, FOPT is the optimal value of the objective function and CONVERGED is the flag indicating if the optimal solution satisfies first order optimality conditions. NORM is the norm of free positions (not blocked) and COUNTER_ITER
is the effective number of iterations. NITER is an alias for COUNTER_ITER, GC_FRACTION is the fraction of iterations that used a conjugate-gradient direction and N_GC_USED is the number of such iterations.
