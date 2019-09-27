# WallE
Bounding Box Optimizer for large problems where the optimal solution lies on the boundary.

BExample

```julia
    using WallE

    # Objective function
    function f(x) 
       (x[1]+2*x[2]-7)^2 + (2*x[1]+x[2]-5)^2 
    end

    # Gradient 
    function df(x)
        df1 = 2*(2*x[2]+x[1]-7)+4*(x[2]+2*x[1]-5)
        df2 = 4*(2*x[2]+x[1]-7)+2*(x[2]+2*x[1]-5)
        return [df1 ; df2]
    end

    # Initial point
    x0 = [-5.0 ; -5.0]

    # Side constraints
    ci =  -Inf*ones(2)
    cs =  [0.5 ; 2.0]

    # Call optimizer
    options = WallE.Init()
    options["NITER"] = 1000
    output = WallE.Solve(f,df,x0,ci,cs,options)

    # Recovering solution
    x_opt = output["RESULT"]
    flag_converged = output["CONVERGED"]

```

Default input options are

```
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

and the output options are

```
    "RESULT"
    "FINI"
    "FOPT"
    "CONVERGED"
```
