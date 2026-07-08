#
# Unconstrained optimization with side constraints
#
module WallE

  using LinearAlgebra, ProgressMeter, Dates

  export Solve, Init

  # 
  # Generate the dictionary with default values (optional arguments)
  # 
  function Init()

        inputs = Dict()
        push!(inputs,"NITER"=>1000)
        push!(inputs,"TOL_NORM"=>1E-6)
        push!(inputs,"SHOW"=>true)
        push!(inputs,"ARMIJO_C"=>0.1)
        push!(inputs,"ARMIJO_TAU"=>0.5)
        push!(inputs,"LS_ALPHA_INI"=>100.0)
        push!(inputs,"LS_ALPHA_MIN"=>1E-12)
        push!(inputs,"LS_SIGMA"=>0.9)
        push!(inputs,"LS_STRONG"=>false)
        push!(inputs,"GC"=>true)
        push!(inputs,"BETA"=>"HS_proj_sec")
        push!(inputs,"GATE"=>"non_conservative")

        # "Hidden" option :)
        push!(inputs,"LS_TYPE"=>"Armijo")
        
        return inputs

  end

  # 
  # Generate the dictionary with the outputs
  # 
  function Outputs(x::Array{T},
                   f0::T1,fn::T1,flag_conv::Bool,
                   norm_D::Float64,counter::Int64,
                   lists) where{T,T1}

      outputs = Dict()
      push!(outputs,"RESULT"=>x) 
      push!(outputs,"FINI"=>f0) 
      push!(outputs,"FOPT"=>fn)
      push!(outputs,"CONVERGED"=>flag_conv)
      push!(outputs,"NORM"=>norm_D)
      push!(outputs,"COUNTER_ITER"=>counter)
      push!(outputs,"lists"=>lists) 

      return outputs

  end    

       
  

  #
  # Main function
  #
  #
 """
  WallE.Solve 

  Solve the problem

  Min f(x)

  where x ∈ ℜ^n and x ∈ [ci, cs]. 

  The inputs are:

  f::Function         -> Objective function     -> f(x)->Float64  <br/>
  df::Function        -> Gradient of f(x)       -> df(x)->Array{Float64,1}  <br/>
  x0::Array{Float64}  -> Initial point  <br/>
  ci::Array{Float64}  -> Lower side constraints  <br/>
  cs::Array{Float64}  -> Upper side constraints  <br/>

  Optional (with default values) inputs are defined in a dictionary
  with keys (and default values)<br/>

   "NITER"=>1000  <br/>
   "TOL_NORM"=>1E-6  <br/>
   "SHOW"=>true  <br/>
   "ARMIJO_C"=>0.1  <br/>
   "ARMIJO_TAU"=>0.5  <br/>
   "LS_ALPHA_INI"=>100.0  <br/>
   "LS_ALPHA_MIN"=>1E-12  <br/>
   "LS_SIGMA"=>0.9  <br/>
   "LS_STRONG"=>false  <br/>
   "GC"=>true  <br/>
   "BETA"=>"HS_proj_sec"  <br/>
   "GATE"=>"non_conservative"  <br/>
   "LS_TYPE"=>"Armijo"  <br/>

  
where NITER is the number of iterations, TOL_NORM is the (relative) 
tolerance of the norm with respect to the objective function, 
SHOW enables progress and a summary at the end of the optimization, 
ARMIJO_C is the constant associated to the expected decrease of the
objective function (first  Wolfe condition), LS_ALPHA_INI is the 
initial step in Armijo's Backtracking line search, LS_ALPHA_MIN is
the minimum allowable step, LS_SIGMA is the parameter associated to
the expected decrease in curvature (second Wolfe condition) that is 
used only if LS_STRONG is true. GC enables the constrained conjugate
gradient. BETA selects the beta formula (SD, FR, PRP, HS_zero,
HS_proj_e or HS_proj_sec), GATE selects the GC acceptance rule
(conservative or non_conservative), and LS_TYPE selects Armijo or Wall line
search. If GC cannot be used in some iteration, the program
automatically switches to steepest descent.<br/>

Outputs are returned in another dictionary with keys <br/>

   "RESULT"  <br/>
   "FINI"  <br/>
   "FOPT"  <br/>
   "CONVERGED"  <br/>
   "NORM" <br/>
   "COUNTER_ITER" <br/>
   "NITER" <br/>
   "GC_FRACTION" <br/>
   "N_GC_USED" <br/>

where RESULT is the vector of optimal design variables, 
FINI is the initial value of the objective function,
FOPT is the optimal value of the objective function and 
CONVERGED is the flag indicating if the optimal solution 
satisfies first order optimality conditions. NORM is the 
norm of free positions (not blocked) and COUNTER_ITER
is the effective number of iterations. NITER is an alias for
COUNTER_ITER, GC_FRACTION is the fraction of iterations that used a
conjugate-gradient direction and N_GC_USED is the number of such
iterations.

Example:  
  
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
"""
  function Solve(f::Function,df::Function,
                 xini::Array{T},
                 ci=T[],
                 cs=T[],
                 inputs=Dict()) where T



  # Size of the problem
  n = length(xini)

  # If ci or cs are empty, we pass them to ±∞
  if isempty(ci)
    ci = -Inf*ones(T,n)
  end

  if isempty(cs)
    cs = Inf*ones(T,n)
  end

  # If inputs is empty, we use default parameters
  if isempty(inputs)
     inputs = Init()
  end

  # Check and extract the consistence of the inputs
  nmax_iter,tol_norm,flag_show,armijo_c,cut_factor,α_ini,α_min,σ,STRONG,ENABLE_GC,BETA_VAR,GATE_VAR = Check_inputs(f,df,xini,ci,cs,inputs)

  # Internal flag to select the GC for constrained/unconstrained problems
  constrained = true
  if ( sum(ci.==-Inf)==length(xini) && sum(cs.==Inf)==length(xini) )
   constrained = false
  end

  # Just a little remainder to the user
  if STRONG 
     println("STRONG does not improves the solution in our tests. So, the use is not advisable.")
  end
 
  # Used to track limit values of α during optimization
  αs = 0.0
  αi = maxintfloat(Float64)

  # Make a copy to unlink initial point with the caller, otherwise 
  # we modify it in the caller, leading to potential problems.
  x0 = copy(xini)

  # List with all variables
  lvar = 1:n

  # First thing..Evaluate initial function value
  f0 = f(x0)
  fn = f0

  # Allocate some vectors we use a lot
  # We start evaluating ∇f here, since it is evaluated
  # in the LS and returned to this function
  D = df(x0)
  d = zeros(T,n)

  # Lists with function values and norms (D)
  functions = zeros(nmax_iter)
  norms     = zeros(nmax_iter)
  steps     = zeros(nmax_iter)

  # Some arrays we want to show after the main loop
  free_x = Int64[]
  last_free_x = Int64[]
  active_r = Int64[]
  active_r_ci = Int64[]
  active_r_cs = Int64[]
  α_I = Float64[]
  delta_m = Int64[]  
  delta_M = Int64[]
  last_x = zeros(T,n)
  last_d = zeros(T,n)
  last_D = zeros(T,n)


  # Counter for GC
  counter_gc = 0
  used_gc = false
  n_gc_used = 0

  # Norm (Gradient, free positions)
  norm_D = 0.0

  # Flag of convergence
  flag_conv = false

  # Number of effective iterations
  counter = 0

  # Step in LS
  α = α_ini

  # We can now enter in the main loop (Steepest)
  tempo = @elapsed  begin
  Prg = flag_show ? Progress(nmax_iter; dt=1, desc="Minimizing objective function...") : nothing
  for iter=1:nmax_iter

    # Increment counter
    counter += 1

    # Store function value
    functions[iter] = fn

    # Search direction. Default is Steepest Descent
    d .= -D

    # Track whether GC actually produced the direction for this iteration.
    used_gc_iter = false

    # If we intend to use GC
    gate_ok = if GATE_VAR == "non_conservative"
      ENABLE_GC && iter>1 && counter_gc <= n
    else
      ENABLE_GC && iter>1 && counter_gc <= n && free_x == last_free_x
    end

    if gate_ok
      flag_gc = GC_projected!(d,last_d,D,last_D,active_r,α_I,
                              x0,last_x,BETA_VAR,GATE_VAR,free_x,
                              ci,cs,α_ini) 
      if flag_gc
        counter_gc += 1
        used_gc = true
        used_gc_iter = true
      end
    else 
      counter_gc  = 0
    end
    used_gc_iter && (n_gc_used += 1)

    # Line search
    if get(inputs,"LS_TYPE","Armijo")=="Armijo"
       xn, fn, dfn, active_r, active_r_ci, active_r_cs, α, α_I, flag_success = Armijo_Projected!(f,df,x0,fn,D,d,ci,cs,constrained,armijo_c,cut_factor,α_ini,α_min,σ,STRONG)
    elseif get(inputs,"LS_TYPE","Armijo")=="Wall"
       xn, fn, dfn, active_r, active_r_ci, active_r_cs, α, α_I, flag_success = Wall_Seach_Projected!(f,df,x0,fn,D,d,ci,cs,constrained,α,armijo_c,cut_factor,α_ini,α_min,σ,STRONG)
    else
       error("WallE::Solve::Hidden option LS_TYPE should be Armijo or Wall")   
    end

    # keep track of αs and αi
    αs = max(αs,α)
    αi = min(αi,α)

    # Copy the new derivative and store the old one
    last_D   .= D
    D        .= dfn

    # Free positions
    last_free_x = copy(free_x)
    free_x = filter(x-> !(x in active_r),lvar)

    # Norm of free positions
    norm_D = norm(D[free_x])  

    # Store the norm (d)
    norms[iter] = norm_D

    # Store the step
    steps[iter] = α

    # Rollover Bethoven
    last_x          .= x0
    last_d          .= d
    x0              .= xn

    # Blocked by below. They must be positive
    delta_m = D[active_r_ci]

    # Blocked by above. They must be negative
    delta_M = D[active_r_cs]

    # Breaking condition when function doesn't improve
    if !flag_success 
        printstyled("\nWallE.Solve::The solution cannot be improved during the line-search. ", color=:red)
        if  norm_D<=tol_norm*(1+abs(fn)) && (all(delta_m .>= 0.0)||isempty(delta_m)) &&
                                            (all(delta_M .<= 0.0)||isempty(delta_M))
          printstyled("\nWallE.Solve::But first order conditions are satisfied.", color=:green)

          flag_conv = true 
        else
          printstyled("\nWallE.Solve::Not all first order conditions are satisfied, proceed with care. ", color=:red)
        end
        break
    end

    # We need to fulfil all the first order conditions..
    if flag_success && iter>2 && norm_D<=tol_norm*(1+abs(fn)) && (all(delta_m .>= 0.0)||isempty(delta_m)) &&
       (all(delta_M .<= 0.0)||isempty(delta_M))
      # Convergence assessed by first order condition. Set the flag and
      # skip the main loop
      flag_conv = true
      break
    end # first order conditions

      
    # Fancy report for the mob :)
    flag_show && ProgressMeter.next!(Prg; showvalues = [
                      (:Iteration,counter), 
                      (:Counter_gc,counter_gc),
                      (:Enable_GC,ENABLE_GC),
                      (:GC,used_gc),
                      (:Norm,norm_D), 
                      (:Target,tol_norm*(1+abs(fn))),
                      (:"Current Step",α),
                      (:"Smaller Step",αi),
                      (:"Larger  Step",αs),
                      (:Objective,fn), 
                      (:ci,length(active_r_ci)),
                      (:cs,length(active_r_cs)),
                      (:(Grad(Max)),maximum(D)), (:(Grad(Min)),minimum(D)),
                      (:Lower,all(delta_m .>= -tol_norm)||isempty(delta_m)),
                      (:Upper,all(delta_M .<= tol_norm)||isempty(delta_M))],
                      valuecolor = :yellow)


    end # iter
    end # block for timing


  # Final report
  if flag_show
    println("\n********************************************************")
    println("End of the main optimization Loop")
    println("Method                 : ",ifelse(ENABLE_GC,"Conjugate gradient","Steepest descent"))
    if STRONG
     println("Using strong L.S       : Yes, with $σ")
    end
     if ENABLE_GC 
       println("GC                     : ",ifelse(used_gc,"used","not used"))
     end
     println("Type of problem        : ",ifelse(constrained,"constrained","unconstrained"))
     println("Number of variables    : $(n)")
     println("Initial objective      : ", f0)
     println("Final objective        : ", fn)
     if f0!=0.0 && fn!=0.0
      println("% of minimization.     : ", 100*(fn-f0)/f0)
    end
    println("Free variables         : ", length(free_x))
    println("Blocked variables      : ", length(active_r),": ",  length(active_r_ci)," for lower bound ",length(active_r_cs)," for upper bound")
    println("Number of iterations   : ", counter , " of ",nmax_iter)
    println("First order conditions : ", flag_conv, " ", all(delta_m .>= -tol_norm)||isempty(delta_m),
                                      " ", all(delta_M .<=  tol_norm)||isempty(delta_M))
    println("Norm(free positions)   : ", norm_D," Reference ",tol_norm*(1+abs(fn)))
    println("Smaller step           : ", αi)
    println("Larger  step           : ", αs)
    println("Total time             : ", canonicalize(Dates.CompoundPeriod(Dates.Second(floor(Int64,tempo)))))
    println("********************************************************")
  end


  # Create the output dictionary
  output = Outputs(x0,f0,fn,flag_conv,norm_D,counter,[functions[1:counter], norms[1:counter], steps[1:counter]])
  push!(output,"NITER"=>counter)
  push!(output,"GC_FRACTION"=> counter>0 ? n_gc_used/counter : 0.0)
  push!(output,"N_GC_USED"=>n_gc_used)

  # Return the optimal point, initial and final value of the obj
  # function and the list of objectives/norm and αs for each iteration
  return output

  end


  #
  # Check if the inputs are consistent
  #
  function Check_inputs(f::Function,df::Function,
                        x0::Array{T},
                        ci::Array{T},
                        cs::Array{T},
                        inputs::Dict) where T

                 

    #
    # First thing is to extract the input parameters 
    #
    nmax_iter  = inputs["NITER"]
    tol_norm   = inputs["TOL_NORM"]
    flag_show  = inputs["SHOW"]
    armijo_c   = inputs["ARMIJO_C"]
    cut_factor = inputs["ARMIJO_TAU"]
    α_ini      = inputs["LS_ALPHA_INI"]
    α_min      = inputs["LS_ALPHA_MIN"]
    σ          = inputs["LS_SIGMA"]
    STRONG     = inputs["LS_STRONG"]
    ENABLE_GC  = inputs["GC"]
    BETA_VAR   = get(inputs,"BETA","HS_proj_sec")
    GATE_VAR   = get(inputs,"GATE","non_conservative")

    # Hidden option
    LS_TYPE    = get(inputs,"LS_TYPE","Armijo")

    # Check if the length of x0, ci and cs are the same
    @assert length(x0)==length(ci)==length(cs) "Solve::Check_inputs:: length of ci, cs and x0 must be the same"

    # Check if x0 is inside the bounds
    @assert  sum(ci .<= x0 .<= cs)==length(x0) "Solve::Check_inputs:: x0 must be inside the bounds ci and cs" 

    # Check if nmax_iter is positive
    @assert  nmax_iter > 0 "Solve::Check_inputs:: NITER must be larger than zero "

    # Check if tol_norm is in (0,1)
    @assert 0.0<tol_norm<1.0 "Solve::Check_inputs:: TOL_NORM must be in (0,1)"

    # Check if armijo_c is in (0,0.5)
    @assert 0.0<armijo_c<0.5 "Solve::Check_inputs:: ARMIJO_C must be in (0,0.5)"

    # Check if cut_factor (τ) is in (0,1)
    @assert 0.0<cut_factor<1.0 "Solve::Check_inputs:: ARMIJO_TAU must be in (0,1)"

    # Check if α_ini is positive
    @assert 0.0<α_ini "Solve::Check_inputs:: LS_ALPHA_INI must larger than zero"

    # Check if α_min is << 1.0 and > 0. At least smaller than α_ini
    @assert  0.0<α_min<α_ini   "Solve::Check_inputs:: LS_ALPHA_MIN must be in (0,LS_ALPHA_INI)"

    # Check if σ is in armijo_c <= \sigma < 1.0
    @assert armijo_c <= σ < 1.0 "Solve::Check_inputs:: LS_SIGMA must be in [ARMIJO_C,1)"

    # Check the hidden option
    @assert (LS_TYPE=="Armijo" || LS_TYPE=="Wall") "Solve::Check_inputs:: LS_TYPE must be Armijo OR Wall"

    @assert BETA_VAR in ("SD","FR","PRP","HS_zero","HS_proj_e","HS_proj_sec") "Solve::Check_inputs:: BETA must be SD, FR, PRP, HS_zero, HS_proj_e OR HS_proj_sec"

    @assert GATE_VAR in ("conservative","non_conservative") "Solve::Check_inputs:: GATE must be conservative OR non_conservative"

    # Finally, we cannot assert anything on using GC and Wall, so we revert to Steepest
    if LS_TYPE=="Wall" && ENABLE_GC
       println("WallE::Solve::GC cannot be used with Wall LS. Disabling")
       ENABLE_GC = false
    end


    # Return input parameters to the main routine
    return nmax_iter,tol_norm,flag_show,armijo_c,cut_factor,α_ini,α_min,σ,STRONG,ENABLE_GC,BETA_VAR,GATE_VAR

  end


  #
  # Return a localization vector 
  #
  function Localization(n::Int64,pos::Int64)
    v = zeros(n)
    @inbounds v[pos] = 1.0
    return v
  end


  #
  # Return a vector with just one position 
  #
  function Extract_as_vector(v::Array{T},pos::Int64) where T
    vv = zero(v)
    @inbounds vv[pos] = v[pos]
    return vv
  end

  #
  # Return a scalar
  #
  function Extract_as_scalar(v::Array{T},pos::Int64) where T
    @inbounds v[pos]
  end




  #
  # Given a point, a search direction and a step
  # return the projected point and the list of
  # effective blocks
  #
  function Project(α::Float64,x0::Array{T},d::Array{T},ci::Array{T},cs::Array{T}) where T



    # Length 
    n = size(x0,1)

    # Next point, without projections
    xn = x0 .+ α*d

    #
    # This is the mathematical form of applying the 
    # projections, as explained in the companion text.
    # 
    #
    # For each direction we look for violations, apply the corrections
    # and evaluate effective step in this direction
    #
    active_r_ci = Int64[]
    active_r_cs = Int64[]
    active_r    = Int64[]
    α_I = Float64[]

    @inbounds for i in LinearIndices(xn)


      # Depending on the search direction, we can test for lower OR upper
      # violations. If violated, store in the arrays
      if d[i] < zero(T)

        # Possible violation 
        violation = ci[i] - xn[i]

        if violation >= zero(T) 
        
         # Infeasible part of the step after hitting the lower bound.
         αS = (ci[i] - x0[i]) / d[i]
         αI = max(0.0, α - αS)

         # Keep on the boundary
         xn[i] = ci[i]

         # Store 
         push!(active_r_ci,i)
         push!(active_r,i)
         push!(α_I,αI)
       end   

     elseif d[i] > zero(T)

        # Possible violation 
        violation =  xn[i] - cs[i]

        if violation >= zero(T)
        
           # Infeasible part of the step after hitting the upper bound.
           αS = (cs[i] - x0[i]) / d[i]
           αI = max(0.0, α - αS)

           # Keep on the boundary
           xn[i] = cs[i]

           # Store 
           push!(active_r_cs,i)
           push!(active_r,i)
           push!(α_I,αI)
       end   

    end

  end

  return xn, active_r, active_r_ci, active_r_cs, α_I

  end # Project




  #
  # Modified Line Search (Armijo). Search direction is modified  (scaled)
  # in this subroutine
  #
  function Armijo_Projected!(f::Function,df::Function,x0::Array{T},
                             f0::Float64,
                             D::Array{T},
                             d::Array{T},
                             ci::Array{T},
                             cs::Array{T},
                             constrained::Bool,
                             c::Float64=0.1,
                             τ::Float64=0.5,
                             α_ini::Float64=1.0,
                             α_min::Float64=1E-12,
                             σ::Float64=0.95,
                             strong::Bool=true) where T


  # "optimal" value
  fn = f0

  # Local vectors
  xn = copy(x0)
  Δx = zero(x0)

  # Local lists to be returned
  active_r = Int64[]
  active_r_ci = Int64[]
  active_r_cs = Int64[]
  α_I = Float64[]

  # Initial step
  α = α_ini

  # Derivative on (next) point
  dfn = copy(D) #zero(x0)

  # Flag (success)
  flag_success = false

  # Normalize search direction
  d .= d./norm(d)    

  # Main Loop
  while true

    # Candidate point (xn)
    xn, active_r, active_r_ci, active_r_cs, α_I = Project(α,x0,d,ci,cs)

    # Effective delta x
    Δx .= xn .- x0 

    # Projection can fully block a direction, producing no effective step.
    if norm(Δx)==0.0
      α = α*τ
      if α<=α_min
        break
      end
      continue
    end

    # Effective slope
    m = dot(D,Δx) 

    # Normalized slope (to help set a proper limit to skip GC) 
    nm = m/(norm(D)*norm(Δx))

    # If we are facing a constrained problem
    # not every initial search direction will
    # lead to an effective projected step. In 
    # this case, we must revert to steepest
    # to make a robust algorithm until we 
    # set a proper direction in GC
    if nm>=-1E-3 && constrained

       d .= -D / norm(D)
       xn, active_r, active_r_ci, active_r_cs, α_I = Project(α,x0,d,ci,cs)
       Δx .= xn .- x0 
       m = dot(D,Δx) 

    end 

    # We just test for this point if the slope is negative
    if m < 0.0 

      # Left side
      fn = f(xn)

      # Rigth side
      right = f0 + c*m

      # First Wolfe condition
      if fn <= right

        # We evaluate derivative anyway, since we 
        # must return it to the main function
        dfn .= df(xn)

        # Check if we must evaluate second (strong) Wolfe condition
        Δnorm = Δx / norm(Δx)
        if !strong || (strong && dot(dfn,Δnorm) >= σ*dot(D,Δnorm)) 
            flag_success = true
            break
        end

      end #fn <= right
    end # m<=0

    # Otherwise, decrease step    
    α = α*τ

    # Check for minimum step
    if α<=α_min
      break
    end

  end #while true


  # return 
  return xn, fn, dfn, active_r, active_r_ci, active_r_cs, α, α_I, flag_success


  end #Armijo_Projected


  #
  # Check the projected descent condition for the current point and direction.
  #
  function Projected_descent(D::Array{T},d::Array{T},x::Array{T},
                             ci::Array{T},cs::Array{T},
                             α::Float64) where T

    dnorm = norm(d)
    if dnorm == 0.0
      return false
    end

    direction = d ./ dnorm
    xn, active_r, active_r_ci, active_r_cs, α_I = Project(α,x,direction,ci,cs)

    lhs = α*dot(D,direction)
    rhs = 0.0
    @inbounds for r in LinearIndices(active_r)
      pos = active_r[r]
      rhs += α_I[r]*direction[pos]*D[pos]
    end

    return lhs < 0.0 && lhs <= rhs

  end


  #
  # Evaluate the deflection for GC
  #
  function GC_projected!(d::Array{T},last_d::Array{T},
                         D::Array{T},last_D::Array{T},
                         active_r::Array{Int64},α_I::Array{Float64},
                         x::Array{T},last_x::Array{T},
                         BETA_VAR::String,GATE_VAR::String,
                         free_x::Array{Int64},
                         ci::Array{T},cs::Array{T},
                         α_ini::Float64) where T

  n = length(D)
  β = 0.0

  if BETA_VAR == "SD"

    β = 0.0

  elseif BETA_VAR == "FR"

    den = dot(last_D,last_D)
    β = den>0.0 ? dot(D,D)/den : 0.0

  elseif BETA_VAR == "PRP"

    den = dot(last_D,last_D)
    β = den>0.0 ? dot(D,D .- last_D)/den : 0.0

  elseif BETA_VAR == "HS_zero"

    y  = D .- last_D
    yf = zeros(T,n)
    df = zeros(T,n)
    @inbounds for i in free_x
      yf[i] = y[i]
      df[i] = last_d[i]
    end
    den = dot(yf,df)
    β = abs(den)>1E-300 ? dot(yf,D)/den : 0.0

  elseif BETA_VAR == "HS_proj_e"

    y = D .- last_D
    @inbounds for r in LinearIndices(active_r)
      pos = active_r[r]
      y .= y .+ α_I[r]*last_d[pos].*Extract_as_vector(last_d,pos)
    end
    den = dot(y,last_d)
    β = abs(den)>1E-300 ? dot(y,D)/den : 0.0

  elseif BETA_VAR == "HS_proj_sec"

    y = D .- last_D
    @inbounds for r in LinearIndices(active_r)
      pos = active_r[r]
      dxr = x[pos] - last_x[pos]
      Anr = abs(dxr)>1E-300 ? (D[pos]-last_D[pos])/dxr : 0.0
      y .= y .+ α_I[r]*last_d[pos]*Anr.*Extract_as_vector(last_d,pos)
    end
    den = dot(y,last_d)
    β = abs(den)>1E-300 ? dot(y,D)/den : 0.0

  end

  if isnan(β) || β<0.0
    β = 0.0
  end

  @inbounds d .= -D .+ β*last_d

  flag_success = true

  if GATE_VAR == "non_conservative"

    if β==0.0 || !Projected_descent(D,d,x,ci,cs,α_ini)
      flag_success = false
      d .= -D
    else
      d ./= norm(d)
    end

  else

    m = dot(d,D)/(norm(d)*norm(D))
    if m >=-1E-3 || β==0.0
      flag_success = false
      d .= -D
    end

  end

  return flag_success

  end




  ##################################################################################
  ############################### HIDDEN FUNCTION ##################################
  ##################################################################################


  #
  # Not fair L.S. Search direction is modified  (scaled)
  # in this subroutine
  #
  function Wall_Seach_Projected!(f::Function,df::Function,x0::Array{T},
                                 f0::Float64,
                                 D::Array{T},
                                 d::Array{T},
                                 ci::Array{T},
                                 cs::Array{T},
                                 constrained::Bool,
                                 last_α::Float64,
                                 c::Float64=0.1,
                                 τ::Float64=0.5,
                                 α_ini::Float64=1.0,
                                 α_min::Float64=1E-12,
                                 σ::Float64=0.95,
                                 strong::Bool=true) where T


  # "optimal" value
  fn = f0

  # Local vectors
  xn = copy(x0)

  # Local lists to be returned
  active_r = Int64[]
  active_r_ci = Int64[]
  active_r_cs = Int64[]
  α_I = Float64[]

  # Initial step
  α = 2*last_α
  if α>α_ini
     α = α_ini
  end

  # Derivative on (next) point
  dfn = copy(D) 

  # Flag (success)
  flag_success = false

  # Normalize search direction
  d .= d./norm(d)    

  # Main Loop
  while true

    # Candidate point (xn)
    xn, active_r, active_r_ci, active_r_cs, α_I = Project(α,x0,d,ci,cs)

    # Function at this next point
    fn = f(xn)

    # if improved we bail out
    if fn < f0 

        # We are done here
        flag_success = true 

        # We evaluate derivative anyway, since we 
        # must return it to the main function
        dfn .= df(xn)

        # Skip the loop
        break

    end 

    # Otherwise, decrease step    
    α = α*τ

    # Check for minimum step
    if α<=α_min
      break
    end

  end #while true


  # return 
  return xn, fn, dfn, active_r, active_r_ci, active_r_cs, α, α_I, flag_success


  end # Dirty LS


end # module
