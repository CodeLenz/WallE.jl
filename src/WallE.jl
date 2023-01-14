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

  
where NITER is the number of iterations, TOL_NORM is the (relative) 
tolerance of the norm with respect to the objective function, 
SHOW enables a summary at the end of the optimization, 
ARMIJO_C is the constant associated to the expected decrease of the
objective function (first  Wolfe condition), LS_ALPHA_INI is the 
initial step in Armijo's Backtracking line search, LS_ALPHA_MIN is
the minimum allowable step, LS_SIGMA is the parameter associated to
the expected decrese in curvature (second Wolfe condition) that is 
used only if LS_STRONG is true. LS_GC enables the (experimental) 
constrained conjugate gradient. If it cannot be used in some iteration,
the program automatically switch to steepest descent.<br/>

Outputs are returned in another dictionary with keys <br/>

   "RESULT"  <br/>
   "FINI"  <br/>
   "FOPT"  <br/>
   "CONVERGED"  <br/>
   "NORM" <br/>
   "COUNTER_ITER" <br/>

where RESULT is the vector of optimal design variables, 
FINI is the initial value of the objective function,
FOPT is the optimal value of the objective function and 
CONVERGED is the flag indicating if the optimal solution 
satisfies first order optimality conditions. NORM is the 
norm of free positions (not blocked) and COUNTER_ITER
is the effective number of iterations.

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
  nmax_iter,tol_norm,flag_show,armijo_c,cut_factor,α_ini,α_min,σ,STRONG,ENABLE_GC = Check_inputs(f,df,xini,ci,cs,inputs)

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
  Prg = Progress(nmax_iter, 1, "Minimizing objective function...")
  for iter=1:nmax_iter

    # Increment counter
    counter += 1

    # Store function value
    functions[iter] = fn

    # Search direction. Default is Steepest Descent
    d .= -D

    # If we intend to use GC
    if ENABLE_GC && iter>1 && counter_gc <= n && free_x == last_free_x
      flag_gc = GC_projected!(d,last_d,D,last_D,active_r,α_I) 
      if flag_gc
        counter_gc += 1
        used_gc = true
      end
    else 
      counter_gc  = 0
    end

    # Line search
    if inputs["LS_TYPE"]=="Armijo"
       xn, fn, dfn, active_r, active_r_ci, active_r_cs, α, α_I, flag_success = Armijo_Projected!(f,df,x0,fn,D,d,ci,cs,constrained,armijo_c,cut_factor,α_ini,α_min,σ,STRONG)
    elseif inputs["LS_TYPE"]=="Wall"
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
    ProgressMeter.next!(Prg; showvalues = [
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


  # Create the output using the OWall type
  output = Outputs(x0,f0,fn,flag_conv,norm_D,counter,[functions[1:counter], norms[1:counter], steps[1:counter]])

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

    # Hidden option
    LS_TYPE    = inputs["LS_TYPE"]

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

    # Finally, we cannot assert anything on using GC and Wall, so we revert to Steepest
    if LS_TYPE=="Wall" && ENABLE_GC
       println("WallE::Solve::GC cannot be used with Wall LS. Disabling")
       ENABLE_GC = false
    end


    # Return input parameters to the main routine
    return nmax_iter,tol_norm,flag_show,armijo_c,cut_factor,α_ini,α_min,σ,STRONG,ENABLE_GC

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
        
         # Effective α_I
         αI = α - violation/d[i]

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
        
           # Effective α_I
           αI = α - violation/d[i]

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
  # Evaluate the deflection for GC
  #
  #
  #
  function GC_projected!(d::Array{T},last_d::Array{T},
                         D::Array{T},last_D::Array{T},
                         active_r::Array{Int64},α_I::Array{Float64}) where T

  #
  # Lets evaluate the left term of both dot products
  #

  # It starts with the difference in gradient
  y = D .- last_D

  # Loop over last (effectively) projected variables
  @inbounds for r in LinearIndices(active_r)

     # Projected variable
     pos = active_r[r]

     # Correction Assuming An_r = d_r n_r
     y .= y .+ α_I[r]*last_d[pos].*Extract_as_vector(last_d,pos)

  end # r

  # Now we can evaluate beta 
  β = dot(y,D)/dot(y,last_d)

  # Avoid a very unfortunate corner case
  if isnan(β) || β<0.0 
   β = 0.0
  end

  # New search direction
  @inbounds d .= -D .+ β*last_d

  # Let's avoid further problems in the L.S
  # m should be -1 for steepest or close
  # and should be > 0 (or a -δ to avoid problems in the L.S)
  # This is the cos of the angle between d and D
  m = dot(d,D)/(norm(d)*norm(D)) 

  flag_success = true
  if m >=-1E-3 || β==0.0
    flag_success = false
    d .= -D
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
