#
# Unconstrained optimization with side constraints
#
module WallE

  using LinearAlgebra, ProgressMeter, Dates

  export Wall_E2

#
# Main function
#
#
  """
  Wall_E2 

  Solve the problem

  Min f(x)

  where x ∈ ℜ^n and x ∈ [ci, cs]. 

  The inputs for this function are:

      f::Function         -> Objective function     -> f(x)->Float64
      df::Function        -> Gradient of f(x)       -> df(x)->Array{Float64,1}
      x0::Array{Float64}  -> Initial point
      ci::Array{Float64}  -> Lower side constraints
      cs::Array{Float64}  -> Upper side constraints

  Optional (with default values) inputs

      nmax_niter::Int64   -> Maximum number of iterations
      tol_norm::Float64   -> Tolerance for the norm
      flag_show::Bool     -> Enable/Disable printing
      cut_factor::Float64 -> Factor to decrease the step length
      α_ini::Float64      -> Initial step length
      α_min::Float64      -> Minimum value for the step length

  Special argument
     
      ENABLE_GC::Bool=false

  Outputs: 

       x0::Array{Float64} -> Optimal point
       f0::Float64        -> Initial objective function
       fn::Float64        -> Final objective function
       flag_conv::Bool    -> Satisfy/Do not satisfy the first order conditions
       [functions, norms, steps] -> lists with values for each iteration
 
  """
function Wall_E2(f::Function,df::Function,
               xini::Array{Float64},
               ci::Array{Float64},
               cs::Array{Float64},
               nmax_iter::Int64=100,
               tol_norm::Float64=1E-6,
               flag_show::Bool=true,
               cut_factor::Float64=0.5,
               α_ini::Float64=10.0,
               α_min::Float64=1E-12; ENABLE_GC::Bool=false)

    
    # TODO -> está inconsistente o estudo dos alphas limites e 
    #         o GC (pois está considerando a direção de steepest)

    # Check the consistence of the inputs
    Check_inputs(f,df,xini,ci,cs,nmax_iter,tol_norm,flag_show,
                 cut_factor,α_ini,α_min,ENABLE_GC)

    # Internal flag to select the GC for constrained/unconstrained problems
    constrained = true
    if ( sum(ci.==-Inf)==length(xini) && sum(cs.==Inf)==length(xini) )
       constrained = false
    end

    #                                                              #
    #                     A little message to our customers        #
    #                     (in case of constrained problems)        #
    #                                                              #
    if ENABLE_GC && constrained 
       println("The actual implementation can lead to a huge improvement in computational time for \nunconstrained problems, but is still in development for constrained problems. Use with care!")
    end



    # Make a copy to unlink initial point with the caller, otherwise 
    # we modify it in the caller, leading to potential problems.
    x0 = copy(xini)

    # Size of the problem
    n = length(x0)

    # List with all variables
    lvar = 1:n

    # First thing..Evaluate initial function value
    f0 = f(x0)
    fn = f0

    # Allocate some vectors we use a lot
    D = zeros(n)
    d = zeros(n)

    # Lists with function values and norms (D)
    functions = zeros(nmax_iter)
    norms     = zeros(nmax_iter)
    steps     = zeros(nmax_iter)

    # Some arrays we whant to show after the main loop
    free_x = Int64[]
    active_r = Int64[]
    active_r_ci = Int64[]
    active_r_cs = Int64[]
    α_I = Float64[]
    delta_m = Int64[]  
    delta_M = Int64[]
    last_x = zeros(n)
    last_d = zeros(n)
    last_D = zeros(n)

  
    # Counter for GC
    counter_gc = 0

    # Norm (Gradient, free positions)
    norm_D = 0.0

    # Flag of convergence
    flag_conv = false

    # Number of effective iterations
    counter = 0
    
    # Step in LS
    α = 0.0

    # We can now enter in the main loop (Steepest)
    tempo = @elapsed  begin
     Prg = Progress(nmax_iter, 1, "Minimizing...")
    for iter=1:nmax_iter

        # Increment counter
        counter += 1

        # Store function value
        functions[iter] = fn

        # Gradient
        D .= df(x0)

        # Search direction. Default is Steepest Descent
        d .= -D

        # If we intend to use GC
        if ENABLE_GC && iter>1 && counter_gc <= n
           flag_gc = GC_projected!(d,last_d,D,last_D,active_r,α_I) 
           if flag_gc
              counter_gc += 1
           end
        else 
           counter_gc  = 0
        end

        # Line search
        xn, fn, active_r, active_r_ci, active_r_cs, α, α_I, flag_success = Armijo_Projected!(f,x0,fn,D,d,ci,cs,constrained)

        # Free positions
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
        last_D          .= D
       

        # Blocked by below. They must be positive
        delta_m = D[active_r_ci]

        # Blocked by above. They must be negative
        delta_M = D[active_r_cs]

        # 


        # We need to fulfil all the first order conditions..
        if iter>2 && norm_D<=tol_norm*(1+abs(fn)) && (all(delta_m .>= 0.0)||isempty(delta_m)) &&
                                                     (all(delta_M .<= 0.0)||isempty(delta_M))

            # Convergence assessed by first order condition. Set the flag and
            # skip the main loop
            flag_conv = true
            break
        end # first order conditions
    
        if !flag_success
          println("WallE2::The solution cannot be improved during the line-search. Bailing out.")
          break
        end


        # Fancy report for the mob :)
        ProgressMeter.next!(Prg; showvalues = [
                          (:Iteration,counter), 
                          (:Counter_gc,counter_gc),
                          (:Norm,norm_D), 
                          (:Target,tol_norm*(1+abs(fn))),
                          (:Objective,fn), 
                          (:(Grad(Max)),maximum(D)), (:(Grad(Min)),minimum(D)),
                          (:Lower,all(delta_m .>= -tol_norm)||isempty(delta_m)),
                          (:Upper,all(delta_M .<= tol_norm)||isempty(delta_M))],
                          valuecolor = :yellow)


    end # iter
    end


     # Final report
    if flag_show
      println("\n********************************************************")
      println("End of the main optimization Loop")
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
      println("Total time             : ", canonicalize(Dates.CompoundPeriod(Dates.Second(floor(Int64,tempo)))))
      println("********************************************************")
    end



    # Return the optimal point, initial and final value of the obj
    # function and the list of objectives/norm and αs for each iteration
    return x0, f0, fn,
           flag_conv, [functions[1:counter], norms[1:counter], steps[1:counter]]

end


#
# Check if the inputs are consistent
#
function Check_inputs(f::Function,df::Function,
               x0::Array{Float64},
               ci::Array{Float64},
               cs::Array{Float64},
               nmax_iter::Int64,
               tol_norm::Float64,
               flag_show::Bool,
               cut_factor::Float64,
               α_ini::Float64,
               α_min::Float64, 
               ENABLE_GC::Bool)


    # Check if the length of x0, ci and cs are the same
    @assert length(x0)==length(ci)==length(cs) "WallE2::Check_inputs:: length of ci, cs and x0 must be the same"

    # Check if x0 is inside the bounds
    @assert  sum(ci .<= x0 .<= cs)==length(x0) "WallE2::Check_inputs:: x0 must be inside the bounds ci and cs" 

    # Check if nmax_iter is positive
    @assert  nmax_iter > 0 "WallE2::Check_inputs:: nmax_iter must be larger than zero "    

    # Check if tol_norm is in (0,1)
    @assert 0.0<tol_norm<1.0 "WallE2::Check_inputs:: tol_norm must be in (0,1)"

    # Check if cut_factor (τ) is in (0,1)
    @assert 0.0<cut_factor<1.0 "WallE2::Check_inputs:: cut_factor must be in (0,1)"

    # Check if α_ini is positive
    @assert 0.0<α_ini "WallE2::Check_inputs:: α_ini must larger than zero"
  
    # Check if α_min is << 1.0 and > 0. At least smaller than α_ini
    @assert  0.0<α_min<α_ini   "WallE2::Check_inputs:: α_min must be in (0,α_ini)"
    
    
end


#
# Return a localization vector 
#
function Localization(n::Int64,pos::Int64)
  v = zeros(n)
  v[pos] = 1.0
  return v
end


#
# Return a vector with just one position 
#
function Extract_as_vector(v::Array{Float64},pos::Int64)
  vv = zero(v)
  vv[pos] = v[pos]
  return vv
end

#
# Return a scalar
#
function Extract_as_scalar(v::Array{Float64},pos::Int64)
  v[pos]
end




#
# Given a point, a search direction and a step
# return the projected point and the list of
# effective blocks
#
function Project(α::Float64,x0::Array{Float64},d::Array{Float64},ci::Array{Float64},cs::Array{Float64})
                  


    # Length 
    n = size(x0,1)

    # Next point, without projections
    xn = x0 + α*d

    #
    # This is the mathematical form of appying the 
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

    for i in LinearIndices(xn)


      # Depending on the seach direction, we can test for lower OR upper
      # violations. If violated, store in the arrays
      if d[i]<0.0
         
        # Possible violation 
        violation = ci[i] - xn[i]

        if violation >= 0.0 
           # Effective α_I
           αI = α - violation/d[i]

           # Keep on the boundary
           xn[i] = ci[i]

           # Store 
           push!(active_r_ci,i)
           push!(active_r,i)
           push!(α_I,αI)
        end   

      elseif d[i]>0.0

        # Possible violation 
        violation =  xn[i] - cs[i]

        if violation >= 0.0 
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
function Armijo_Projected!(f::Function,x0::Array{Float64},
                          f0::Float64,
                          D::Array{Float64},
                          d::Array{Float64},
                          ci::Array{Float64},
                          cs::Array{Float64},
                          constrained::Bool,
                          c::Float64=0.1,
                          τ::Float64=0.5,
                          α_ini::Float64=10.0,
                          α_min::Float64=1E-12)


    # "optimal" value
    fn = 0.0

    # Local vectors
    xn = zero(x0)
    Δx = zero(x0)

    # Local lists to be returned
    active_r = Int64[]
    active_r_ci = Int64[]
    active_r_cs = Int64[]
    α_I = Float64[]

    # Initial step
    α = α_ini

    # Flag (sucess)
    flag_sucess = false

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

        # If we are facing a constrained problem
        # not every initial search direction will
        # lead to an effective projected step. In 
        # this case, we must revert to steepest
        # to make a robust algorithm until we 
        # set a proper direction in GC
        if m>=0.0 && constrained
 
           d .= -D
           xn, active_r, active_r_ci, active_r_cs, α_I = Project(α,x0,d,ci,cs)
           Δx .= xn .- x0 
           m = dot(D,Δx) 

        end 

        if m<0.0 

            # Left side
            fn = f(xn)

            # Rigth side
            rigth = f0 + c*m

            # Armijo condition
            if fn <= rigth 
                flag_sucess= true
                break
            end
        end
        
        # Otherwise, decrease step    
        α = α*τ

        # Check for minimum step
        if α<=α_min
             break
        end

    end #while true

  
    # return 
    return xn, fn, active_r, active_r_ci, active_r_cs, α, α_I, flag_sucess


end #Armijo_Projected




#
# Evaluate the deflection for GC
#
#
#
function GC_projected!(d::Array{Float64},last_d::Array{Float64},
                       D::Array{Float64},last_D::Array{Float64},
                       active_r::Array{Int64},α_I::Array{Float64})

         #
         # Lets evaluate the left term of both dot products
         # 

         # It starts with the difference in gradient
         y = D .- last_D

         # Loop over last (effectivelly) projected variables
         for r in LinearIndices(active_r)

                 # Projected variable
                 pos = active_r[r]

                 # Correction Assuming Ae = d_r e_r
                 y .= y .+ α_I[r]*last_d[pos].*Extract_as_vector(last_d,pos)

         end # r
                
         # Now we can evaluate beta 
         β = dot(y,D)/dot(y,last_d)

         # Avoid a very unfortunate corner case
         if isnan(β) || β<0.0
             β = 0.0
         end

         # New search direction
         d .= -D .+ β*last_d

         # Let's avoid further problems in the L.S
         m = dot(d,D) 
         if m>=0.0
            d .= -D
         end

         # Return a flag true if we succed
         flag_sucess = true
         if m >=0 || β==0.0
            flag_sucess = false
         end

          return flag_sucess

end


end # module
