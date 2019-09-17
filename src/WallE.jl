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
               α_min::Float64=1E-12; ENABLE_GC::Bool=false, ENABLE_QN::Bool=false)

    
    # TODO -> está inconsistente o estudo dos alphas limites e 
    #         o GC (pois está considerando a direção de steepest)

    # Check the consistence of the inputs
    Check_inputs(f,df,xini,ci,cs,nmax_iter,tol_norm,flag_show,
                 cut_factor,α_ini,α_min,ENABLE_GC,ENABLE_QN)

    #                                                              #
    #                     A little message to our customers        #
    #                     (in case of constrained problems)        #
    #                                                              #
    if ENABLE_GC && ( sum(ci.==-Inf)!=length(xini) || sum(cs.==Inf)!=length(xini) ) 
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
    blocked_x = Int64[]
    I_blk_m = Int64[]
    I_blk_M = Int64[]
    delta_m = Int64[]  
    delta_M = Int64[]

    # We must keep track of some variables to evaluate β
    α_limit          = Float64[]
    last_x           = zeros(n)
    list_r           = Int64[]
    last_free_x      = Int64[]
    last_d           = zeros(n)
    last_D           = zeros(n)
    last_α           = 0.0
    last_α_limit     = Float64[]
    last_list_r      = Int64   

    # Norm (Gradient, free positions)
    norm_D = 0.0

    # Flag of convergence
    flag_conv = false

    # Number of effective iterations
    counter = 0
    
    # Number of iterations with GC
    counter_gc = 0

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

        # Search direction
        d .= -D

        # If some variable is at the boundary and 
        # there is a tendency to violate, we can 
        # supress the variable from the search direction
        I_blk_m, I_blk_M = Project_d!(x0,d,ci,cs)

        # Join the list of blocked variables
        blocked_x = sort(vcat(I_blk_m,I_blk_M))

        # List of free variables
        free_x = filter(x-> !(x in blocked_x),lvar)

        # Make a copy here  
        last_α           = α
        last_α_limit     = copy(α_limit)
        last_list_r      = copy(list_r)   

        # Norm (to scale search direction)
        norm_d = norm(d)

        # Norm of the free positions
        norm_D = norm(D[free_x])

        # Store the norm (d)
        norms[iter] = norm_d
  
        # Scale
        if length(blocked_x)!=n
           d .= d./norm_d
        end

        # Find limit alphas and constrained elements
        α_limit, list_r = Find_limit_alphas(x0,d,ci,cs)

        # Line search
        α, xn, fn, last_d, counter_gc, flag_sucess = Armijo_Projected_GC(f,x0,
                                                             fn,D,last_D,
                                                             d,last_d,ci,cs,α_limit,last_α_limit,
                                                             list_r,last_list_r,
                                                             iter,counter_gc,ENABLE_GC)

        # Store the step
        steps[iter] = α

        # Rollover Bethoven
        last_x          .= x0
        x0              .= xn
        last_free_x      = copy(free_x)
        last_D          .= D
       

        # Blocked by below. They must be positive
        delta_m = D[I_blk_m]

        # Blocked by above. They must be negative
        delta_M = D[I_blk_M]

        # We need to fulfil all the first order conditions..
        if norm_D<=tol_norm*(1+abs(fn)) && (all(delta_m .>= 0.0)||isempty(delta_m)) &&
                                           (all(delta_M .<= 0.0)||isempty(delta_M))

            # Convergence assessed by first order condition. Set the flag and
            # skip the main loop
            flag_conv = true
            break
        end # first order conditions
    
        # Fancy report for the mob :)
        ProgressMeter.next!(Prg; showvalues = [
                          (:Iteration,counter), 
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
      println("Blocked variables      : ", length(blocked_x),": ",  length(I_blk_m)," for lower bound ",length(I_blk_M)," for upper bound")
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
               ENABLE_GC::Bool,
               ENABLE_QN::Bool)


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
    
    # We cannot set both ENABLE_GC and ENABLE_QN at the same time
    @assert (!(ENABLE_GC&ENABLE_QN)) "WallE2::Check_inputs:: You can only enable GC OR QN"

end

#
# Find the maximum feasible step for each variable.
#
function Find_limit_alphas(x0::Array{Float64},d::Array{Float64},
                           ci::Array{Float64},cs::Array{Float64})

    #
    # Empty lists with limit steps and variables that can be 
    # blocked from below or from above
    #
    alpha_limit = Float64[]
    list_r      = Int64[]

    # For each variable, find the limit step
    for i in LinearIndices(x0)

        if d[i]<0.0
           
           s_min = (ci[i]-x0[i])/d[i]
           push!(alpha_limit,s_min)
           push!(list_r,i)

        elseif d[i]>0.0   

           s_max = (cs[i]-x0[i])/d[i]
           push!(alpha_limit,s_max)
           push!(list_r,i)

        end # d<0

    end #i

    # Return the lists
    return alpha_limit, list_r

end # Find_limit_alphas

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
function Project!(α::Float64,x0::Array{Float64},xn::Array{Float64},
                  d::Array{Float64},ci::Array{Float64},cs::Array{Float64},
                  alpha_limit::Array{Float64}, 
                  list_r::Array{Int64},dirty::Bool=true)


    # Length 
    n = size(x0,1)

    # Next point, without projections
    xn .= x0 + α*d


    # Fast and dirty projection. The result is the same 
    # as the mathematical form, commented bellow
    if dirty
       xn .= max.(ci,min.(cs,xn))

    else
      #
      # This is the mathematical form of appying the 
      # projections, as explained in the companion text.
      # 
      #
      # for each candidate constraint and for α > alpha_limit_r
      # apply the projection
      for r in LinearIndices(list_r)

          # Effective step
          α_eff = max(0.0, α - alpha_limit[r])

          if α_eff > 0.0

              # This is the mathematical projection to the boundary
              xn .= xn .- (α_eff  .* Extract_as_vector(d,list_r[r]))
   

          end 

      end
    end #dirty

end # Project


#
# Evaluate if we can simplify the search direction. This can 
# also be done if the variable is "very close" to the Wall
# and the search direction leads to the Wall.
#
function Project_d!(x0::Array{Float64},d::Array{Float64},
                   ci::Array{Float64},cs::Array{Float64})

     # Constrained sets
     I_blk_m = Int64[]
     I_blk_M = Int64[] 

     for i in LinearIndices(x0)

         # Tolerance for the bound
         # TODO
         # This should be linked to minimum step in LS.
         tol = 1E-6
         if cs[i]!=Inf && ci[i]!=-Inf
            tol = (cs[i]-ci[i])/1E6 
         end

         if d[i]<0 && x0[i]<=(ci[i] + tol)
          x0[i] = ci[i]
          d[i] = 0.0
          push!(I_blk_m,i)
         elseif d[i]>0 && x0[i]>=(cs[i] - tol)
          x0[i] = cs[i]
          d[i] = 0.0
          push!(I_blk_M,i)
         end  
     end #i

    return I_blk_m, I_blk_M
end # Project_d


#
# Modified Line Search (Armijo)
#
function Armijo_Projected(f::Function,x0::Array{Float64},
                          f0::Float64,
                          D::Array{Float64},
                          d::Array{Float64},
                          ci::Array{Float64},
                          cs::Array{Float64},
                          alpha_limit::Array{Float64},
                          list_r::Array{Int64},
                          c::Float64=0.1,
                          τ::Float64=0.5,
                          α_ini::Float64=10.0,
                          α_min::Float64=1E-12)


    # "optimal" value
    fn = 0.0

    # Local vectors
    xn = zero(x0)
    Δx = zero(x0)

    # Initial step
    α = α_ini

    # Flag (sucess)
    flag_sucess = false

    # Main Loop
    while true

        # Candidate point (xn)
        Project!(α,x0,xn,d,ci,cs,alpha_limit,list_r)
 
        # Effective delta x
        Δx .= xn .- x0 

        # Left side
        fn = f(xn)

        # Rigth side
        rigth = f0 + c*dot(D,Δx)

        # Armijo condition
        if fn <= rigth 
            flag_sucess= true
            break
        else
          α = α*τ
          if α<=α_min
             break
          end 
        end

    end #while true

    # return effective step, next point, function value at this
    # point and the flag
    return α, xn, fn, flag_sucess


end #Armijo_Projected



#
# Modified Line Search (Armijo) - with GC
#
function Armijo_Projected_GC(f::Function,x0::Array{Float64},
                             f0::Float64,
                             D::Array{Float64},
                             last_D::Array{Float64},
                             d::Array{Float64},
                             last_d::Array{Float64},
                             ci::Array{Float64},
                             cs::Array{Float64},
                             alpha_limit::Array{Float64},
                             last_alpha_limit::Array{Float64},
                             list_r::Array{Int64},
                             last_list_r::Array{Int64},
                             iter::Int64,
                             counter_gc::Int64,
                             ENABLE_GC::Bool=true,
                             c::Float64=0.1,
                             τ::Float64=0.5,
                             α_ini::Float64=10.0,
                             α_min::Float64=1E-12)


    # "optimal" value
    fn = 0.0

    # Problem size
    n = length(x0)

    # Local vectors
    xn = zeros(n)
    Δx = zeros(n)
    d_eff = zeros(n)

    # Initial step
    α = α_ini

    # Flag (sucess)
    flag_sucess = false

    # Main Loop
    while true

      #
      # Evaluate deflection (limit β)
      #
      if ENABLE_GC && iter > 1 && counter_gc <= n && list_r==last_list_r


          cima = α*dot(D,D)
          baixo = α*dot(last_D,last_D)

          if  length(list_r)>0 

              for r in LinearIndices(list_r)

                # Effective step
                α_eff = max(0.0, α - alpha_limit[r])
                lα_eff = max(0.0, α - last_alpha_limit[r])

                if α_eff > 0.0 || lα_eff > 0.0

                  # Correct both terms
                  D_pos  = Extract_as_scalar(D,list_r[r])
                  lD_pos = Extract_as_scalar(last_D,list_r[r])
                  cima = cima - α_eff  * D_pos^2
                  baixo = baixo -  lα_eff  * lD_pos^2
          
                end 
              end
          end

      # Deflection
      β = max(0.0,cima/baixo)

      

    
    else # Steepest Descent

      β = 0.0
      
    end 
  
        # Search direction
        d_eff .= -D .+ β*last_d

     
        # Normalize
        d_eff .= d_eff / norm(d_eff)

        # Candidate point (xn)
        Project!(α,x0,xn,d_eff,ci,cs,alpha_limit,list_r)
 
        # Effective delta x
        Δx .= xn .- x0 

        # Is this a minimizer direction ?
        m = dot(D,Δx)

        if m>=0 
           println("Not minimizing...") 
        end

        # Left side
        fn = f(xn)

        # Rigth side
        rigth = f0 + c*dot(D,Δx)

        # Armijo condition
        if fn <= rigth 
            flag_sucess= true
            break
        else
          α = α*τ
          if α<=α_min
             break
          end 
        end

    end #while true

    if ENABLE_GC && counter_gc<n
      counter_gc += 1
    elseif ENABLE_GC && counter_gc>=n
      counter_gc = 0
    end

    # return effective step, next point, function value at this
    # point and the flag
    return α, xn, fn, d_eff, counter_gc, flag_sucess


end #Armijo_Projected_GC






#
# Evaluate the deflection for GC
#
# This is not working, and it is under heavy 
# development. You have been warned !
#
#
function GC_projected!(d::Array{Float64},D::Array{Float64},last_D::Array{Float64},
                       last_α::Float64,last_list_r::Array{Int64},last_α_limit::Array{Float64},
                       last_d::Array{Float64},free_x::Array{Int64},
                       x::Array{Float64},
                       last_x::Array{Float64})

         #
         # Lets evaluate the left term of both dot products
         # 

         # It starts with the difference in gradient
         y = D .- last_D

         # Denominator of the Rank-1 update
         ys = dot(y, x.-last_x) + 1E-8

         # Loop over last (effectivelly) projected variables
         L = copy(y)
         for r in LinearIndices(last_list_r)

             # Effective alfa 
             effective_α = max(0.0, last_α - last_α_limit[r])

             # If positive, it was projected in the last iteration
             if effective_α > 0.0

                 # Projected variable
                 pos = last_list_r[r]

                 # Correction Assuming A = I, such that (d⋅e_r)(A*e_r) = d_r e_r
                 #L .= L .+ effective_α.*Extract_as_vector(last_d,pos)

                 # Correction Assuming that A*e_r is (∇f(x^k) - ∇f(x^{k-1})) / ||Δx|| [pos] 
                 L .= L .+ effective_α.*last_d[pos]*(y.*y[pos]/ys)  


             end # effective_α

         end # r
                
         # Now we can evaluate beta 
         β = dot(L,D)/dot(L,last_d)

         # Avoid a very unfortunate corner case
         if isnan(β) #|| β<0.0
             β = 0.0
         end

         # New search direction
         d .= -D .+ β*last_d

         # Assert that this is a descent direction
         if dot(d,D)>=0.0 
            #println("Reverting to Steepest...")
            d .= -D
         else
             #println("Using $β ...")
             nothing
         end

end



end # module
