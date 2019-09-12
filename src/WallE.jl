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
               x0::Array{Float64},
               ci::Array{Float64},
               cs::Array{Float64},
               nmax_iter::Int64=100,
               tol_norm::Float64=1E-6,
               flag_show::Bool=true,
               cut_factor::Float64=0.5,
               α_ini::Float64=10.0,
               α_min::Float64=1E-12; ENABLE_GC::Bool=false)

    
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
    last_free_x      = Int64[]
    last_d           = zeros(n)
    last_D           = zeros(n)
    last_α           = 0.0
    last_alpha_limit = Float64[]
    last_list_r      = Int64   

    # Norm (Gradient, free positions)
    norm_D = 0.0

    # Flag of convergence
    flag_conv = false

    # Number of effective iterations
    counter = 0
    
    # Number of iterations with GC
    cont_gc = 1

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


        # Evaluate deflection for the free variables.
        # It is not working :o(
        if ENABLE_GC && iter>1 && last_free_x == free_x && cont_gc <= n
          GC_projected!(d,D,last_D,last_d,free_x)
          cont_gc += 1
        else
            cont_gc = 0 
        end

        # Norm (to scale search direction)
        norm_d = norm(d)

        # Norm of the free positions
        norm_D = norm(D[free_x])

        # Store the norm (d)
        norms[iter] = norm_d
  
        # Scale
        d .= d./norm_d

        # Find limit alphas and constrained elements
        alpha_limit, list_r = Find_limit_alphas(x0,d,ci,cs)

        # Line search
        α, xn, fn, flag_sucess = Armijo_Projected(f,x0,
                                                fn,D,
                                                d,alpha_limit,
                                                list_r)
        # Store the step
        steps[iter] = α

        # Rollover Bethoven
        x0 .= xn
        last_free_x   = copy(free_x)
        last_d          .= d
        last_D          .= D
        last_α           = α
        last_alpha_limit = copy(alpha_limit)
        last_list_r      = copy(list_r)   

        # Blocked by below. They must be positive
        delta_m = D[I_blk_m]

        # Blocked by above. They must be negative
        delta_M = D[I_blk_M]

        # We need to fulfil all the first order conditions..
        if norm_D<=tol_norm*(1+abs(fn)) && (all(delta_m .>= 0.0)||isempty(delta_m)) &&
                                           (all(delta_M .<= 0.0)||isempty(delta_M))

            # Convergence assessed by first order condition. Set the flag and
            # skip the main loop
            #println("Convergence:: ",norm_D," ",tol_norm*(1+abs(fn))," ",
            #        (all(delta_m .>= 0.0)||isempty(delta_m))," ",
            #        (all(delta_M .<= 0.0)||isempty(delta_M)) )
            flag_conv = true
            break
        end # first order conditions
    
        # Fancy report for the mob :)
        ProgressMeter.next!(Prg; showvalues = [(:Norma,norm_D), (:Objetivo,fn), 
                          (:(Grad(Max)),maximum(D)), (:(Grad(Min)),minimum(D)),
                          (:Inferior,all(delta_m .>= -tol_norm)||isempty(delta_m)),
                          (:Superior,all(delta_M .<= tol_norm)||isempty(delta_M))],
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
function Extract_and_vector(v::Array{Float64},pos::Int64)
  vv = zero(v)
  vv[pos] = v[pos]
  return vv
end


#
# Given a point, a search direction and a step
# return the projected point and the list of
# effective blocks
#
function Project!(α::Float64,x0::Array{Float64},xn::Array{Float64},
                  d::Array{Float64},
                  alpha_limit::Array{Float64}, list_r::Array{Int64})


    # Length 
    n = size(x0,1)

    # Next point, without projections
    xn .= x0 + α*d

    # for each candidate constraint and for α > alpha_limit_r
    # apply the projection
    for r in LinearIndices(list_r)

        # Effective step
        α_eff = max(0.0, α - alpha_limit[r])

        if α_eff > 0.0

            xn .= xn .- (α_eff  .* Extract_and_vector(d,list_r[r]))

        end 

    end


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
         tol = 0.0 #1E-6
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
        Project!(α,x0,xn,d,alpha_limit,list_r)
 
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
# Evaluate the deflection for GC
#
function GC_projected!(d::Array{Float64},D::Array{Float64},last_D::Array{Float64},
                     last_d::Array{Float64},free_x::Array{Int64})

         #
         # Modified F & R proposed by Zhang, Zhou and Li
         #
         β = dot(D[free_x],D[free_x])/dot(last_D[free_x],last_D[free_x])

         # Difference in gradient
         y = D[free_x] .- last_D[free_x] 
        
         # 
         # Correction for the free terms
         #
         d[free_x] +=  -β*(dot(last_d[free_x],y)/dot(D[free_x],D[free_x]))*D[free_x]
                       +β*last_d[free_x]

         # Assert that this is a descent direction
         if dot(d,D)>=0.0 
          println("Reverting to SD...")
          d .= -D
         end

end



end # module
