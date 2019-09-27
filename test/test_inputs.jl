@testset "Test inputs" begin


    # Silly f(x) and df(x) - Not used since we are forcing errors
    function f(x)
        sum(x)
    end

    function df(x)
        one(x)
    end

    # Default inputs
    options = WallE.Init()
    

    #
    # First test - Length of x0,ci and cs must be the same.
    # we are also testing if we can call Solve without input 
    # parameters (Solve must use default parameters)
    #
    x0 = ones(10); ci = zeros(5);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs)

    x0 = ones(1); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs)

    x0 = ones(10); ci = zeros(10);  cs = 2*ones(5)
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs)    
 
    #
    # Second test - Check if x0 is inside the bounds
    #
    x0 = -1*ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs)

    x0 = 5*ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs)

    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    x0[5] = 10.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs)
   
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    x0[2] = -10.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs)
   

    #
    # Third test - nmax_iter > 0
    #
    options["NITER"]=-1
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["NITER"]=0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
   
    # Restore default parameters
    options = WallE.Init()

    #
    # Fourth test - Check if tol_norm is in (0,1)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["TOL_NORM"]=0.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["TOL_NORM"]=1.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["TOL_NORM"]=10.2
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    
    # Restore default parameters
    options = WallE.Init()

    #
    # New test - Check if armijo_c (c) is in (0,0.5)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["ARMIJO_C"]=-1.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["ARMIJO_C"]=0.5
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["ARMIJO_C"]=0.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["ARMIJO_C"]=2.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
        

    # Restore default parameters
    options = WallE.Init()

    #
    # Fifth test - Check if cut_factor (τ) is in (0,1)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["ARMIJO_TAU"]=-1.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["ARMIJO_TAU"]=0.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["ARMIJO_TAU"]=1.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["ARMIJO_TAU"]=2.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
           
    # Restore default parameters
    options = WallE.Init()

    #
    # Sixth test - Check if α_ini is > 0.0
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["LS_ALPHA_INI"]=-1.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["LS_ALPHA_INI"]=0.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    
    # Restore default parameters
    options = WallE.Init()

    #
    # Seventh test - Check if α_min is << 1.0 and > 0. At least smaller than α_ini
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["LS_ALPHA_MIN"]=-1.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["LS_ALPHA_MIN"]=0.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["LS_ALPHA_MIN"]=10.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["LS_ALPHA_INI"]=1.0
    options["LS_ALPHA_MIN"]=1.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)

    # Restore default parameters
    options = WallE.Init()

    #
    # New test - σ must be in [τ,1.0)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["LS_SIGMA"]=1.0
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
    options["LS_TAU"]=0.1
    options["LS_SIGMA"]=0.09
    @test_throws AssertionError WallE.Solve(f,df,x0,ci,cs,options)
  
end