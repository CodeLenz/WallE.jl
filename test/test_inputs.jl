@testset "Test inputs" begin


    # Silly f(x) and df(x) - Not used since we are forcing errors
    function f(x)
        sum(x)
    end

    function df(x)
        one(x)
    end

    #
    # First test - Length of x0,ci and cs must be the same
    #
    x0 = ones(10); ci = zeros(5);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs)

    x0 = ones(1); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs)

    x0 = ones(10); ci = zeros(10);  cs = 2*ones(5)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs)    
 
    #
    # Second test - Check if x0 is inside the bounds
    #
    x0 = -1*ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs)

    x0 = 5*ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs)

    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    x0[5] = 10.0
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs)
   
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    x0[2] = -10.0
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs)
   

    #
    # Third test - nmax_iter > 0
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,-1)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,0)
   

    #
    # Fourth test - Check if tol_norm is in (0,1)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,0.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,10.2)
    
      
    #
    # New test - Check if armijo_c (c) is in (0,0.5)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,-1.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.5)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,2.0)


    #
    # Fifth test - Check if cut_factor (τ) is in (0,1)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,-1.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,1.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,2.0)


    #
    # Sixth test - Check if α_ini is > 0.0
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,-1.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,0.0)
    

    #
    # Seventh test - Check if α_min is << 1.0 and > 0. At least smaller than α_ini
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,10.0,-1.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,10.0,0.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,10.0,10.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,10.0,100.0)

    #
    # New test - σ must be in [τ,1.0)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,10.0,1E-6,0.4)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,10.0,1E-6,1.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,10.0,1E-6,2.0)
    @test_throws AssertionError WallE.Wall_E2(f,df,x0,ci,cs,100,1E-8,true,0.1,0.5,10.0,1E-6,-1.0)
   

end