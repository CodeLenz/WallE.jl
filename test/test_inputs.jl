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
 
end