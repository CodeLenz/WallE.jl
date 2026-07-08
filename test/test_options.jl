@testset "Options and output contract" begin

    function quadratic_f(x)
        return (x[1]-2.0)^2 + 3.0*(x[2]+1.0)^2
    end

    function quadratic_df(x)
        return [2.0*(x[1]-2.0); 6.0*(x[2]+1.0)]
    end

    x0 = [-4.0; 5.0]
    ci = [-Inf; -Inf]
    cs = [Inf; Inf]

    for beta in ("SD","FR","PRP","HS_zero","HS_proj_e","HS_proj_sec")
        options = WallE.Init()
        options["SHOW"] = false
        options["NITER"] = 500
        options["BETA"] = beta

        output = WallE.Solve(quadratic_f,quadratic_df,x0,ci,cs,options)

        @test isapprox(output["RESULT"],[2.0; -1.0],rtol=1E-4,atol=1E-4)
        @test output["CONVERGED"]
        @test haskey(output,"NORM")
        @test haskey(output,"COUNTER_ITER")
        @test haskey(output,"NITER")
        @test haskey(output,"GC_FRACTION")
        @test haskey(output,"N_GC_USED")
        @test output["COUNTER_ITER"] == output["NITER"]
        @test 0.0 <= output["GC_FRACTION"] <= 1.0
        @test output["N_GC_USED"] <= output["COUNTER_ITER"]
    end

    options = WallE.Init()
    options["SHOW"] = false
    options["NITER"] = 500
    options["GATE"] = "non_conservative"

    output = WallE.Solve(quadratic_f,quadratic_df,x0,ci,cs,options)
    @test isapprox(output["RESULT"],[2.0; -1.0],rtol=1E-4,atol=1E-4)
    @test output["CONVERGED"]

    options = WallE.Init()
    options["SHOW"] = false
    delete!(options,"BETA")
    delete!(options,"GATE")
    delete!(options,"LS_TYPE")

    output = WallE.Solve(quadratic_f,quadratic_df,x0,ci,cs,options)
    @test isapprox(output["RESULT"],[2.0; -1.0],rtol=1E-4,atol=1E-4)
    @test output["CONVERGED"]

end

@testset "Projected descent gate" begin

    D = [1.0; 1.0]
    x = [0.0; 0.0]
    ci = [0.0; -Inf]
    cs = [Inf; Inf]

    @test WallE.Projected_descent(D,[-1.0; -1.0],x,ci,cs,1.0)
    @test !WallE.Projected_descent(D,[-1.0; 0.2],x,ci,cs,1.0)
    @test !WallE.Projected_descent(D,[0.0; 0.0],x,ci,cs,1.0)

end

@testset "Line search blocked direction" begin

    f(x) = x[1]^2
    df(x) = [2.0*x[1]]
    x0 = [0.0]
    D = [1.0]
    d = [-1.0]
    ci = [0.0]
    cs = [Inf]

    xn, fn, dfn, active_r, active_r_ci, active_r_cs, alpha, alpha_I, flag_success =
        WallE.Armijo_Projected!(f,df,x0,0.0,D,d,ci,cs,true,0.1,0.5,1.0,1E-12,0.9,false)

    @test !flag_success
    @test xn == x0
    @test dfn == D
    @test alpha <= 1E-12

end

@testset "Projection alpha split" begin

    xn, active_r, active_r_ci, active_r_cs, alpha_I =
        WallE.Project(0.25,[2.0; 0.0],[-40.0; 24.0],[-2.0; -Inf],[Inf; Inf])

    @test xn == [-2.0; 6.0]
    @test active_r == [1]
    @test active_r_ci == [1]
    @test isempty(active_r_cs)
    @test isapprox(alpha_I,[0.15],atol=eps())

    xn, active_r, active_r_ci, active_r_cs, alpha_I =
        WallE.Project(1.0,[0.0; 3.0],[2.0; -600.0],[-Inf; 0.5],[0.8; Inf])

    @test xn == [0.8; 0.5]
    @test active_r == [1,2]
    @test active_r_ci == [2]
    @test active_r_cs == [1]
    @test isapprox(alpha_I,[0.6,0.9958333333333333],rtol=0.0,atol=10eps())

end
