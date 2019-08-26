using WallE
using Test


#
# Function with no side constraints
#
@testset "No side constraints" begin

    # First test
    #
    # min (x-3)^2 + (y-5)^2
    #
    # optimal solution in (3.0 , 5.0)

    println("\n\n############\n  Test 1.1\n############")

    function f(x) 
       (x[1]-3)^2 + (x[2]-5)^2
    end

    function df(x)
    	df1 = 2*(x[1]-3)
    	df2 = 2*(x[2]-5)
    	return [df1 ; df2]
    end

    # Ponto inicial
    x0 = 10*rand(2)

    # Restrições laterais
    ci = -Inf*ones(2)
    cs =  Inf*ones(2)

    # Chama o otimizador
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,true,true,100)

    # The test
    @test isapprox(x_opt,[3.0 ; 5.0],rtol=1E-2)

    println("\n","# Resultado #")
    show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")


    # Second test - Rosenbrock function
    #
    # ∑_{i=1}^{n-1} [ 100( x_{i+1} - x_{i}^{2} )^{2} + ( 1 - x_{i} )^{2} ]
    #
    # Optimal solution f(1,...,1)=0.0

    println("\n\n############\n  Test 1.2\n############")

    function f(x)
        a = 0.0
        for i=1:length(x)-1
            a += 100*(x[i+1]-x[i]^2)^2 + (1-x[i])^2
        end
        return a
    end

    function df(x)
        dx = zero(x)
        for i=1:length(x)-1
            dx[i] += -400*x[i]*(x[i+1]-x[i]^2)-2*(1-x[i])
            dx[i+1] += 200*(x[i+1]-x[i]^2)
        end
        return dx
    end

    # Declaring problem size
    n = 1500

    # Initial point
    x0 = 10*ones(n)

    # Side constraints
    ci = -Inf*ones(n)
    cs =  Inf*ones(n)

    # Calling the optimizer
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,true,true,1000000)

    # The test
    @test isapprox(x_opt,ones(n),rtol=1E-5)

    println("\n","# Resultado #")
    show(IOContext(stdout, :compact => false, :limit => true), "text/plain", [x_opt ones(n)])
    println("\n")

end # testset
