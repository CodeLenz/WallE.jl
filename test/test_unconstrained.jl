#
# Function with no side constraints
#
@testset "No side constraints" begin


    # First test, optimal solution = (3.0 , 5.0)
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
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,100)
    x_opt_GC, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,100,ENABLE_GC=true)
 

    # The test
    @test isapprox(x_opt,[3.0 ; 5.0],rtol=1E-2)
    @test isapprox(x_opt_GC,[3.0 ; 5.0],rtol=1E-2)

    #println("\n","# Resultado #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")




   println("\n\n############\n  Test 1.2\n############")
   # Função de Booth
   # Dica de ponto inicial = (-5,-5)
   # minimo (1,3), valor 0.0

    function f(x) 
       (x[1]+2*x[2]-7)^2 + (2*x[1]+x[2]-5)^2 
    end

    function df(x)
        df1 = 2*(2*x[2]+x[1]-7)+4*(x[2]+2*x[1]-5)
        df2 = 4*(2*x[2]+x[1]-7)+2*(x[2]+2*x[1]-5)
        return [df1 ; df2]
    end

    # Ponto inicial
    x0 = -5*rand(2)

    # Restrições laterais
    ci = -Inf*ones(2)
    cs =  Inf*ones(2)

    # Chama o otimizador
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,100)
    x_opt_GC, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,100,ENABLE_GC=true)
 

    # The test
    @test isapprox(x_opt,[1.0 ; 3.0],rtol=1E-2)
    @test isapprox(x_opt_GC,[1.0 ; 3.0],rtol=1E-2)

    #println("\n","# Resultado #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")


end # testset