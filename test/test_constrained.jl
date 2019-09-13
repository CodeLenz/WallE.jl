#
# Function with many side constraints 
#
@testset "Many side constraints" begin



    # First test
    #  min f(x,y) = (x*(y-3))^2 + 4*x
    #
    #   -2 <= x <= ∞
    #    0 <= y <= ∞
    #
    #  solution in (-2.0, 3.0)
    #
    println("\n\n############\n  Test 2.1\n############")
    
    function f(x) 
        return (x[1]*(x[2]-3))^2 + 4*x[1]
    end

    function df(x)
        return [2*x[1]*(x[2]-3)^2 + 4 ; 2*x[1]^2*(x[2]-3)]
    end

   
    # Restrições laterais
    ci = [-2.0 ; 0.0]
    cs = [Inf ; Inf]

    # Ponto inicial
    x0 = [2.0 ; 0.0]


    # Chama o otimizador

    x_opt, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000)

    # The test
    @test isapprox(x_opt,[-2.0 ; 3.0],atol=1E-6)

    println("\n","# Results #")
    show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [-2.0 ; 3.0]])

    x_opt_GC, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000,ENABLE_GC=true)

    # The test
    @test isapprox(x_opt,[-2.0 ; 3.0],atol=1E-6)
    @test isapprox(x_opt_GC,[-2.0 ; 3.0],atol=1E-6)


    # Second test
    println("\n\n\n############\n  Test 2.2\n############")

    function f(x)
        valor = 0.0
        for i=1:100
            valor += (x[i]-i)^2
        end
        return valor
    end

    function df(x)
        d = zeros(100)
        for i=1:100
            d[i] = 2*(x[i]-i)
        end
        return d
    end


    # Restrições laterais
    ci = 10*ones(100)
    cs = 50*ones(100)

    # Ponto inicial
    x0 = max.(ci,30*rand(100))


    # Chama o otimizador
    x_opt, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000,1E-8)
    x_opt_GC, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000,1E-8,ENABLE_GC=true)


    # The test
    @test isapprox(x_opt,[10*ones(10) ; 11:49 ; 50*ones(51)],rtol=1E-4)
    @test isapprox(x_opt_GC,[10*ones(10) ; 11:49 ; 50*ones(51)],rtol=1E-4)


    #println("\n","# Results #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [10*ones(10) ; 11:49 ; 50*ones(51)]])


    # Third test
    println("\n\n\n############\n  Test 2.3\n############")

    function f(x)
        valor = 0.0
        for i=1:100
            if isodd(i)
                valor += (x[i]-i)^1
            else
                valor += (x[i]-i)^2
            end
        end
        return valor
    end

    function df(x)
        d = zeros(100)
        for i=1:100
            if isodd(i)
                d[i] = 1.0
            else
                d[i] = 2*(x[i]-i)
            end
        end
        return d
    end


    # Restrições laterais
    ci = 10*ones(100)
    cs = 100*ones(100)

    # Ponto inicial
    x0 = max.(ci,60*rand(100))


    # Chama o otimizador
    x_opt, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000)
    x_opt_GC, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000,ENABLE_GC=true)

    # The test
    resp = 10*ones(100); [resp[i]=i for i in 12:2:100]
    @test isapprox(x_opt,resp,rtol=1E-4) 
    @test isapprox(x_opt_GC,resp,rtol=1E-4)

    #println("\n","# Results #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt resp])
    println("\n")


   println("\n\n############\n  Test 2.4\n############")
   # Função de Booth com restr
   # x1 <= 0.5 ; x2 <= 2.0
   # Dica de ponto inicial = (-5,-5)
   # minimo (0.5 , 2.0)

    function f(x) 
       (x[1]+2*x[2]-7)^2 + (2*x[1]+x[2]-5)^2 
    end

    function df(x)
        df1 = 2*(2*x[2]+x[1]-7)+4*(x[2]+2*x[1]-5)
        df2 = 4*(2*x[2]+x[1]-7)+2*(x[2]+2*x[1]-5)
        return [df1 ; df2]
    end

    # Ponto inicial
    x0 = [-5.0 ; -5.0]

    # Restrições laterais
    ci =  -Inf*ones(2)
    cs =  [0.5 ; 2.0]

    # Chama o otimizador
    x_opt, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000)
    x_opt_GC, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000,ENABLE_GC=true)
 

    # The test
    @test isapprox(x_opt,[0.5 ; 2.0],rtol=1E-2)
    @test isapprox(x_opt_GC,[0.5 ; 2.0 ],rtol=1E-2)

    #println("\n","# Resultado #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")





   println("\n\n############\n  Test 2.5\n############")
   # Rosenbrook
   # partida [0,3] solução sem bloqueio é (1,1)
   # x1 <= 0.8 e x2 >= 0.5 -> sol em (0.8 , 0.64)

    function f(x) 
        100*(x[2]-x[1]^2)^2+(x[1]-1)^2
    end

       
    function df(x)
        df1 = 2.0*(x[1]-1)-400*x[1]*(x[2]-x[1]^2)
        df2 = 200.0*(x[2]-x[1]^2)
        return [df1 ; df2]
    end

    # Ponto inicial
    x0 = [0.0 ; 3.0]

    # Restrições laterais
    ci = [-Inf ; 0.5]
    cs = [0.8 ; Inf] 

    # Chama o otimizador
    x_opt, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,10_000)
    x_opt_GC, f0, fopt, flag, _ = WallE.Wall_E2(f,df,x0,ci,cs,1000,ENABLE_GC=true)
 

    # The test
    @test isapprox(x_opt,[0.8 ; 0.64],rtol=1E-2)
    @test isapprox(x_opt_GC,[0.8 ; 0.64],rtol=1E-2)

    #println("\n","# Resultado #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")

 
end