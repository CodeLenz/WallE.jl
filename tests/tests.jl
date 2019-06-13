using WallE
using Test


#
# Function with no side constraints
#
@testset "No side constraints" begin


    # First test, optimal solution = (3.0 , 5.0)

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
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,false,100)
 
    # The test
    @test isapprox(x_opt,[3.0 ; 5.0],rtol=1E-2)



end # testset



#
# Function with many side constraints 
#
@testset "Many side constraints" begin


    #  min f(x,y) = x(y-3)^2 + 4*x
    #
    #  -2 <= x <= ∞
    #   0 <= y <= ∞ 
    #
    # solution in (0.0, 0.0)
    #
    function f(x) 
       return x[1]*(x[2]-3)^2 + 4*x[1]
    end

    function df(x)
        [(x[2]-3)^2 + 4 ; 2*x[1]*(x[2]-3)]
    end

   
    # Restrições laterais
    ci = [-2.0 ; 0]
    cs = [Inf ; Inf]

    # Ponto inicial
    x0 = [2.0 ; 0.0]

       
    # Chama o otimizador
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,false,500)

    @test isapprox(x_opt,[0.0;0.0],atol=1E-6)



    # Second test

    function f(x) 
       valor = 0.0   
       for i=1:100
         valor +=  (x[i]-i)^2
        end
     return valor
    end

    function df(x)
        df = zeros(100)

        for i=1:100
            df[i] = 2*(x[i]-i)
        end
        return df
    end


    # Restrições laterais
    ci = 10*ones(100)
    cs = 50*ones(100)

    # Ponto inicial
    x0 = max.(ci,30*rand(100))


    # Chama o otimizador
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,false)

    @test isapprox(x_opt,[10*ones(10) ; 11:49 ; 50*ones(51)],rtol=1E-2)


end
