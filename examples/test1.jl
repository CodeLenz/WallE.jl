using LinearAlgebra
using WallE

#
# Minimiza uma função sem que restrições laterais
#
function f(x) 
   (x[1]-3)^2 + (x[2]-5)^2
end

function df(x)
	df1 = 2*(x[1]-3)
	df2 = 2*(x[2]-5)
	return [df1 ; df2]
end

function Test1()

    # Ponto inicial
    x0 = 10*rand(2)

    # Restrições laterais
    ci = -Inf*ones(2)
    cs = Inf*ones(2)

    # Chama o otimizador
    x_opt, flag, norma = Wall_E(f,df,x0,ci,cs,true,100)

end
