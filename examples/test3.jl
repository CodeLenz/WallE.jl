using LinearAlgebra
using WallE


#
# Minimiza uma função  com restrições laterais
#  f(x,y) = x(y-3)^2 + 4*x
#
#  -2 <= x <= ∞
#   0 <= y <= ∞ 
#
# Solução em [0,0]
#
function f(x) 
   return x[1]*(x[2]-3)^2 + 4*x[1]
end

function df(x)
    [(x[2]-3)^2 + 4 ; 2*x[1]*(x[2]-3)]
end

function Test3()

    # Restrições laterais
    ci = [-2.0 ; 0]
    cs = [Inf ; Inf]

    # Ponto inicial
    x0 = [2.0 ; 0.0]

   
    # Chama o otimizador
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,true,true,500)

end
