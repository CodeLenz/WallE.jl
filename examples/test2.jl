using LinearAlgebra
#using WallE

#
# Minimiza uma função sem com restrições laterais
# Vamos gerar algo muito doido
#
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

function Test2()

    # Restrições laterais
    ci = 10*ones(100)
    cs = 50*ones(100)

    # Ponto inicial
    x0 = max.(ci,30*rand(100))

   

    # Chama o otimizador
    x_opt, flag, norma = WallE.Wall_E2(f,df,x0,ci,cs,true)

end
