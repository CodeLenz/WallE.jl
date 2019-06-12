#
# Otimização por projeção no bounding-box
#
module WallE

using LinearAlgebra

export Wall_E, Wall_E2


function Wall_E(f::Function, df::Function, x0::Array{Float64},
                ci::Array{Float64}, cs::Array{Float64},
                flag_show::Bool=true,
                niter=0,  tol_norm=1E-6,
                passo_inicial=5.0, 
                fator_corte = 0.5,
                passo_minimo = 1E-10,
                limite_movel_inicial = 0.2,
                limite_movel_minimo = 0.001,
                fator_aumento_limite_movel = 1.1,
                fator_diminuicao_limite_movel = 0.7)


  # Numero de variáveis de projeto
  nx = length(x0)

  # Se o número de iterações não for informado, utilizamos um múltiplo
  # do número de variáveis
  if niter==0
     niter = 2*nx
  end
  
  if flag_show
     println("Wall_E::número máximo de iterações internas:: ",niter)
  end

  # Testa para ver se as dimensões das restrições laterais estão OK
  @assert length(ci)==nx "Steepest_Descent_Proj:: size(ci)!=size(x0)"
  @assert length(cs)==nx "Steepest_Descent_Proj:: size(cs)!=size(x0)"

  # Assert para verificar se  ci .<= cs
  for i=1:nx
      @assert ci[i]<=cs[i] "Steepest_Descent_Proj:: ci deve ser .<= cs "
  end

  # Não podemos aceitar qualquer x0 que viole as restrições laterais
  for i=1:nx
      @assert x0[i]>=ci[i] "Steepest_Descent_Proj:: x0[$i] deve ser maior do que ci[$i] $(x0[i]) $(ci[i])"
      @assert x0[i]<=cs[i] "Steepest_Descent_Proj:: x0[$i] deve ser menor do que cs[$i] $(x0[i]) $(cs[i])"
  end


  # Incializa variáveis internas

  # Copia o ponto de entrada para um vetor de estimativa de 
  # próximo ponto (a ser utilizado no line-search)
  xe = copy(x0)

  # Flag que indica se o otimizador saiu pela tolerância da norma2
  flag_conv = false

  # Flag que indica que a direção não é mais de minimização
  flag_minimizacao = true

  # Norma 2 (começa no valor máximo para Float64)
  norma = maxintfloat(Float64)
  norm_blocked = maxintfloat(Float64)

  # Norma do passo anterior
  norma_anterior = norma
 
  # Valor do passo (line-search)
  # Aqui usamos uma versão muito simples do backtracking, então
  # o passo inicial deve ser "elevado"
  passo = passo_inicial

  # Limites móveis. Todos iniciam no valor máximo
  limite_movel = limite_movel_inicial*ones(nx)

  # Vamos manter uma cópia dos valores das variáveis
  # de projeto nas últimas 3 iterações 
  x1 = zeros(nx) #Array{Float64}(undef,nx)
  x2 = zeros(nx) #Array{Float64}(undef,nx)

  # Limites móveis efetivos
  x_min = zeros(nx)
  x_max = zeros(nx)

  # Valor de referência. f0 será utilizado durante o line-search
  # e objetivo_inicial será utilizado para avaliarmos o quanto
  # de melhoria tivemos durante todo o processo. f1 será o valor
  # da estimativa no line-search
  f0 = f(x0)
  f1 = f0
  objetivo_inicial = f0

  # Conjunto de variáveis sendo bloqueadas no gradiente (block constrainted)
  Ilast  = Int64[]
  Iblock = Int64[]

  # Número de variáveis bloqueadas em ci ou em cs
  nblock_sup = 0
  nblock_inf = 0
 
  # Número de iterações efetivas
  contador = 0

  # Lista com todas as variáveis. Será utilizado para gerarmos 
  # uma lista de variáveis livres (complemento das bloqueadas)
  lvar = 1:nx

  # Direção de busca 
  d = zeros(nx)

  # Direção de busca anterior
  da = zeros(nx)

  # Vetor gradiente
  D = zeros(nx)


  # Número de iterações seguidas com gradientes conjugados
  number_fletcher = 1

  # Indica se usou GC em algum momento e quantas vezes
  usou_fletcher = 0

  ####################### LOOP PRINCIPAL #######################
  tempo = @elapsed for iter=1:niter

      # Loop do LS com bloqueio. Será verdadeiro se a tolerância
      # na norma do gradiente for satisfeita
      flag_conv_interna = false

      # Calcula a derivada de f, chamando df(x)
      D .= df(x0)

      # Atualiza o contador de iterações
      contador += 1

      # Cópia para análise de direção de minimização
      Doriginal = copy(D)

      # Calcula os limites móveis
      Moving_Limits!(limite_movel, x_min,x_max, x0, x1, x2, 
                     nx,iter,fator_aumento_limite_movel,
                     fator_diminuicao_limite_movel,limite_movel_minimo,
                     limite_movel_inicial,ci,cs)
       
      

      # Testa para ver se alguma variável de projeto que já 
      # está no box será levada para fora do box (por causa
      # do gradiente). Neste caso, levamos as variáveis para o 
      # box e zeramos as componentes do gradiente. Também
      # geramos uma lista de elementos bloqueados.
      Ilast = copy(Iblock)
      Iblock, nblock_inf, nblock_sup, norm_blocked = Select_Sets!(D,x0,x_min,x_max)

      # Calcula a norma atual das posições não bloqueadas do gradiente
      norma_anterior = norma
      norma = norm(D)

      # Se a tolerância da norma for satisfeita, setamos 
      # o flag_conv como verdadeiro e saimos do loop iter
      if norma<=tol_norm #|| norm_blocked <= tol_norm
         flag_conv = true
         break
      end

      # Cópia da direção de busca anterior
      da .= d 

      # Dependendo da situação, calculamos a direção de minimização
      # por Steepest Descent ou por Gradientes Conjugados (Fletcher and Reeves)
      # Steepest é selecionado se:
      #    1) Iter==1 (primeira iteração)
      #    2) Houve uma alteração na lista de elementos bloqueados
      #    3) numero de iterações por Gradientes Conjugados atingiu o limite
      #
      # Do contrário, utilizamos gradientes conjugados. 
      #
      method="Steepest"
      if iter==1 || Iblock != Ilast || mod(number_fletcher,nx)==0
         d .= -D/norma
         number_fletcher = 1
      else
         method="Fletcher"
         usou_fletcher += 1
         number_fletcher += 1
         beta = (norma/norma_anterior)^2
         d .= -D .+ beta*d
         d .= d/norm(d)
      end

       # Testa a direção de minimização. Sabemos que D⋅d <=0 
       # Se o produto interno não for negativo, podemos sair
       # do loop principal
       produto_busca = dot(Doriginal,d)
       if produto_busca > 0.0
           #print("\r Direção bloqueada não é mais minimizante $(number_fletcher)                                         ")
           flag_conv_interna = true
           flag_minimizacao = false
           break
       end

       ###################### LINE SEARCH #######################

       # Seta o flag de convergência para false (sem convergência)
       flag_conv_interna = false

       # Passo no começo da busca
       passo0 = passo

       
       # Testa passos até que o valor fique muito pequeno
       while passo > passo_minimo

           # Estimativa do próximo ponto
           xe .= x0 .+ passo*d

           # Projeta a estimativa no bounding box e gera uma lista
           # de posições bloqueadas. x1 é modificado no processo.
           # Esta lista será fundida com  Iblock. A rotina também 
           # indica se a direção projetada é minimizante.
           # Se a direção projetada for minimizante, então 
           # flag_direction será verdadeiro
           # 
           Iwall, flag_direction = Wall!(xe,x_min,x_max,Doriginal,x0)
         
           # Vamos fundir as listas de bloqueio
           Iblock = sort(unique([Iwall;Iblock])) 


           # Verifica o valor de f neste novo pto, mas só se a direção
           # for minimizante.
           if flag_direction
              f1 = f(xe)
           end
           
           # Condição de descida bem simples. Se melhorou,
           # então aceitamos o ponto. Estou fazendo isto devido
           # ao fato de estarmos projetando as variáveis a cada passo
           # o que muda efetivamente a direção de busca. Conforme já
           # provado (dissertação do Gustavo) estas direções ainda são
           # de minimização.
           if  f1<f0  && flag_direction 

             # Aceita o ponto e sai do loop de LS,
             # atualizando o valor de f0 (referência)
             # e do ponto atual
             x0 .= xe
             f0 = f1

             # Aumenta o passo ligeiramente, mas só se for menor
             # do que o passo inicial
             if passo0 < passo_inicial
                passo = passo0*2.0
             end

             # Setamos flag_conv_interna para verdadeiro 
             # e saimos do while (line-search)
             flag_conv_interna = true
             break

          else

             # Não minimizamos o valor da função. Diminui o passo
             passo = passo*fator_corte

          end #if f1<f0
 
       end #while

       #print(" Iteração interna $(iter) | norma $(norma) | metodo $(method)::$(number_fletcher) | nblock $(length(Iblock)) | passo $(passo) | $(produto_busca)                                \r")
       #flush(stdout)
       
       

     # Se chegamos aqui, então devemos testar pela
     # convergência da norma.
     if flag_conv
        # Temos uma convergência por norma do gradiente. Podemos sair
        break
     end

    # Se flag_conv_interna for falso, então o LS não teve sucesso
    if !flag_conv_interna 
       #if flag_show
       #    println("\nLS_Proj:: Não tivemos convergência interna - Passo muito pequeno $passo")
       #end
       flag_conv = false
       break
   end

  end #iter


  # Tivemos uma solução
  if flag_show
      println("\n********************************************************")
      println("Final do Steepest com projeção")
      println("Objetivo inicial   : ", objetivo_inicial)
      println("Objetivo final     : ", f0)
      if objetivo_inicial!=0.0 && f0!=0.0
         println("% de min.          : ", 100*(f0-objetivo_inicial)/objetivo_inicial )
      end
      println("Bloqueios          : ", nblock_inf," ",nblock_sup)
      println("Numero de iterações: ", contador , " de ",niter)
      println("Converg. por norma : ", flag_conv)
      println("Direção de min.    : ", flag_minimizacao)
      println("Passo final        : ", passo)
      println("Usou GC            : ", usou_fletcher)
      println("fator móvel mínimo : ", minimum(limite_movel))
      println("fator móvel máximo : ", maximum(limite_movel))
      println("Limite móvel mínimo: ", minimum(x_min))
      println("Limite móvel máximo: ", maximum(x_max))
      println("Normas (free/block): ", norma," ",norm_blocked)
      
      
      println("Tempo total [min]  : ", tempo/60.0)
      println("********************************************************")
  end

  # Returns the solution, convergence flag and last norm
  return x0, flag_conv, norma, norm_blocked

end


#
# Selects the set containing the elements that have been blocked
# Returns the set, and the number of elements being blocked (down and up)
# AND MODIFIES the gradient and the design variables
#
function Select_Sets!(D::Array{Float64},x::Array{Float64},
                      ci::Array{Float64},cs::Array{Float64})

      # Empty set
      Iblock = Int64[]; sizehint!(Iblock,length(x))

      # Norm of blocked directions. If we do not block any
      # direction, than we must return something really large
      norm_blocked = 0.0

      # Test for any box constraint
      nblock_sup = 0
      nblock_inf = 0
      @inbounds for i in LinearIndices(D)
          if D[i]>=0.0 &&  x[i]<=ci[i]
             x[i] = ci[i]
             norm_blocked += D[i]^2
             D[i] = 0.0
             nblock_inf += 1
             push!(Iblock,i)
          elseif  D[i]<=0.0 && x[i]>=cs[i]
             x[i] = cs[i]
             norm_blocked += D[i]^2
             D[i] = 0.0
             nblock_sup += 1
             push!(Iblock,i)
          end
      end

      # Avoid a null norm when there is no variable 
      # being blocked.
      if nblock_sup >0 || nblock_inf >0
         norm_blocked = sqrt(norm_blocked)
      else
         norm_blocked = maxintfloat(Float64)
      end
 
      return Iblock, nblock_inf, nblock_sup, norm_blocked
    end


    #
    # Applies the "wall" x .= max.(ci,min.(x,cs)) --> Modifies x
    #
    function Wall!(x::Array{Float64},ci::Array{Float64},cs::Array{Float64},
                   D::Array{Float64},x0::Array{Float64})

          # List of Blocked variables
          Iblock_x = Int64[]

          for i in LinearIndices(x)
            if x[i]<=ci[i] 
               x[i] = ci[i]
               push!(Iblock_x,i)
            elseif  x[i]>=cs[i] 
               x[i] = cs[i]
               push!(Iblock_x,i)
            end
          end  

          # Test if the projected direction is a minimizer 
          Δ = x-x0
          Δ /= norm(Δ)

          direction = dot(Δ,D)
          flag_direction = true
          if direction > 0.0
             flag_direction = false
             #println("Wall!::projected direction is not a minimizer direction $direction")
          end

        
          return Iblock_x, flag_direction

    end



#
# Atualiza limite_movel, x_min e x_max
#
function Moving_Limits!(limite_movel::Array{Float64}, x_min::Array{Float64},x_max::Array{Float64},
                        x0::Array{Float64}, x1::Array{Float64},x2::Array{Float64},
                        nx::Int64, iter::Int64, 
                        fator_aumento_limite_movel::Float64,
                        fator_diminuicao_limite_movel::Float64,limite_movel_minimo::Float64,
                        limite_movel_inicial::Float64,ci::Array{Float64},cs::Array{Float64} )


       # Se a iteração for maior do que 3, então podemos começar a 
       # ajustar os limites móveis pelo histórico de cada variável

      # Vetor unitário
      UM = ones(nx)

       if iter>3

          for i=1:nx

             # Variações
             Δ1 = x0[i]-x1[i] 
             Δ2 = x1[i]-x2[i]

             # Ajustes
             if Δ1*Δ2 > 0.0
                # Boas vibrações...podemos aumentar o lm
                limite_movel[i] *= fator_aumento_limite_movel

             elseif Δ1*Δ2 < 0.0
                # Bad vibes...
                limite_movel[i] *= fator_diminuicao_limite_movel

             end

             # Testa limites dos ajustes
             limite_movel[i] = max(limite_movel_minimo,
                               min(limite_movel_inicial,limite_movel[i]))

          end # for i  
       end # if iter > 3

       #
       # Vamos calcular os limites móveis, sempre respeitando
       # as restrições laterais. Nas primeiras 2 iterações
       # o limite móvel é o inicial (que é o máximo)
       #
       x_min .= max.(ci, (UM.-limite_movel).*x0 ) 
       x_max .= min.(cs, (UM.+limite_movel).*x0 )

       # Atualiza x1 e x2
       x2 .= x1
       x1 .= x0

end
 

##########################################################################
##########################################################################
##########################################################################
#
# ALTERAÇÕES DA DOCUMENTAÇÃO
#
    #
    # Applies the "wall" x .= max.(ci,min.(x,cs)) --> Modifies x
    # and indicates if any variable is projected
    #
    function Wall2!(x::Array{Float64},ci::Array{Float64},cs::Array{Float64})

          # List of Blocked variables
          Iblock_x = Int64[]

          for i in LinearIndices(x)
            if x[i]<=ci[i] 
               x[i] = ci[i]
               push!(Iblock_x,i)
            elseif  x[i]>=cs[i] 
               x[i] = cs[i]
               push!(Iblock_x,i)
            end
          end  

          return Iblock_x 
   end

   #
   # Armijo's Bactracking LS ove f(x)
   #
   function Modified_Armijo(x0::Array{Float64},f0::Float64,D::Array{Float64},
                            ci::Array{Float64},cs::Array{Float64},
                            f::Function)

      # Fixed parameters
      c = 0.1
      τ = 0.5

      # Normalize D if its not yet normalized
      D = D./norm(D)

      # First value of α
      α = 10.0

      # "Optimal" point and function value
      xn = copy(x0) 
      fn = f0

      # Blocked variables
      Iblock = Int64[] 

      # Counter
      iter = 0

      # The rigth hand side of the inequality is variable 
      # in this version. So we will skip the loop when
      # the search condition is true
      while true
 
        iter += 1

        # Evaluate the canditate point
        xn = x0 .- α*D

        # Projects the point into the boundary δS, modifying
        # xn 
        Iblock = Wall2!(xn,ci,cs)

        # The effective step is then 
        # (remember that we already projected xn into the box)
        Δx = xn .- x0
        
        # Such that m is given by 
        m = dot(D,Δx)

        # That should be negative
        @assert m < 0.0 "Armijo::not a descent direction"

        # And the condition is 
        fn = f(xn)
        #@show fn, f0 + c*m, α, length(Iblock)
        if  fn < (f0 + c*m)
           break 
        else
           α *= τ
        end

        # Check for a lower bound for α
        # TODO
        if α < 1E-8
          break
        end

        last_f = fn
        last_x = copy(xn)

      end # while

      # We should have a better point by now
      return xn, fn, Iblock
 
   end


  #
   # Crude LS ove f(x)
   #
   function Crude_LS(x0::Array{Float64},f0::Float64,D::Array{Float64},
                     ci::Array{Float64},cs::Array{Float64},
                     f::Function)

      # Fixed parameters (decrease of the step length)
      τ = 0.5

      # Normalize D if its not yet normalized
      D = D./norm(D)

      # First value of α
      α = 10.0

      # "Optimal" point and function value
      xn = copy(x0) 
      fn = f0

      # Blocked variables
      Iblock = Int64[] 

      # Counter
      iter = 0

      # Lets keep track of the minimization
      last_f = fn
      last_x = copy(xn)

      # Flag of improvement (or not)
      improved = true

      while true
 
        iter += 1

        # Evaluate the canditate point
        xn = x0 .- α*D

        # Projects the point into the boundary δS, modifying
        # xn 
        Iblock = Wall2!(xn,ci,cs)

        # The effective step is then 
        # (remember that we already projected xn into the box)
        Δx = xn .- x0
        
        # Such that m is given by 
        m = dot(D,Δx)

        # That should be negative
        @assert m < 0.0 "Crude_LS::not a descent direction"

        # And the condition is 
        fn = f(xn)
        #@show fn, last_f, iter
        if  fn < last_f
           α *= τ
        else
           # leave the while..Lets just take care of the 
           # first iteration
           if iter==1
             last_f = fn
             last_x .= xn
           end
           break
        end

        # Check for a lower bound for α
        if α < 1E-12
          improved = false
          break
        end

        last_f = fn
        last_x = copy(xn)

      end # while

      # We should have a better point by now
      return last_x, last_f, improved, Iblock
 
   end



function Wall_E2(f::Function, df::Function, x0::Array{Float64},
                 ci::Array{Float64}, cs::Array{Float64},
                 flag_show::Bool=true,
                 niter=2000,  tol_norm=1E-6,
                 passo_inicial=5.0, 
                 fator_corte = 0.5,
                 passo_minimo = 1E-10,
                 limite_movel_inicial = 0.2,
                 limite_movel_minimo = 0.001,
                 fator_aumento_limite_movel = 1.1,
                 fator_diminuicao_limite_movel = 0.7)


  # Numero de variáveis de projeto
  nx = length(x0)

  # Se o número de iterações não for informado, utilizamos um múltiplo
  # do número de variáveis
  if niter==0
     niter = 2*nx
  end
  
  if flag_show
     println("Wall_E2::número máximo de iterações internas:: ",niter)
  end

  # Lista com todas as variáveis. Será utilizado para gerarmos 
  # uma lista de variáveis livres (complemento das bloqueadas)
  lvar = 1:nx
  Iblock = Int64[]

  # Testa para ver se as dimensões das restrições laterais estão OK
  @assert length(ci)==nx "Wall_E2:: size(ci)!=size(x0)"
  @assert length(cs)==nx "Wall_E2:: size(cs)!=size(x0)"

  # Assert para verificar se  ci .<= cs
  for i=1:nx
      @assert ci[i]<=cs[i] "Wall_E2:: ci deve ser .<= cs "
  end

  # Não podemos aceitar qualquer x0 que viole as restrições laterais
  for i=1:nx
      @assert x0[i]>=ci[i] "Wall_E2:: x0[$i] deve ser maior do que ci[$i] $(x0[i]) $(ci[i])"
      @assert x0[i]<=cs[i] "Wall_E2:: x0[$i] deve ser menor do que cs[$i] $(x0[i]) $(cs[i])"
  end


  # Incializa variáveis internas

  # Copia o ponto de entrada para um vetor de estimativa de 
  # próximo ponto (a ser utilizado no line-search)
  xe = copy(x0)

  # Flag que indica se o otimizador saiu pela tolerância da norma2
  flag_conv = false


  # Norma 2 (começa no valor máximo para Float64)
  norma = maxintfloat(Float64)

 
  # Valor do passo (line-search)
  # Aqui usamos uma versão muito simples do backtracking, então
  # o passo inicial deve ser "elevado"
  passo = passo_inicial

  # Limites móveis. Todos iniciam no valor máximo
  limite_movel = limite_movel_inicial*ones(nx)

  # Vamos manter uma cópia dos valores das variáveis
  # de projeto nas últimas 3 iterações 
  x1 = zeros(nx) #Array{Float64}(undef,nx)
  x2 = zeros(nx) #Array{Float64}(undef,nx)

  # Limites móveis efetivos
  x_min = zeros(nx)
  x_max = zeros(nx)

  # Valor de referência. f0 será utilizado durante o line-search
  # e objetivo_inicial será utilizado para avaliarmos o quanto
  # de melhoria tivemos durante todo o processo. f1 será o valor
  # da estimativa no line-search
  f0 = f(x0)
  objetivo_inicial = f0
 
  # Número de iterações efetivas
  contador = 0

  # Vetor gradiente
  D = zeros(nx)

  # Será verdadeiro se a tolerância
  # na norma do gradiente for satisfeita
  flag_conv_interna = false


  ####################### LOOP PRINCIPAL #######################
  tempo = @elapsed for iter=1:niter

      # Calcula a derivada de f, chamando df(x)
      D .= df(x0)
 
      # Norma de D
      norma = norm(D) 

      # It iter > 1, than we can consider just the 
      # free (not blocked) variables to evaluate the norm
      #if iter>1
      #   free_x = filter(x->!(x in Iblock),lvar)
      #   norma = norm(D[free_x])
      #end

      # Se a tolerância da norma for satisfeita, setamos 
      # o flag_conv como verdadeiro e saimos do loop iter
      if norma<=tol_norm 
         flag_conv = true
         break
      end

      # Atualiza o contador de iterações
      contador += 1

      # Calcula os limites móveis
      Moving_Limits!(limite_movel, x_min,x_max, x0, x1, x2, 
                     nx,iter,fator_aumento_limite_movel,
                     fator_diminuicao_limite_movel,limite_movel_minimo,
                     limite_movel_inicial,ci,cs)
       

      # Armijo com próximo ponto, novo valor do objetivo e variáveis
      # bloqueadas
      x0, f0, improved, Iblock = Crude_LS(x0,f0,D,x_min,x_max,f)

      if !improved
        println("The solution cannot be improved during the line-search")
        break
      end

  end # for interno
       

  # Tivemos uma solução
  if flag_show
      println("\n********************************************************")
      println("Final do Steepest com projeção")
      println("Objetivo inicial   : ", objetivo_inicial)
      println("Objetivo final     : ", f0)
      if objetivo_inicial!=0.0 && f0!=0.0
         println("% de min.          : ", 100*(f0-objetivo_inicial)/objetivo_inicial )
      end
      #println("Bloqueios          : ", nblock_inf," ",nblock_sup)
      println("Numero de iterações: ", contador , " de ",niter)
      println("Converg. por norma : ", flag_conv)
      #println("Direção de min.    : ", flag_minimizacao)
      #println("Passo final        : ", passo)
      #println("Usou GC            : ", usou_fletcher)
      println("fator móvel mínimo : ", minimum(limite_movel))
      println("fator móvel máximo : ", maximum(limite_movel))
      println("Limite móvel mínimo: ", minimum(x_min))
      println("Limite móvel máximo: ", maximum(x_max))
      println("Norma              : ", norma)
      
      
      println("Tempo total [min]  : ", tempo/60.0)
      println("********************************************************")
  end

  # Returns the solution, convergence flag and last norm
  return x0, flag_conv, norma

end




end # module
