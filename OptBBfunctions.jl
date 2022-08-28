

 ##################################################################
 # Create Tree Structure for Branch-cut-and-price algorithm

   #type TreeNode
    struct TreeNode
      parent::Int
      children::Vector{Int}
      addedsolutions::Vector{Int}
      addedcenarios::Vector{Int}
      varxfixedsolutionzero::Vector{Vector{Int}}
      varxfixedsolutionone::Vector{Vector{Int}}
      vartfixedsolutionzero::Vector{Vector{Int}}
      vartfixedsolutionone::Vector{Vector{Int}}
      tsolution::Vector{Vector{Float64}}
      zsolution::Vector{Float64}
      varzfixedsolutionzero::Vector{Int}
      restr::Vector{ConstraintRef}
    end 

#
function solveSTODROPCEXACT(V,A,DESTC,PROBC,REWV,CAPV,PROBOD,REWOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    
    global firstpassnode1nonewscenario 
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vartsolutionfixedzero = []
    vartsolutionfixedone = []
    varzsolutionfixedzero = []
    restractive = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vartsolutionfixedzero = [vartsolutionfixedzero;tree[node].vartfixedsolutionzero]
      vartsolutionfixedone = [vartsolutionfixedone;tree[node].vartfixedsolutionone]
      varzsolutionfixedzero = [varzsolutionfixedzero;tree[node].varzfixedsolutionzero]
      restractive = [restractive;tree[node].restr]
      node = tree[node].parent
      node == 0 && break
    end
    #println("parent= ", tree[current].parent)
    #println("cenariosused =",cenariosused)
    #read(STDIN,Char)
    #println("routesused =",routesused)
    #read(STDIN,Char)
    #println("varxsolutionfixedzero =",varxsolutionfixedzero)
    #println("varxsolutionfixedone =",varxsolutionfixedone) 
    #println("vartsolutionfixedzero =",vartsolutionfixedzero)
    #println("vartsolutionfixedone =",vartsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    Gomorycutsrounds = 0
    #Now formulate problem (relaxed problem now) 
    z = Vector{Variable}()
    t = Vector{Vector{Variable}}()
    y = Vector{Vector{Vector{Vector{Variable}}}}()
    const1 = Vector{ConstraintRef}(length(cenariosused))
    #const3[w in 1:length(cenariosused),i in 1:length(DESTC)]
    const3= Vector{Vector{ConstraintRef}}()
    #const4[w in 1:length(cenariosused),i in 1:length(DESTC)]
    const4= Vector{Vector{ConstraintRef}}()
    #const5[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j]
    const5 = Vector{Vector{Vector{Vector{ConstraintRef}}}}()
    #const6[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC)]
    const6 = Vector{Vector{Vector{ConstraintRef}}}()
    #const7[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC)]
    const7 = Vector{Vector{Vector{ConstraintRef}}}()
    #const8[w in 1:length(cenariosused),i in 1:length(DESTC); scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0]
    const8 = Vector{Vector{ConstraintRef}}()
    #const9[w in 1:length(cenariosused),i in 1:length(DESTC); scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0]
    const9 = Vector{Vector{ConstraintRef}}()



    const11 = Vector{ConstraintRef}()
    #const12 = Vector{ConstraintRef}(length(scenario))

    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=3600))
    #Reserve first fixed indexes for continuous variables
    @variable(model,  s)
    @variable(model,  u[1:length(DESTC)])
    @variable(model,  v[1:length(DESTC)])
    for i in 1:length(DESTC)
      setupperbound(u[i],100000000)
      setupperbound(v[i],100000000)
    end
      
    
    for r in 1:length(routesused)
      if in(routesused[r],varzsolutionfixedzero)
        push!(z,@variable(model, lowerbound = 0, upperbound =0, basename="z[$r]"))
        JuMP.fix(z[end], 0)
      else
        push!(z,@variable(model, lowerbound = 0, basename="z[$r]"))
      end
    end

    #for i in 1:length(tree[current].zsolution)
    #     setvalue(z[i],tree[current].zsolution[i])  
    #end
  
    @variable(model,  x[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j]>=0)
    for i in 1:length(DESTC1)
      for j in 1:length(DESTC1)
        if i != j
          setupperbound(x[i,j],1)
        end
      end 
    end
   
    for w in 1:length(cenariosused)
      push!(t,[])
      for i in 1:length(DESTC)
        push!(t[end], @variable(model, lowerbound = 0, basename="t[$w][$i]"))
      end 
    end

    
    for w in 1:length(cenariosused)
     push!(y,[])
     for r in 1:length(routesused)
       push!(y[end], [] )
       for i in 1:length(DESTC1)
         push!(y[end][end], [] )
         for j in 1:length(DESTC1)
           if i == j
             push!(y[end][end][end], @variable(model, lowerbound = 0, upperbound =0, basename="y[$w][$r][$i][$j]") )
             JuMP.fix(y[end][end][end][j], 0)
           else
             #push!(y[end][end][end], @variable(model, lowerbound = 0, upperbound =1, basename="y[$w][$r][$i][$j]") )
             push!(y[end][end][end], @variable(model, lowerbound = 0, basename="y[$w][$r][$i][$j]") )
           end
         end
       end
     end 
    end

  

    for i in varxsolutionfixedzero
      setlowerbound(x[i[1],i[2]],0)
      setupperbound(x[i[1],i[2]],0) 
    end
    for i in varxsolutionfixedone
      setlowerbound(x[i[1],i[2]],1)
      setupperbound(x[i[1],i[2]],1)
    end

    for i in vartsolutionfixedzero
      JuMP.fix(t[i[1]][i[2]], 0)
      #setlowerbound(t[i[1]][i[2]],0)
      #setupperbound(t[i[1]][i[2]],0) 
    end
    for i in vartsolutionfixedone
      JuMP.fix(t[i[1]][i[2]], 1.0)
      #setlowerbound(t[i[1]][i[2]],1)
      #setupperbound(t[i[1]][i[2]],1)
    end

    
    
    @objective(model, Min, s + sum(PROBC[i]*u[i] + PROBOD[i]*v[i] for i in 1:length(DESTC))  )
    
    for w in 1:length(cenariosused)
      const1[w] = @constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] 
+ scenario[cenariosused[w]][length(DESTC)+i]*v[i]  for i in 1:length(DESTC))  -
sum( 1.0*PRICEOD[i]*t[w][i] for i in 1:length(DESTC)) - sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y[w][r][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= 0)
    end

    @constraint(model, const2[i in 1:length(DESTC)], sum(B1[routesused[r]][i]*z[r] for r in 1:length(routesused)) == 1)
      
    @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r]][i][j]*z[r] for r in 1:length(routesused))==0)

    for w in 1:length(cenariosused)
      push!(const3,[])
      for i in 1:length(DESTC)
        push!(const3[end],@constraint(model,  t[w][i] <= 1-scenario[cenariosused[w]][i]))
      end
    end

    for w in 1:length(cenariosused)
      push!(const4,[])
      for i in 1:length(DESTC)
        push!(const4[end],@constraint(model,  t[w][i] <= scenario[cenariosused[w]][length(DESTC)+i] ))
      end
    end
   
    for w in 1:length(cenariosused)
      push!(const5,[])
      for  r in 1:length(routesused)
        push!(const5[end],[])
        for i in 1:length(DESTC)
          push!(const5[end][end],[])
          for j in 1:length(DESTC)
            if i != j
              push!(const5[end][end][end],@constraint(model, y[w][r][i][j] - z[r] +t[w][i] +t[w][j] - sum( t[w][k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scenario[cenariosused[w]][i] +1 - scenario[cenariosused[w]][j] + B2[routesused[r]][i][j] + sum( scenario[cenariosused[w]][k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 ))
            else
              push!(const5[end][end][end],@constraint(model,0==0))
            end
          end
        end
      end
    end


    for w in 1:length(cenariosused)
      push!(const6,[])
      for  r in 1:length(routesused)
        push!(const6[end],[])
        for i in 1:length(DESTC)
          push!(const6[end][end],@constraint(model, y[w][r][length(DESTC1)][i]- z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const7,[])
      for  r in 1:length(routesused)
        push!(const7[end],[])
        for i in 1:length(DESTC)
          push!(const7[end][end],@constraint(model, y[w][r][i][length(DESTC1)] - z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const8,[])
      for i in 1:length(DESTC)
        if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
          push!(const8[end],@constraint(model, t[w][i] + sum(sum( y[w][r][i][j] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
        else
          push!(const8[end],@constraint(model,0==0))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const9,[])
      for i in 1:length(DESTC)
        if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
          push!(const9[end],@constraint(model, t[w][i] + sum(sum( y[w][r][j][i] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
        else
          push!(const9[end],@constraint(model,0==0))
        end
      end
    end
 
    const11 = [const11;restractive] #active cutting planes from parent nodes

    
   
    
    while true
       
      status = solve(model)
      result = getobjectivevalue(model)
      resultz = getvalue(z)  
      global TimeMain += toq()
      tic()
      println("Solved Node ",current," problem, result =", result, "; Time is ",TimeMain, " and Queue length is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #println(model)
      #read(STDIN,Char)

      stop = false
      mincost = false
      if status != :Optimal
        #prune   
        #pos = findfirst(x -> x == current, Queue)
        #deleteat!(Queue,pos)
        #delete!(Queue, current) 
        println("Infeasible...prune ", status)
        #read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolutiont
      global globalincumbentsolutionz
     
      #ysol = getvalue(y)
      xsol1 = getvalue(x)
      xsol = zeros(length(DESTC1),length(DESTC1))

      for i=1:length(DESTC1)
        for j=1:length(DESTC1)
          if i != j    
            if xsol1[i,j]<=0.00001 
              xsol[i,j] = 0
            elseif  xsol1[i,j] >=0.99999
              xsol[i,j]= 1
            else
              xsol[i,j]=xsol1[i,j]
            end
          end
        end
      end
      tsol1 = getvalue(t)
      tsol = zeros(length(cenariosused),length(DESTC))
      for w=1:length(cenariosused)
        for i=1:length(DESTC)
            
            if tsol1[w][i]<=0.00001 
              tsol[w,i] = 0
            elseif  tsol1[w][i] >=0.99999
              tsol[w,i]= 1
            else
              tsol[w,i]=tsol1[w][i]
            end
          
        end
      end
      
      #First execute pricing

      #Collect duals necessary for route reduced cost calculation 
      if current != 1 || firstpassnode1nonewscenario #pricing for node 1 only if no more scenario separation in first pass
      p = getdual(const1)
      #println("p= ",p)
      #readline() 
      alpha = getdual(const2)
      #println("alpha = ",alpha)
      #readline() 
      neta= getdual(const10)
      #println("neta= ",neta)
      #readline()
      gamma2 = getdual(const3)
      #println("gamma2= ",gamma2)
      #readline()
      gamma1 = getdual(const4)
      #println("gamma1= ",gamma1)
      #readline()
      beta = getdual(const5)
      #println("beta")
      #println("beta= ",beta)
      #readline()
      beta0i = getdual(const6)
      #println(beta0i)
      #readline()
      betai0 = getdual(const7)
      #println(betai0)
      #readline()
      gamma3 = getdual(const8)
      #println(gamma3)
      #readline()
      gamma4 = getdual(const9)
      #println(gamma4)
      #readline()
      #println(p)
      #println(alpha)
      #println(neta)
      #Calculate parameters K2 and K3 for dynamic prog algorithm
      K3 = 1.0*zeros(length(cenariosused),length(DESTC1),length(DESTC1))
      for w in 1:length(cenariosused)
        for i in 1:length(DESTC1)
          for j in 1:length(DESTC1)
            if i != j
              K3[w,i,j]=1.0*REWV*A[DESTC1[i],DESTC1[j]]*p[w]
              if i != length(DESTC1) && scenario[cenariosused[w]][i]==0 && scenario[cenariosused[w]][length(DESTC)+i]==1
                K3[w,i,j] -= gamma3[w][i]
              end
              if j != length(DESTC1) && scenario[cenariosused[w]][j]==0 && scenario[cenariosused[w]][length(DESTC)+j]==1
               K3[w,i,j] -= gamma4[w][j]
              end
            end
          end
        end
      end
   

      #Calculate K2[w,i]
      K2 = 1.0*zeros(length(cenariosused), length(DESTC)) 
      for w in 1:length(cenariosused)
        for k in 1:length(DESTC)
          K2[w,k] = 1.0*PRICEOD[k]*p[w] - gamma1[w][k] - gamma2[w][k] - sum(sum( beta[w][r][k][j] for j in 1:length(DESTC) if j !=k) for r in 1:length(routesused)) - sum(sum( beta[w][r][j][k] for j in 1:length(DESTC) if j !=k) for r in 1:length(routesused)) - sum( beta0i[w][r][k] for r in 1:length(routesused)) - sum( betai0[w][r][k] for r in 1:length(routesused))
        
          if scenario[cenariosused[w]][k]==0 && scenario[cenariosused[w]][length(DESTC)+k]==1
            K2[w,k] -= gamma3[w][k]
            K2[w,k] -= gamma4[w][k]  
          end
        end
      end
      
      #for w in 1:length(cenariosused)
      #  for k in 1:length(DESTC)
      #    k2iw= 1.0*PRICEOD[k]*p[w] - gamma1[w][k] - gamma2[w][k]
      #    if scenario[cenariosused[w]][k]==0 && scenario[cenariosused[w]][length(DESTC)+k]==1
      #      k2iw -= gamma3[w][k]
      #      k2iw -= gamma4[w][k]  
      #    end
      #    temp = 0
      #    for r in 1:length(routesused)
      #      for j in 1:length(DESTC)
      #        if k != j
      #          temp += beta[w][r][k][j]
      #          temp += beta[w][r][j][k]
      #        end
      #      end
      #      temp += betai0[w][r][k] + betai0[w][r][k]
      #      temp -= sum( 1.0*B2[routesused[r]][k][j]*betai0[w][r][j] for j in 1:length(DESTC) if j != k)
      #      temp -= sum( 1.0*B2[routesused[r]][j][k]*beta0i[w][r][j] for j in 1:length(DESTC) if j != k)
      #      temp -= sum( sum(1.0*B3[routesused[r]][i][j][k]*beta[w][r][i][j] for i in 1:length(DESTC) if i!=j && i !=k) for j in 1:length(DESTC) if j !=k)  
      #    end
      #    println("Dif reducedcost 106 original ",w," ,",k, ", ",temp-k2iw, ", ",getdual(t[w][k]))
          #if temp - k2iw >=0.1
          #   println("Erro dif reducedcost 106 original ",w," ,",k, ", ",temp," ,", k2iw, ", ",tsol[w,k]," ,",getdual(t[w][k]))
          #   #readline()
          #end
      #  end
      #end
     
      #readline()
      #Start dynamic programming algorithm

      #Scale and round up parameters for capacity constraint
      RDemand = round.(100*(1-PROBC))
      RCAPV=100*CAPV
      #println(RDemand, RCAPV)
      minRDemand =- minimum(RDemand)+1
      M = Array{Any}(floor(Int,RCAPV + minRDemand),length(DESTC))
      #Initialize cells
      for i in 1: floor(Int,(RCAPV + minRDemand))
        for j in 1:length(DESTC)
          if i == RDemand[j]+ minRDemand
            reducedcost = alpha[j]
            temp = 1.0*zeros(1:length(cenariosused))
            for w in 1:length(cenariosused)
              #add arcs to new j customer
              #First term of min 
              if scenario[cenariosused[w]][j] == 0 
                 temp[w] += K3[w,j,length(DESTC1)]
                 temp[w] += K3[w,length(DESTC1),j]
              end
              #Second term of min
              mink = Inf
              for k in 1:length(DESTC)
                temp2 = K2[w,k]
                temp2 += temp[w]
                if k == j  && scenario[cenariosused[w]][j] == 0
                  temp2 -= K3[w,k,length(DESTC1)]
                  temp2 -= K3[w,length(DESTC1),k]
                end
                if temp2 < mink
                  mink = temp2
                end
              end
              reducedcost -= min(temp[w], mink)
            end
            M[i,j]= [[j],reducedcost,alpha[j],0,temp]
          else
            M[i,j]= [[],-Inf,0,0,zeros(1:length(cenariosused))]
          end
        end
      end


      #fill up Matrix  label (route pathsofar, costsofar, costalphasofar, costnetasofar; costK3(w)*B8sofar)
      for d in 1: floor(Int,(RCAPV + minRDemand))
        for i in 1:length(DESTC)
          for j in 1:length(DESTC)
            if i != j  
              if length(M[d,i][1]) != 0 && findfirst(x -> x == j, M[d,i][1]) == 0 && M[d,i][1][end] == i && d - minRDemand + RDemand[j] <= RCAPV   
                #calculate reduced cost for extended visit from i to j with demand d
                reducedcost = M[d,i][3] - M[d,i][4]
                reducedcost += alpha[j]
                reducedcost -= neta[i,j]
                temp = 1.0*zeros(1:length(cenariosused))
                for w in 1:length(cenariosused)
                  temp[w] =  M[d,i][5][w]
                  #Eliminate previous links between last available customer towards depot and
                  #add arcs to new j customer 
                  if scenario[cenariosused[w]][j] == 0
                    for k in M[d,i][1]
                      pos = findfirst(x -> x == k, M[d,i][1])
                      if scenario[cenariosused[w]][k] == 0 && (k== M[d,i][1][end] || prod(scenario[cenariosused[w]][M[d,i][1][pos+1:end]] == 1) )
                        temp[w] -= K3[w,k,length(DESTC1)]
                        temp[w] += K3[w,k,j]
                        temp[w] += K3[w,j,length(DESTC1)]
                        break
                      end
                    end
                  end
                  #Second term of min
                  mink = Inf
                  for k in 1:length(DESTC)
                    temp2 = K2[w,k]+temp[w]
                    for k2 in vcat(M[d,i][1],[j])
                      if in(k,vcat(M[d,i][1],[j])) && k2 != k && scenario[cenariosused[w]][k2] == 0 && scenario[cenariosused[w]][k] == 0
                        pos = findfirst(x -> x == k, vcat(M[d,i][1],[j]))
                        pos2 = findfirst(x -> x == k2, vcat(M[d,i][1],[j]))
                        if pos < pos2 
                          if prod(scenario[cenariosused[w]][M[d,i][1][pos+1:pos2-1]]) == 1
                            temp2 -= K3[w,k,k2]
                          end 
                        elseif  pos > pos2
                          if prod(scenario[cenariosused[w]][M[d,i][1][pos2+1:pos-1]]) == 1
                            temp2 -= K3[w,k2,k]
                          end  
                        end
                      end
                    end
                    if temp2 < mink
                      mink = temp2
                    end
                  end 
                  reducedcost -= min(temp[w], mink)
                end
                if reducedcost > M[floor(Int,d+ RDemand[j]),j][2]
                 M[floor(Int,d+ RDemand[j]),j][1]= vcat(M[d,i][1],[j])
                 M[floor(Int,d+ RDemand[j]),j][2]= reducedcost
                 M[floor(Int,d+ RDemand[j]),j][3]= M[d,i][3]+ alpha[j]
                 M[floor(Int,d+ RDemand[j]),j][4]= M[d,i][4]+ neta[i,j]
                 M[floor(Int,d+ RDemand[j]),j][5]= temp
                end 
              end
            end
          end
        end
      end
      
      temproutes=[]
      temp = []
      maxreducedcost = -Inf
      for i in 1:length(DESTC)
        for d =floor(Int,(RCAPV + minRDemand)):-1:1
           #println("d= ",d,"; i= ",i," : ", M[d,i])           
           if M[d,i][2] != -Inf && M[d,i][2]>0.01
             ins=true
             for ind in 1:length(routesused)
               if length(M[d,i][1])==length(route[routesused[ind]]) && setdiff(M[d,i][1],route[routesused[ind]]) == []
                 ins=false
                 break
               end
             end
             for ind in 1:length(temproutes)
               if length(M[d,i][1])==length(temproutes[ind]) && setdiff(M[d,i][1],temproutes[ind]) == []
                 ins=false
                 break
               end
             end
 
             if  ins &&   M[d,i][2] > maxreducedcost
               maxreducedcost = M[d,i][2]
               temp =  M[d,i][1]
               #push!(temproutes, M[d,i][1])
               break
             end
           end
        end
        #readline()
      end
      if temp != []
        push!(temproutes,temp)
      end
      println(temproutes," ,",maxreducedcost)
      #routes with positive reduced  cost should be inserted
      if length(temproutes) != 0
        mincost = true
        for i = 1:length(temproutes)
          newroute = temproutes[i]
          #println("new routes= ", temproutes[i])
          push!(route,newroute)
          push!(routesused,size(route,1))  
          #Update Current node information on routes used
          push!(tree[current].addedsolutions,size(route,1))  
          global ROUTCOUNT += 1
      
          #Update B parameters
          push!(B1,1.0*zeros(1:length(DESTC)))
          push!(B2,[])
          for j in 1:length(DESTC1) 
            push!(B2[end],[])
            for j1 in 1:length(DESTC1)
              push!(B2[end][j],0 )
            end
          end
          push!(B3,[])
          for j in 1:length(DESTC) 
            push!(B3[end],[])
            for j1 in 1:length(DESTC)
              push!(B3[end][j],[])
              for j2 in 1:length(DESTC)
                push!(B3[end][j][j1],0 )
              end
            end
          end
          push!(B4,[])
          for j in 1:length(DESTC) 
            push!(B4[end],[])
            for j1 in 1:length(DESTC)
              push!(B4[end][j],[])
              for j2 in 1:length(DESTC)
                push!(B4[end][j][j1],0 )
              end
            end
          end
          push!(B5,[])
          for j in 1:length(DESTC) 
            push!(B5[end],[])
            for j1 in 1:length(DESTC)
              push!(B5[end][j],[])
              for j2 in 1:length(DESTC)
                push!(B5[end][j][j1],0 )
              end
            end
          end
          push!(B6,[])
          for j in 1:length(DESTC1) 
            push!(B6[end],[])
            for j1 in 1:length(DESTC1)
              push!(B6[end][j],0)
            end
          end
          push!(B7,zeros(1:length(DESTC)))
          for i in 1:length(scenario)
            push!(B8[i],[]) 
            for j1 in 1:length(DESTC1) 
              push!(B8[i][end],[])
              for j2 in 1:length(DESTC1)
                push!(B8[i][end][j1],0 )
              end
            end
          end
          r = size(route,1)
          for i in route[r]
            B7[r][i] = 1.0
            B1[r][i] += 1.0
          end
          B6[r][route[r][end]][length(DESTC1)] = 1.0
          B6[r][length(DESTC1)][route[r][1]] = 1.0
          for i in 1:length(DESTC)
            if findfirst(x -> x==i, route[r]) != 0
              B2[r][length(DESTC1)][i] = 1.0
              B2[r][i][length(DESTC1)] = 1.0
              if i != route[r][end]
                B6[r][i][route[r][findfirst(x -> x==i, route[r])+1]] = 1.0
                #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
              end
              #println("B1[ ",r,",",i," ]= ", B1[r,i])
              #read(STDIN,Char)
            end 
            for j in 1:length(DESTC)
              if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
                B2[r][i][j] = 1.0
                #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
                #read(STDIN,Char)
              end
              #println()
              for k in 1:length(DESTC)
                if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
                  B3[r][i][j][k] = 1.0
                  #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
                  #read(STDIN,Char)
                end
                #println()
                if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
                  B4[r][i][j][k] = 1.0
                  #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
                  #read(STDIN,Char)
                end
                #println()
                if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
                  B5[r][i][j][k] = 1.0
                  #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
                  #read(STDIN,Char)
                end
                #println()
              end
            end
          end
          for  w in 1:length(scenario),i in 1:length(DESTC),j in 1:length(DESTC)
            B8[w][r][i][j]  = (1- scenario[w][i])*(1 - scenario[w][j])*B2[r][i][j]*prod(scenario[w][k]*B3[r][i][j][k]+B4[r][i][j][k]+B5[r][i][j][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != j && k != i)
          end
          for  w in 1:length(scenario),i in 1:length(DESTC)
            B8[w][r][length(DESTC1)][i ]  = (1 - scenario[w][i])*prod(scenario[w][k]*B2[r][k][i]+B2[r][i][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != i)  
            B8[w][r][i][length(DESTC1)] = ( 1- scenario[w][i] )*prod(scenario[w][k]*B2[r][i][k]+B2[r][k][i] + 1-B7[r][k] for k in 1:length(DESTC) if k != i) 
          end 
          #Define new variables and constraints
          tes = [1.0*B1[r][1]]
          tes1 = [const2[1]]
          for i in 2:length(DESTC)
            push!(tes,1.0*B1[r][i])
            push!(tes1,const2[i])
          end
          for i in 1:length(DESTC1)
            for j in 1:length(DESTC1)
              if i != j
                push!(tes, -1.0*B6[r][i][j])
                push!(tes1, const10[i,j])
              end
            end
          end
          push!(z,@variable(model,  lowerbound = 0,  basename="z[$r]",objective = 0.0,inconstraints = tes1, coefficients = tes ))
      
   
          for w in 1:length(cenariosused)
            push!(y[w], [] )
            for i in 1:length(DESTC1)
              push!(y[w][end], [] )
              for j in 1:length(DESTC1)
                if i == j
                  push!(y[w][end][i], @variable(model, lowerbound = 0, upperbound =0, basename="y[$w][$r][$i][$j]") )
                  JuMP.fix(y[w][end][i][j], 0)
                else
                  tes1 = [const1[w]]
                  tes = [-1.0*REWV*A[DESTC1[i],DESTC1[j]]]
                  if i != length(DESTC1) && scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
                    push!(tes1, const8[w][i])
                    push!(tes, 1.0)
                  end
                  if j != length(DESTC1) && scenario[cenariosused[w]][length(DESTC)+j] == 1 && scenario[cenariosused[w]][j]==0
                    push!(tes1, const9[w][j])
                    push!(tes, 1.0)             
                  end
                  #push!(y[w][end][i], @variable(model, lowerbound = 0, upperbound =1, basename="y[$w][$r][$i][$j]",objective = 0.0, inconstraints = tes1, coefficients = tes) )
                  push!(y[w][end][i], @variable(model, lowerbound = 0, basename="y[$w][$r][$i][$j]",objective = 0.0, inconstraints = tes1, coefficients = tes) )
                end
              end
            end
          end
  
          r = length(routesused)
          for w in 1:length(cenariosused)
            push!(const5[w],[])
            for i in 1:length(DESTC)
              push!(const5[w][end],[])
              for j in 1:length(DESTC)
                if i != j
                  push!(const5[w][end][end],@constraint(model, y[w][r][i][j] - z[r] +t[w][i] +t[w][j] - sum( t[w][k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scenario[cenariosused[w]][i] +1 - scenario[cenariosused[w]][j] + B2[routesused[r]][i][j] + sum( scenario[cenariosused[w]][k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 ))
                else
                  push!(const5[w][end][end],@constraint(model,0==0))
                end
              end
            end
          end

          for w in 1:length(cenariosused)
            push!(const6[w],[])
            for i in 1:length(DESTC)
              push!(const6[w][end],@constraint(model, y[w][r][length(DESTC1)][i]- z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
            end
          end
 

          for w in 1:length(cenariosused)
            push!(const7[w],[])
            for i in 1:length(DESTC)
              push!(const7[w][end],@constraint(model, y[w][r][i][length(DESTC1)] - z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
            end
          end
        end
      end #end temproutes insert
      global TimePric += toq()
      tic()
      end  #end current != 1

      #println("mincost =", mincost)
      if !mincost
        #println("after mincost not true no column inserted")
        if  result > globalincumbentvalue + 0.0001
         #prune  since a lower bound is already worse then incumbent 
         #pos = findfirst(x -> x == current, Queue)
         #deleteat!(Queue,pos)
         #delete!(Queue, current)
         #println("Prune by Bound =", result, ", ", globalincumbentvalue)
         #read(STDIN,Char) 
         stop = true
         break
        end

        
        ssol = getvalue(s)
        usol= getvalue(u)
        vsol= getvalue(v)
        if current == 1 || (maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(tsol-floor.(tsol)) <= 0.00001)  #heuristic inserts all scenarios in node 1
          #new integer solution
          #println("Integer Solution Found")
          ##############################################
          #Formulate Scenario Separation Problem and Solve 
          #Get First Stage results
          for i in 1:length(resultz) #check for error and display
            if  resultz[i] > 0.0001 && resultz[i] < 0.99
              if current != 1
                println("SOLUCAO NAO INTEIRA!! ",i, ",",resultz[i])
                read(STDIN,Char)
              end
            end
          end


          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))
                          
          #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:2*length(DESTC)], Bin)
          @variable(model2, t2[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1)
          @variable(model2, obj1)
          @variable(model2, obj2) 
 
          #@objective(model2, Min, - ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC))  -  sum(scen[length(DESTC)+i]*vsol[i] for i in 1:length(DESTC)) + sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) + sum(PRICEOD[i]*t2[i] for i in 1:length(DESTC)))
          @objective(model2, Min, obj1 + obj2)
          @constraint(model2, obj1a, obj1 >=   ssol + sum(scen[i]*usol[i] for i in 1:length(DESTC))  +  sum(scen[length(DESTC)+i]*vsol[i] for i in 1:length(DESTC)) )
          @constraint(model2, obj2a, obj2 >=   sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) + sum(1.0*PRICEOD[i]*t2[i] for i in 1:length(DESTC)))
          @constraint(model2, obj3a, obj1 <= obj2)

          @constraint(model2, const2b1[i in 1:length(DESTC)], t2[i] <= scen[length(DESTC)+i])
          @constraint(model2, const2b2[i in 1:length(DESTC)], t2[i] <= 1-scen[i])
          @constraint(model2, const2b3[i in 1:length(DESTC)], scen[length(DESTC)+i]+scen[i] <=1)

          @constraint(model2, const2a1[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= resultz[r] )
          @constraint(model2, const2a2[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a3[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j] - t2[j])  )
          @constraint(model2, const2a4[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= B2[routesused[r]][i][j] )
          @constraint(model2,  const2a5[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC), k in 1:length(DESTC); i != j && k !=i && k !=j  && resultz[r] > 0.001], y2[r,i,j] <=  (scen[k]+ t2[k])*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1-B7[routesused[r]][k] )
          @constraint(model2,  const2a6[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,i,j] - resultz[r] +t2[i] +t2[j] -
sum( t2[k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scen[i] +1 - scen[j] + B2[routesused[r]][i][j] + sum( scen[k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 )

          @constraint(model2, const2a7[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= resultz[r] )
          @constraint(model2, const2a8[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a9[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= B7[routesused[r]][i]  )
          @constraint(model2,  const2a10[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,length(DESTC1),i] <=  (scen[j]+ t2[j])*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1-B7[routesused[r]][j]) 
          @constraint(model2,  const2a11[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] - resultz[r] +t2[i] - sum( t2[j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scen[i]  + B7[routesused[r]][i] + sum( scen[j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 )

          @constraint(model2, const2a12[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= resultz[r] )
          @constraint(model2, const2a13[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a14[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= B7[routesused[r]][i]  )
          @constraint(model2,  const2a15[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <=  (scen[j]+ t2[j])*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1-B7[routesused[r]][j] )
          @constraint(model2,  const2a16[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] - resultz[r] +t2[i] - sum( t2[j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scen[i]  + B7[routesused[r]][i] + sum( scen[j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 )

           
          #println("Will solve Separation Problem for integer solution")
          
          solve(model2)
          #println("solve model 2 ",getobjectivevalue(model2))
          #println(getvalue(obj1))
          #println(getvalue(obj2), ", ", result)
          #println(getvalue(obj1)-getvalue(obj2))
          #println(getvalue(scen))
          #read(STDIN,Char)
          scenarionew = round.(getvalue(scen))
          flag=false 
         
          #if getobjectivevalue(model2) <= -0.0001 && flag == false#-0.05  <= -0.0001
          if getvalue(obj1)-getvalue(obj2) <= -0.0001 #&& flag == false#-0.05  <= -0.0001
            println("ADD NEW SCENARIO ",getvalue(scen), " get= ", getobjectivevalue(model2))
            #read(STDIN,Char)
            #Add new scenario and continue
  
            push!(scenario,scenarionew)
            #Update Current node information on scenarios used
            push!(tree[current].addedcenarios,size(scenario,1))
            push!(cenariosused,size(scenario,1))
            #Create new variables
            w= length(cenariosused)
            push!(t,[])
            for i in 1:length(DESTC)
              push!(t[end], @variable(model, lowerbound = 0,  basename="t[$w][$i]"))
            end 
            push!(y,[])
            for r in 1:length(routesused)
              push!(y[end], [] )
              for i in 1:length(DESTC1)
                push!(y[end][r], [] )
                for j in 1:length(DESTC1)
                  if i == j
                    push!(y[end][r][i], @variable(model, lowerbound = 0, upperbound =0, basename="y[$w][$r][$i][$j]") )
                    JuMP.fix(y[end][r][i][j], 0)
                  else
                    push!(y[end][r][i], @variable(model, lowerbound = 0, basename="y[$w][$r][$i][$j]") )
                  end
                end
              end
            end 
  
            #Create new constraints
            push!(const1,@constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] 
+ scenario[cenariosused[w]][length(DESTC)+i]*v[i]  for i in 1:length(DESTC))  -
sum( 1.0*PRICEOD[i]*t[w][i] for i in 1:length(DESTC)) - sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y[w][r][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= 0) )

    
            push!(const3,[])
            for i in 1:length(DESTC)
              push!(const3[end],@constraint(model,  t[w][i] <= 1-scenario[cenariosused[w]][i]))
            end
        
            push!(const4,[])
            for i in 1:length(DESTC)
              push!(const4[end],@constraint(model,  t[w][i] <= scenario[cenariosused[w]][length(DESTC)+i] ))
            end
    
   
   
            push!(const5,[])
            for  r in 1:length(routesused)
              push!(const5[end],[])
              for i in 1:length(DESTC)
                push!(const5[end][r],[])
                for j in 1:length(DESTC)
                  if i != j
                    push!(const5[end][r][i],@constraint(model, y[w][r][i][j] - z[r] +t[w][i] +t[w][j] - sum( t[w][k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scenario[cenariosused[w]][i] +1 - scenario[cenariosused[w]][j] + B2[routesused[r]][i][j] + sum( scenario[cenariosused[w]][k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 ))
                  else
                    push!(const5[end][r][i],@constraint(model,0==0))
                  end
                end
              end
            end

            push!(const6,[])
            for  r in 1:length(routesused)
              push!(const6[end],[])
              for i in 1:length(DESTC)
                push!(const6[end][r],@constraint(model, y[w][r][length(DESTC1)][i]- z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
              end
            end
              
            push!(const7,[])
            for  r in 1:length(routesused)
              push!(const7[end],[])
              for i in 1:length(DESTC)
                push!(const7[end][r],@constraint(model, y[w][r][i][length(DESTC1)] - z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
              end
            end
  
            push!(const8,[])
            for i in 1:length(DESTC)
              if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
                push!(const8[end],@constraint(model, t[w][i] + sum(sum( y[w][r][i][j] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
              else
                  push!(const8[end],@constraint(model,0==0))
              end
            end


            push!(const9,[])
            for i in 1:length(DESTC)
              if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
                push!(const9[end],@constraint(model, t[w][i] + sum(sum( y[w][r][j][i] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
              else
                push!(const9[end],@constraint(model,0==0))
              end
            end
            
            #Calculate additional parameters B8
            global B8    #zeros(length(scenario),length(route),length(DESTC1),length(DESTC1))
            push!(B8,[])
            for j in 1:length(route)
              push!(B8[end],[]) 
              for j1 in 1:length(DESTC1) 
                push!(B8[end][j],[])
                for j2 in 1:length(DESTC1)
                  push!(B8[end][j][j1],0 )
                end
              end
            end
            
            
            for r in 1:length(route),i in 1:length(DESTC),j in 1:length(DESTC)
              B8[cenariosused[w]][r][i][j]  = (1- scenario[cenariosused[w]][i])*(1 - scenario[cenariosused[w]][j])*B2[r][i][j]*prod(scenario[cenariosused[w]][k]*B3[r][i][j][k]+B4[r][i][j][k]+B5[r][i][j][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != j && k != i)
            end
            for r in 1:length(route),i in 1:length(DESTC)
              B8[cenariosused[w]][r][length(DESTC1)][i ]  = (1 - scenario[cenariosused[w]][i])*prod(scenario[cenariosused[w]][k]*B2[r][k][i]+B2[r][i][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != i)  
              B8[cenariosused[w]][r][i][length(DESTC1)] = ( 1- scenario[cenariosused[w]][i] )*prod(scenario[cenariosused[w]][k]*B2[r][i][k]+B2[r][k][i] + 1-B7[r][k] for k in 1:length(DESTC) if k != i) 
            end 
              
            stop=false
          elseif  firstpassnode1nonewscenario && (maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(tsol-floor.(tsol)) <= 0.00001) # end getobj < +.05
            #println("No new scenario on integer solution")
            ##############################################
            #Verify if it can be new incubent
            if result < globalincumbentvalue
              #println("New Incumbent Found")
              globalincumbentvalue = result
              globalincumbentsolutiont = tsol 
              globalincumbentsolutionz = resultz
            end 
            #pos = findfirst(x -> x == current, Queue)
            #deleteat!(Queue,pos)
            #delete!(Queue, current)
            #println("Will break after integer solution")
            stop=true
            break
          elseif  current == 1
            if !firstpassnode1nonewscenario
              firstpassnode1nonewscenario = true
              stop = false
            else
              stop = true
              #Insert Gomory custs for variables  "x"
              if  Gomorycutsrounds == 0 && current == 1
                println("...will insert Gomory cut")
                Gomorycutsrounds += 1
                cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
                rowv2 = []

                for i in 1:length(DESTC1)
                  for j in 1:length(DESTC1)
                    if i != j
                      cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
                      #println("test x[",i," ,",j,"]") 
                      if cbasis[linearindex(x[i,j])] == :Basic && xsol[i,j]-floor(xsol[i,j]) > 0.00001  #is basic variable and is fractional, so insert gomory cut
                        #println("x[",i," ,",j,"] is elegible ", cbasis[linearindex(x[i,j])], ",",xsol[i,j],",",in([i;j],varxsolutionfixedzero)," ,",in([i;j],varxsolutionfixedone)) 
                        #read(STDIN,Char)
                        row = 1.0*zeros(length(cbasis))
                        #identify row where variable is basic
                        for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
                          #print("test x[",i," ,",j,"] in constraint ", k)
                          ci = model.internalModel.inner
                          ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
                          #println("(",linearindex(x[i,j]),",",row[linearindex(x[i,j])],",")
                          #if  row[linearindex(x[i,j])] != 0  read(STDIN,Char) end 
                          if row[linearindex(x[i,j])] == 1
                            #println("found x[",i," ,",j,"] in ", k) 
                            break   
                          end
                        end
                        if row == 1.0*zeros(length(cbasis))
                          println("ERROR GOMORY CUT X SECTION!")
                          read(STDIN,Char)
                        else #insert cut
                          #push!(const11,@constraint(model, sum(min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(x[i,j]) && row[k] != 0 ) >= 1) )
                          for k in 1:length(row)
                            if k != linearindex(x[i,j]) && row[k] != 0
                              if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
                                row[k] =  row[k]/(1- xsol[i,j]+floor(xsol[i,j]) )
                              elseif k <= 2*length(DESTC)+1  && row[k] > 0 
                                row[k] =  row[k]/(xsol[i,j]-floor(xsol[i,j]) )
                              elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= xsol[i,j]-floor(xsol[i,j]) 
                                row[k] =  (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) )
                              elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > xsol[i,j]-floor(xsol[i,j]) 
                                row[k] =  (1-row[k]+floor(row[k]))/(1-xsol[i,j]+floor(xsol[i,j]) )
                              else
                                println("ERROR! NOT COVERED GOMMORY X")
                                read(STDIN,Char)
                              end
                              #row[k] = min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )
                            end
                          end 
                          row[linearindex(x[i,j])]=0
                          #push!(row,xsol[i,j]-floor(xsol[i,j]) )
                          push!(rowv2,row)
                          #println("GOMORY CUT FOR VAR X!")
                          #read(STDIN,Char)
                        end 
                      end
                    end
                  end
                end

                       
                #Only at the end insert new constraints
                for v = 1:length(rowv2)
                  push!(const11,@constraint(model, sum( rowv2[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv2[v][k] != 0 ) >= 1 )  )
                  push!(tree[current].restr,const11[end])
                end
              end           #End insert gomory custs        
        
        
              #time for branching
              println("time for branching node 1")
              if maximum(xsol-floor.(xsol)) > 0.00001
                println("1a") 
                if length(vartsolutionfixedzero) != 0 || length(vartsolutionfixedone) != 0
                  #println("d ja estava fixed e veio xij frac")
                  #println(vartsolutionfixedzero)
                  #println(vartsolutionfixedone)
                  #read(STDIN,Char)
                end 
                println("2") 
                a = maximum(xsol-floor.(xsol))
                f=0
                g=0
                for i in 1:length(DESTC1),j in 1:length(DESTC1)
                  if i != j
                    if xsol[i,j]-floor.(xsol[i,j])== a
                      f = i
                      g  = j
                      break
                    end
                  end
                end
                println("3b") 
                #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
                b =[f;g] 
                push!(tree,TreeNode(current,[],[],[],[],[b],[],[],[],resultz,[],const11))
                push!(tree[current].children,length(tree))
                #push!(Queue,length(tree))
                Queue[length(tree)]=result

                println("4b") 
                #push!(tree,TreeNode(current,[],[],[],[b],[],[],[],tsol,resultz,[],const11))
                push!(tree,TreeNode(current,[],[],[],[b],[],[],[],[],resultz,[],const11)) 
                #current = length(tree)
                push!(tree[current].children,length(tree))
                ##push!(Queue,length(tree))
                Queue[length(tree)]=result

                ##pos = findfirst(x -> x == current, Queue)
                ##deleteat!(Queue,pos)
                ##delete!(Queue,current)

                ##current=Queue[end]
          
                ##Gomorycutsrounds = 0
                ##current=length(tree)

                #global NODECOUNT += 1
                #setlowerbound(x[f,g],0)
                #setupperbound(x[f,g],0)
                #push!(varxsolutionfixedzero,b) 
              elseif maximum(tsol-floor.(tsol)) > 0.00001
                println("5a") 
                a = maximum(tsol-floor.(tsol))
                f=0
                g=0
                for i in 1:length(cenariosused),j in 1:length(DESTC)
                  if tsol[i,j]-floor.(tsol[i,j])== a
                    f = i
                    g  = j
                    break
                  end
                end
                #println("6") 
                #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
                b =[f;g]
                #println("Sol d still frac, so split into two nodes using max frac", f)  
                #push!(tree,TreeNode(current,[],[],[],[],[],[],[b],tsol,resultz,[],const11))
                push!(tree,TreeNode(current,[],[],[],[],[],[],[b],[],resultz,[],const11))
                push!(tree[current].children,length(tree))
                #push!(Queue,length(tree))
                Queue[length(tree)]=result
                #push!(tree,TreeNode(current,[],[],[],[],[],[b],[],tsol,resultz,[],const11))
                push!(tree,TreeNode(current,[],[],[],[],[],[b],[],[],resultz,[],const11))
                push!(tree[current].children,length(tree))
                #push!(Queue,length(tree))
                Queue[length(tree)]=result
                #println("7") 
                #pos = findfirst(x -> x == current, Queue)
                #deleteat!(Queue,pos)
                #delete!(Queue,current)



                #current=Queue[end]

         
                #Gomorycutsrounds = 0
                #current=length(tree)

                #global NODECOUNT += 1
                #setlowerbound(d[f],0)
                #setupperbound(d[f],0)
                #push!(vardsolutionfixedzero,f) 
              end #
            end #end firstpassnode1nonewscenario
          end #end current = 1  
        else #Not a Integer solution
          stop = true
          #If root note, verify scenario to be inserted and continue solving master problem
             #but this version ot inserting scenarios for root node
          # Reduced Cost variable fixing for z
          zdual = getdual(z)
          for i in 1:length(zdual)
            if result + zdual[i] > globalincumbentvalue && !in(i,varzsolutionfixedzero)
              println("node ",current," z solution for route ",routesused[i], " will not improve")
              #read(STDIN,Char)
              #if resultz[i] >= 0.1 && resultz[i] <= 0.9
              #  println("erro z[", i,"]")
              #  read(STDIN,Char)
              #end
              #setlowerbound(z[i],0)
              #setupperbound(z[i],0)
              JuMP.fix(z[i], 0) 
              push!(varzsolutionfixedzero,routesused[i])
              push!(tree[current].varzfixedsolutionzero,routesused[i])
              #Should still eliminate constraints not needed here  
            end
          end
          #End Section for reduced Cost variable fixing

          #insert rounded capacity cuts if possible (deterministic)
          #println("verify rounded capacity vi")
          for num1 in 1:20    #try at most 20 times
            #create partition
            setpart = sample(collect(1:length(DESTC)), rand(5:length(DESTC)), replace = false)
            num2 = length(setpart)
            for num3 in 1:num2-3
              setpart = setdiff(setpart,[setpart[1]])
              #verify if partition violates x solution
              capacity = 2*ceil(prod((1-PROBOD[i]) for i in setpart)/CAPV)
              sumtemp=0
              for i in setpart
                sumtemp += sum(xsol[i,j] for j in 1:length(DESTC1))
              end
              if sumtemp < capacity
                #insert rounded capacity
                println("vi inserted")
                push!(const11, @constraint(model, sum(sum(x[i,j] for i in setpart ) for j=1:length(DESTC1) if i != j  ) >= capacity ) )
                push!(tree[current].restr,const11[end])
                break
              end
            end 
          end #end 20 times most

          #Insert Gomory custs for variables  "x"
          if  Gomorycutsrounds == 0 && current == 1
            println("...will insert Gomory cut")
            Gomorycutsrounds += 1
            cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
            rowv2 = []

            for i in 1:length(DESTC1)
              for j in 1:length(DESTC1)
                if i != j
                  cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
                  #println("test x[",i," ,",j,"]") 
                  if cbasis[linearindex(x[i,j])] == :Basic && xsol[i,j]-floor(xsol[i,j]) > 0.00001  #is basic variable and is fractional, so insert gomory cut
                    #println("x[",i," ,",j,"] is elegible ", cbasis[linearindex(x[i,j])], ",",xsol[i,j],",",in([i;j],varxsolutionfixedzero)," ,",in([i;j],varxsolutionfixedone)) 
                    #read(STDIN,Char)
                    row = 1.0*zeros(length(cbasis))
                    #identify row where variable is basic
                    for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
                      #print("test x[",i," ,",j,"] in constraint ", k)
                      ci = model.internalModel.inner
                      ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
                      #println("(",linearindex(x[i,j]),",",row[linearindex(x[i,j])],",")
                      #if  row[linearindex(x[i,j])] != 0  read(STDIN,Char) end 
                      if row[linearindex(x[i,j])] == 1
                        #println("found x[",i," ,",j,"] in ", k) 
                        break   
                      end
                    end
                    if row == 1.0*zeros(length(cbasis))
                      println("ERROR GOMORY CUT X SECTION!")
                      read(STDIN,Char)
                    else #insert cut
                      #push!(const11,@constraint(model, sum(min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(x[i,j]) && row[k] != 0 ) >= 1) )
                      for k in 1:length(row)
                        if k != linearindex(x[i,j]) && row[k] != 0
                          if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
                            row[k] =  row[k]/(1- xsol[i,j]+floor(xsol[i,j]) )
                          elseif k <= 2*length(DESTC)+1  && row[k] > 0 
                            row[k] =  row[k]/(xsol[i,j]-floor(xsol[i,j]) )
                          elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= xsol[i,j]-floor(xsol[i,j]) 
                            row[k] =  (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) )
                          elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > xsol[i,j]-floor(xsol[i,j]) 
                            row[k] =  (1-row[k]+floor(row[k]))/(1-xsol[i,j]+floor(xsol[i,j]) )
                          else
                           println("ERROR! NOT COVERED GOMMORY X")
                            read(STDIN,Char)
                          end
                          #row[k] = min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )
                        end
                      end 
                      row[linearindex(x[i,j])]=0
                      #push!(row,xsol[i,j]-floor(xsol[i,j]) )
                      push!(rowv2,row)
                      #println("GOMORY CUT FOR VAR X!")
                      #read(STDIN,Char)
                    end 
                  end
                end
              end
            end
                       
            #Only at the end insert new constraints
            for v = 1:length(rowv2)
              push!(const11,@constraint(model, sum( rowv2[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv2[v][k] != 0 ) >= 1 )  )
              push!(tree[current].restr,const11[end])
            end
          end           #End insert gomory custs        
        
        
          #time for branching
          println("time for branching")
          if maximum(xsol-floor.(xsol)) > 0.00001
            println("1") 
            if length(vartsolutionfixedzero) != 0 || length(vartsolutionfixedone) != 0
              #println("d ja estava fixed e veio xij frac")
              #println(vartsolutionfixedzero)
              #println(vartsolutionfixedone)
              #read(STDIN,Char)
            end 
            println("2") 
            a = maximum(xsol-floor.(xsol))
            f=0
            g=0
            for i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                if xsol[i,j]-floor.(xsol[i,j])== a
                  f = i
                  g  = j
                  break
                end
              end
            end
            println("3") 
            #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
            b =[f;g] 
            push!(tree,TreeNode(current,[],[],[],[],[b],[],[],[],resultz,[],const11))
            push!(tree[current].children,length(tree))
            #push!(Queue,length(tree))
            Queue[length(tree)]=result

             println("4") 
            #push!(tree,TreeNode(current,[],[],[],[b],[],[],[],tsol,resultz,[],const11))
            push!(tree,TreeNode(current,[],[],[],[b],[],[],[],[],resultz,[],const11)) 
            #current = length(tree)
            push!(tree[current].children,length(tree))
            ##push!(Queue,length(tree))
            Queue[length(tree)]=result

            ##pos = findfirst(x -> x == current, Queue)
            ##deleteat!(Queue,pos)
            ##delete!(Queue,current)

            ##current=Queue[end]
          
            ##Gomorycutsrounds = 0
            ##current=length(tree)

            #global NODECOUNT += 1
            #setlowerbound(x[f,g],0)
            #setupperbound(x[f,g],0)
            #push!(varxsolutionfixedzero,b) 
          elseif maximum(tsol-floor.(tsol)) > 0.00001
            println("5") 
            a = maximum(tsol-floor.(tsol))
            f=0
            g=0
            for i in 1:length(cenariosused),j in 1:length(DESTC)
              
                if tsol[i,j]-floor.(tsol[i,j])== a
                  f = i
                  g  = j
                  break
                end
              
            end
            println("6") 
            #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
            b =[f;g]
          #println("Sol d still frac, so split into two nodes using max frac", f)  
          #push!(tree,TreeNode(current,[],[],[],[],[],[],[b],tsol,resultz,[],const11))
          push!(tree,TreeNode(current,[],[],[],[],[],[],[b],[],resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result
          #push!(tree,TreeNode(current,[],[],[],[],[],[b],[],tsol,resultz,[],const11))
          push!(tree,TreeNode(current,[],[],[],[],[],[b],[],[],resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result
          println("7") 
          #pos = findfirst(x -> x == current, Queue)
          #deleteat!(Queue,pos)
          #delete!(Queue,current)



          #current=Queue[end]

         
          #Gomorycutsrounds = 0
          #current=length(tree)

          #global NODECOUNT += 1
          #setlowerbound(d[f],0)
          #setupperbound(d[f],0)
          #push!(vardsolutionfixedzero,f) 
          end #
        end  
        global TimeSce += toq()
        tic()
      end # !mincost 
      stop && break  
    end #end while true 'solve)
    return model, result 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate  Initial Routes, Scenarios, Prices and Incumbent solution
  ###############################################
  #Create vectors for all possible scenarios
  #global scenarioall = Vector{Vector{Int}}() 
  #fillscenario(scenarioall,zeros(length(DESTC)),length(DESTC),1)
  #println("All possible Scenario vectors created of size ", length(scenarioall) )
  println("Start run of Algorithm")
  #Create vectors for initial routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(DESTC)
    push!(route, [s])
  end
  println("Initial Route vectors created of size ", length(route))
  #################################################
 


  ###############################################
  #Create Vectors for initial scenarios 
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  push!(scenario, vcat(ones(1:length(DESTC)),zeros(1:length(DESTC)) )  )   #Scenario 1 no customer orders. 
  push!(scenario, vcat(zeros(1:length(DESTC)),ones(1:length(DESTC)) )  )   #Scenario 1 all customer orders with ODs.
  ################################################
  #println(scenario)
  #readline()
  println("Initial Scenarios vectors created of size ", length(scenario))

  ################################################
  #Bypass REWOD and define new prices to pay ODs
  global PRICEOD
  PRICEOD = Inf*1.0*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > 1.0*REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = 1.0*REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  ################################################


  println("Calculating first incumbent solution")
  ################################################
  #Initial Incumbent solution
  routetemp = [0]
  for num in DESTC
    minj = 0
    sumj = Inf
    for j in setdiff(collect(1:length(DESTC)), routetemp)
      sumtemp=0
      for i in routetemp
        if i == 0
          sumtemp += 1.0*REWV*A[1,DESTC[j]]*prod(PROBC[routetemp[2:end]])*(1-PROBC[j])
        else
          sumtemp += 1.0*REWV*A[DESTC[i],DESTC[j]]*prod(PROBC[routetemp[findfirst(x -> x==i, routetemp)+1:end]])*(1-PROBC[i])*(1-PROBC[j])
        end 
      end
      if sumtemp<sumj
        minj =j
        sumj = sumtemp
      end 
    end
    push!(routetemp,minj)
  end
  #println(routetemp)
  
  #Break ordering of routetemp by average vehicle capacity to create initial routes
  pos=2
  while true
    avgcap = 0
    routetemp2=[]
    while true
      avgcap += 1-PROBC[routetemp[pos]]
      if avgcap>CAPV  break end
      push!(routetemp2,routetemp[pos])
      #Calculate average capacity next
      pos +=1 
      if  pos == length(routetemp)+1 break end
    end
    #println(routetemp2)
    ins=true
    for ind in 1:length(route)
      if length(routetemp2)==length(route[ind]) && setdiff(routetemp2,route[ind]) == []
        ins=false
        break
      end
    end
    if ins
      push!(route,routetemp2) 
    end
    if pos == length(routetemp)+1 break end
  end
  #println(route)

  #To generate incumbent value for routes above, run first node with routes as integer values and insert all scenarios as row as needed.
  
  #readline()
  
  println("Calculating routes parameters ")
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = [] #zeros(length(route),length(DESTC))
  for i in 1:length(route)
    push!(B1,1.0*zeros(1:length(DESTC)))
  end
  global B2 = [] #zeros(length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(route)
    push!(B2,[])
    for j in 1:length(DESTC1) 
      push!(B2[i],[])
      for j1 in 1:length(DESTC1)
        push!(B2[i][j],0 )
      end
    end
  end

  global B3 = [] # zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B3,[])
    for j in 1:length(DESTC) 
      push!(B3[i],[])
      for j1 in 1:length(DESTC)
        push!(B3[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B3[i][j][j1],0 )
        end
      end
    end
  end
  global B4 = [] #zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B4,[])
    for j in 1:length(DESTC) 
      push!(B4[i],[])
      for j1 in 1:length(DESTC)
        push!(B4[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B4[i][j][j1],0 )
        end
      end
    end
  end
  global B5 = [] #zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B5,[])
    for j in 1:length(DESTC) 
      push!(B5[i],[])
      for j1 in 1:length(DESTC)
        push!(B5[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B5[i][j][j1],0 )
        end
      end
    end
  end
  global B6 = [] #zeros(length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(route)
    push!(B6,[])
    for j in 1:length(DESTC1) 
      push!(B6[i],[])
      for j1 in 1:length(DESTC1)
        push!(B6[i][j],0)
      end
    end
  end
  global B7 = [] #zeros(length(route),length(DESTC))
  for i in 1:length(route)
    push!(B7,zeros(1:length(DESTC)))
  end 
  global B8 = [] #zeros(length(scenario),length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(scenario)
    push!(B8,[])
    for j in 1:length(route)
      push!(B8[i],[]) 
      for j1 in 1:length(DESTC1) 
        push!(B8[i][j],[])
        for j2 in 1:length(DESTC1)
          push!(B8[i][j][j1],0 )
        end
      end
    end
  end
  
  for r in 1:length(route)
    for i in route[r]
      B7[r][i] = 1.0
      B1[r][i] += 1.0
    end
    B6[r][route[r][end]][length(DESTC1)] = 1.0
    B6[r][length(DESTC1)][route[r][1]] = 1.0
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        B2[r][length(DESTC1)][i] = 1.0
        B2[r][i][length(DESTC1)] = 1.0
        if i != route[r][end]
          B6[r][i][route[r][findfirst(x -> x==i, route[r])+1]] = 1.0
          #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
        end
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r][i][j] = 1.0
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r][i][j][k] = 1.0
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r][i][j][k] = 1.0
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r][i][j][k] = 1.0
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end
  for r in 1:length(route), w in 1:length(scenario),i in 1:length(DESTC),j in 1:length(DESTC)
    B8[w][r][i][j]  = (1- scenario[w][i])*(1 - scenario[w][j])*B2[r][i][j]*prod(scenario[w][k]*B3[r][i][j][k]+B4[r][i][j][k]+B5[r][i][j][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != j && k != i)
  end
  for r in 1:length(route), w in 1:length(scenario),i in 1:length(DESTC)
    B8[w][r][length(DESTC1)][i ]  = (1 - scenario[w][i])*prod(scenario[w][k]*B2[r][k][i]+B2[r][i][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != i)  
    B8[w][r][i][length(DESTC1)] = ( 1- scenario[w][i] )*prod(scenario[w][k]*B2[r][i][k]+B2[r][k][i] + 1-B7[r][k] for k in 1:length(DESTC) if k != i) 
  end 
  #####################################################################
  #read(STDIN,Char)
  global firstpassnode1nonewscenario = false
  println("Start Processing nodes")

  routesused = collect(1:length(route))
  #Initialize branch-and-bound tree and queue
  #Queue = Vector{Int}()
  Queue = PriorityQueue()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,collect(1:length(scenario)),[],[],[],[],[],[],[],[])) 
  #push!(Queue,1)
  Queue[1]= 0
  #Start processing queue in a Deep First mode
  global heursol
  global globalincumbentvalue = +Inf
  global globalincumbentsolutiont = []
  global globalincumbentsolutionz = []
    
  masterx = 0
  master = []
  #always process last element of queue
  tic()
  while length(Queue)>0
    #println("11") 
    #current=Queue[end]
    current= dequeue!(Queue)
    (master, masterx) =process(current)
    #println("12") 
    if current == 1
      global ROOTSOL=masterx
      global ROOTTIME= TimeMain +TimePric+TimeSce
    end
    #println("13") 
    global NODECOUNT += 1
    if length(Queue) > 0 && globalincumbentvalue != +Inf
      a,b=peek(Queue)
      if (globalincumbentvalue-b)/globalincumbentvalue < 0.02
        break 
      end
    end
    #println("14") 
    if TimeMain +TimePric+TimeSce >= 36000
      break
    end
  end
  println("Nodes processing has ended")  
  for i in 1:length(globalincumbentsolutionz)
    if globalincumbentsolutionz[i] > 0.0001 && globalincumbentsolutionz[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  #println("result= ",globalincumbentvalue,"dsol= ", globalincumbentsolutiond, "time = ", TimeMain +TimePric+TimeSce)
  global NROUTES = sum(globalincumbentsolutionz)
  #global PERHIGH = sum(globalincumbentsolutiont)/length(globalincumbentsolutiont)
  #print(route)
  return globalincumbentvalue,globalincumbentsolutiont
end #end solveSTODROPCEXACT function



 
#
function solveSTODROPCEXACTALLCOLUMNS(V,A,DESTC,PROBC,REWV,CAPV,PROBOD,REWOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    println("process")
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vartsolutionfixedzero = []
    vartsolutionfixedone = []
    varzsolutionfixedzero = []
    restractive = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vartsolutionfixedzero = [vartsolutionfixedzero;tree[node].vartfixedsolutionzero]
      vartsolutionfixedone = [vartsolutionfixedone;tree[node].vartfixedsolutionone]
      varzsolutionfixedzero = [varzsolutionfixedzero;tree[node].varzfixedsolutionzero]
      restractive = [restractive;tree[node].restr]
      node = tree[node].parent
      node == 0 && break
    end
    #println("parent= ", tree[current].parent)
    #println("cenariosused =",cenariosused)
    #read(STDIN,Char)
    #println("routesused =",routesused)
    #read(STDIN,Char)
    #println("varxsolutionfixedzero =",varxsolutionfixedzero)
    #println("varxsolutionfixedone =",varxsolutionfixedone) 
    #println("vartsolutionfixedzero =",vartsolutionfixedzero)
    #println("vartsolutionfixedone =",vartsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    Gomorycutsrounds = 0
    #Now formulate problem (relaxed problem now) 
    z = Vector{Variable}()
    t = Vector{Vector{Variable}}()
    y = Vector{Vector{Vector{Vector{Variable}}}}()
    const1 = Vector{ConstraintRef}(length(cenariosused))
    #const3[w in 1:length(cenariosused),i in 1:length(DESTC)]
    const3= Vector{Vector{ConstraintRef}}()
    #const4[w in 1:length(cenariosused),i in 1:length(DESTC)]
    const4= Vector{Vector{ConstraintRef}}()
    #const5[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j]
    const5 = Vector{Vector{Vector{Vector{ConstraintRef}}}}()
    #const6[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC)]
    const6 = Vector{Vector{Vector{ConstraintRef}}}()
    #const7[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC)]
    const7 = Vector{Vector{Vector{ConstraintRef}}}()
    #const8[w in 1:length(cenariosused),i in 1:length(DESTC); scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0]
    const8 = Vector{Vector{ConstraintRef}}()
    #const9[w in 1:length(cenariosused),i in 1:length(DESTC); scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0]
    const9 = Vector{Vector{ConstraintRef}}()



    const11 = Vector{ConstraintRef}()
    #const12 = Vector{ConstraintRef}(length(scenario))

    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=3600))
    #Reserve first fixed indexes for continuous variables
    @variable(model,  s)
    @variable(model,  u[1:length(DESTC)])
    @variable(model,  v[1:length(DESTC)])
    for i in 1:length(DESTC)
      setupperbound(u[i],100000000)
      setupperbound(v[i],100000000)
    end
      
    
    for r in 1:length(routesused)
      if in(routesused[r],varzsolutionfixedzero)
        push!(z,@variable(model, lowerbound = 0, upperbound =0, basename="z[$r]"))
        JuMP.fix(z[end], 0)
      else
        push!(z,@variable(model, lowerbound = 0, basename="z[$r]"))
      end
    end

    #for i in 1:length(tree[current].zsolution)
    #   setvalue(z[i],tree[current].zsolution[i])  
    #end
  
    @variable(model,  x[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j]>=0)
    for i in 1:length(DESTC1)
      for j in 1:length(DESTC1)
        if i != j
          setupperbound(x[i,j],1)
        end
      end 
    end
   
    for w in 1:length(cenariosused)
      push!(t,[])
      for i in 1:length(DESTC)
        push!(t[end], @variable(model, lowerbound = 0,  basename="t[$w][$i]"))
      end 
    end

    
    for w in 1:length(cenariosused)
     push!(y,[])
     for r in 1:length(routesused)
       push!(y[end], [] )
       for i in 1:length(DESTC1)
         push!(y[end][end], [] )
         for j in 1:length(DESTC1)
           if i == j
             push!(y[end][end][end], @variable(model, lowerbound = 0, upperbound =0, basename="y[$w][$r][$i][$j]") )
             JuMP.fix(y[end][end][end][j], 0)
           else
             push!(y[end][end][end], @variable(model, lowerbound = 0, basename="y[$w][$r][$i][$j]") )
           end
         end
       end
     end 
    end

  

    for i in varxsolutionfixedzero
      setlowerbound(x[i[1],i[2]],0)
      setupperbound(x[i[1],i[2]],0) 
    end
    for i in varxsolutionfixedone
      setlowerbound(x[i[1],i[2]],1)
      setupperbound(x[i[1],i[2]],1)
    end

    for i in vartsolutionfixedzero
      #setlowerbound(t[i[1]][i[2]],0)
      #setupperbound(t[i[1]][i[2]],0)
      JuMP.fix(t[i[1]][i[2]], 0) 
    end
    for i in vartsolutionfixedone
      JuMP.fix(t[i[1]][i[2]], 1.0) 
      #setlowerbound(t[i[1]][i[2]],1)
      #setupperbound(t[i[1]][i[2]],1)
    end

    
    
    @objective(model, Min, s + sum(PROBC[i]*u[i] + PROBOD[i]*v[i] for i in 1:length(DESTC))  )
    
    for w in 1:length(cenariosused)
      const1[w] = @constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] 
+ scenario[cenariosused[w]][length(DESTC)+i]*v[i]  for i in 1:length(DESTC))  -
sum( PRICEOD[i]*t[w][i] for i in 1:length(DESTC)) - sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y[w][r][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= 0)
    end

    @constraint(model, const2[i in 1:length(DESTC)], sum(B1[routesused[r]][i]*z[r] for r in 1:length(routesused)) == 1)
      
    @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r]][i][j]*z[r] for r in 1:length(routesused))==0)

    for w in 1:length(cenariosused)
      push!(const3,[])
      for i in 1:length(DESTC)
        push!(const3[end],@constraint(model,  t[w][i] <= 1-scenario[cenariosused[w]][i]))
      end
    end

    for w in 1:length(cenariosused)
      push!(const4,[])
      for i in 1:length(DESTC)
        push!(const4[end],@constraint(model,  t[w][i] <= scenario[cenariosused[w]][length(DESTC)+i] ))
      end
    end
   
    for w in 1:length(cenariosused)
      push!(const5,[])
      for  r in 1:length(routesused)
        push!(const5[end],[])
        for i in 1:length(DESTC)
          push!(const5[end][end],[])
          for j in 1:length(DESTC)
            if i != j
              push!(const5[end][end][end],@constraint(model, y[w][r][i][j] - z[r] +t[w][i] +t[w][j] - sum( t[w][k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scenario[cenariosused[w]][i] +1 - scenario[cenariosused[w]][j] + B2[routesused[r]][i][j] + sum( scenario[cenariosused[w]][k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 ))
            else
              push!(const5[end][end][end],@constraint(model,0==0))
            end
          end
        end
      end
    end


    for w in 1:length(cenariosused)
      push!(const6,[])
      for  r in 1:length(routesused)
        push!(const6[end],[])
        for i in 1:length(DESTC)
          push!(const6[end][end],@constraint(model, y[w][r][length(DESTC1)][i]- z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const7,[])
      for  r in 1:length(routesused)
        push!(const7[end],[])
        for i in 1:length(DESTC)
          push!(const7[end][end],@constraint(model, y[w][r][i][length(DESTC1)] - z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const8,[])
      for i in 1:length(DESTC)
        if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
          push!(const8[end],@constraint(model, t[w][i] + sum(sum( y[w][r][i][j] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
        else
          push!(const8[end],@constraint(model,0==0))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const9,[])
      for i in 1:length(DESTC)
        if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
          push!(const9[end],@constraint(model, t[w][i] + sum(sum( y[w][r][j][i] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
        else
          push!(const9[end],@constraint(model,0==0))
        end
      end
    end
 
    const11 = [const11;restractive] #active cutting planes from parent nodes

    
   

    while true
      println("start") 
      status = solve(model)
      result = getobjectivevalue(model)
      resultz = getvalue(z)
      println("a1")  
      global TimeMain += toq()
      tic()
      println("Solved Node ",current," problem, result =", result, "; Time is ",TimeMain, " and Queue length is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #println(model)
      #read(STDIN,Char)

      stop = false
      mincost = false
      if status != :Optimal
        #prune   
        #pos = findfirst(x -> x == current, Queue)
        #deleteat!(Queue,pos)
        #delete!(Queue, current) 
        println("Infeasible...prune ", status)
        #read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolutiont
      global globalincumbentsolutionz
     
      #ysol = getvalue(y)
      xsol1 = getvalue(x)
      xsol = zeros(length(DESTC1),length(DESTC1))

      for i=1:length(DESTC1)
        for j=1:length(DESTC1)
          if i != j    
            if xsol1[i,j]<=0.00001 
              xsol[i,j] = 0
            elseif  xsol1[i,j] >=0.99999
              xsol[i,j]= 1
            else
              xsol[i,j]=xsol1[i,j]
            end
          end
        end
      end
      tsol1 = getvalue(t)
      tsol = zeros(length(cenariosused),length(DESTC))
      for w=1:length(cenariosused)
        for i=1:length(DESTC)
            
            if tsol1[w][i]<=0.00001 
              tsol[w,i] = 0
            elseif  tsol1[w][i] >=0.99999
              tsol[w,i]= 1
            else
              tsol[w,i]=tsol1[w][i]
            end
          
        end
      end
      
      #No pricing 

      #println("mincost =", mincost)
      if true !mincost
        #println("after mincost not true no column inserted")
        if  result > globalincumbentvalue + 0.0001
         #prune  since a lower bound is already worse then incumbent 
         #pos = findfirst(x -> x == current, Queue)
         #deleteat!(Queue,pos)
         #delete!(Queue, current)
         #println("Prune by Bound =", result, ", ", globalincumbentvalue)
         #read(STDIN,Char) 
         stop = true
         break
        end

        
        ssol = getvalue(s)
        usol= getvalue(u)
        vsol= getvalue(v)
        if maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(tsol-floor.(tsol)) <= 0.00001
          #new integer solution
          #println("Integer Solution Found")
          ##############################################
          #Formulate Scenario Separation Problem and Solve 
          #Get First Stage results
          for i in 1:length(resultz) #check for error and display
            if resultz[i] > 0.0001 && resultz[i] < 0.99
              println("SOLUCAO NAO INTEIRA!! ",i, ",",resultz[i])
              read(STDIN,Char) 
            end
          end


          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))
                          
          #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:2*length(DESTC)], Bin)
          @variable(model2, t2[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1)
          @variable(model2, obj1)
          @variable(model2, obj2) 
 
          #@objective(model2, Min, - ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC))  -  sum(scen[length(DESTC)+i]*vsol[i] for i in 1:length(DESTC)) + sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) + sum(PRICEOD[i]*t2[i] for i in 1:length(DESTC)))
          @objective(model2, Min, obj1 + obj2)
          @constraint(model2, obj1a, obj1 >=   ssol + sum(scen[i]*usol[i] for i in 1:length(DESTC))  +  sum(scen[length(DESTC)+i]*vsol[i] for i in 1:length(DESTC)) )
          @constraint(model2, obj2a, obj2 >=   sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) + sum(PRICEOD[i]*t2[i] for i in 1:length(DESTC)))
          @constraint(model2, obj3a, obj1 <= obj2)

          @constraint(model2, const2b1[i in 1:length(DESTC)], t2[i] <= scen[length(DESTC)+i])
          @constraint(model2, const2b2[i in 1:length(DESTC)], t2[i] <= 1-scen[i])
          @constraint(model2, const2b3[i in 1:length(DESTC)], scen[length(DESTC)+i]+scen[i] <=1)

          @constraint(model2, const2a1[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= resultz[r] )
          @constraint(model2, const2a2[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a3[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j] - t2[j])  )
          @constraint(model2, const2a4[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= B2[routesused[r]][i][j] )
          @constraint(model2,  const2a5[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC), k in 1:length(DESTC); i != j && k !=i && k !=j  && resultz[r] > 0.001], y2[r,i,j] <=  (scen[k]+ t2[k])*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1-B7[routesused[r]][k] )
          @constraint(model2,  const2a6[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,i,j] - resultz[r] +t2[i] +t2[j] -
sum( t2[k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scen[i] +1 - scen[j] + B2[routesused[r]][i][j] + sum( scen[k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 )

          @constraint(model2, const2a7[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= resultz[r] )
          @constraint(model2, const2a8[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a9[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= B7[routesused[r]][i]  )
          @constraint(model2,  const2a10[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,length(DESTC1),i] <=  (scen[j]+ t2[j])*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1-B7[routesused[r]][j]) 
          @constraint(model2,  const2a11[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] - resultz[r] +t2[i] - sum( t2[j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scen[i]  + B7[routesused[r]][i] + sum( scen[j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 )

          @constraint(model2, const2a12[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= resultz[r] )
          @constraint(model2, const2a13[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a14[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= B7[routesused[r]][i]  )
          @constraint(model2,  const2a15[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <=  (scen[j]+ t2[j])*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1-B7[routesused[r]][j] )
          @constraint(model2,  const2a16[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] - resultz[r] +t2[i] - sum( t2[j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scen[i]  + B7[routesused[r]][i] + sum( scen[j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 )

           
          #println("Will solve Separation Problem for integer solution")
          
          solve(model2)
          #println("solve model 2 ",getobjectivevalue(model2))
          #println(getvalue(obj1))
          #println(getvalue(obj2), ", ", result)
          #println(getvalue(obj1)-getvalue(obj2))
          #println(getvalue(scen))
          #read(STDIN,Char)
          scenarionew = round.(getvalue(scen))
       
          #if getobjectivevalue(model2) <= -0.0001 && flag == false#-0.05  <= -0.0001
          if getvalue(obj1)-getvalue(obj2) <= -0.0001 #&& flag == false#-0.05  <= -0.0001
            println("ADD NEW SCENARIO ",getvalue(scen), " get= ", getobjectivevalue(model2))
            #read(STDIN,Char)
            #Add new scenario and continue
  
            push!(scenario,scenarionew)
            #Update Current node information on scenarios used
            push!(tree[current].addedcenarios,size(scenario,1))
            push!(cenariosused,size(scenario,1))
            #Create new variables
            w= length(cenariosused)
            push!(t,[])
            for i in 1:length(DESTC)
              push!(t[end], @variable(model, lowerbound = 0, basename="t[$w][$i]"))
            end 
            push!(y,[])
            for r in 1:length(routesused)
              push!(y[end], [] )
              for i in 1:length(DESTC1)
                push!(y[end][r], [] )
                for j in 1:length(DESTC1)
                  if i == j
                    push!(y[end][r][i], @variable(model, lowerbound = 0, upperbound =0, basename="y[$w][$r][$i][$j]") )
                    JuMP.fix(y[end][r][i][j], 0) 
                  else
                    push!(y[end][r][i], @variable(model, lowerbound = 0, basename="y[$w][$r][$i][$j]") )
                  end
                end
              end
            end 
  
            #Create new constraints
            push!(const1,@constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] 
+ scenario[cenariosused[w]][length(DESTC)+i]*v[i]  for i in 1:length(DESTC))  -
sum( PRICEOD[i]*t[w][i] for i in 1:length(DESTC)) - sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y[w][r][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= 0) )

    
            push!(const3,[])
            for i in 1:length(DESTC)
              push!(const3[end],@constraint(model,  t[w][i] <= 1-scenario[cenariosused[w]][i]))
            end
        
            push!(const4,[])
            for i in 1:length(DESTC)
              push!(const4[end],@constraint(model,  t[w][i] <= scenario[cenariosused[w]][length(DESTC)+i] ))
            end
    
   
   
            push!(const5,[])
            for  r in 1:length(routesused)
              push!(const5[end],[])
              for i in 1:length(DESTC)
                push!(const5[end][r],[])
                for j in 1:length(DESTC)
                  if i != j
                    push!(const5[end][r][i],@constraint(model, y[w][r][i][j] - z[r] +t[w][i] +t[w][j] - sum( t[w][k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scenario[cenariosused[w]][i] +1 - scenario[cenariosused[w]][j] + B2[routesused[r]][i][j] + sum( scenario[cenariosused[w]][k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 ))
                  else
                    push!(const5[end][r][i],@constraint(model,0==0))
                  end
                end
              end
            end

            push!(const6,[])
            for  r in 1:length(routesused)
              push!(const6[end],[])
              for i in 1:length(DESTC)
                push!(const6[end][r],@constraint(model, y[w][r][length(DESTC1)][i]- z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
              end
            end
              
            push!(const7,[])
            for  r in 1:length(routesused)
              push!(const7[end],[])
              for i in 1:length(DESTC)
                push!(const7[end][r],@constraint(model, y[w][r][i][length(DESTC1)] - z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
              end
            end
  
            push!(const8,[])
            for i in 1:length(DESTC)
              if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
                push!(const8[end],@constraint(model, t[w][i] + sum(sum( y[w][r][i][j] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
              else
                  push!(const8[end],@constraint(model,0==0))
              end
            end


            push!(const9,[])
            for i in 1:length(DESTC)
              if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
                push!(const9[end],@constraint(model, t[w][i] + sum(sum( y[w][r][j][i] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
              else
                push!(const9[end],@constraint(model,0==0))
              end
            end
            
            #Calculate additional parameters B8
            global B8    #zeros(length(scenario),length(route),length(DESTC1),length(DESTC1))
            push!(B8,[])
            for j in 1:length(route)
              push!(B8[end],[]) 
              for j1 in 1:length(DESTC1) 
                push!(B8[end][j],[])
                for j2 in 1:length(DESTC1)
                  push!(B8[end][j][j1],0 )
                end
              end
            end
            
            
            for r in 1:length(route),i in 1:length(DESTC),j in 1:length(DESTC)
              B8[cenariosused[w]][r][i][j]  = (1- scenario[cenariosused[w]][i])*(1 - scenario[cenariosused[w]][j])*B2[r][i][j]*prod(scenario[cenariosused[w]][k]*B3[r][i][j][k]+B4[r][i][j][k]+B5[r][i][j][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != j && k != i)
            end
            for r in 1:length(route),i in 1:length(DESTC)
              B8[cenariosused[w]][r][length(DESTC1)][i ]  = (1 - scenario[cenariosused[w]][i])*prod(scenario[cenariosused[w]][k]*B2[r][k][i]+B2[r][i][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != i)  
              B8[cenariosused[w]][r][i][length(DESTC1)] = ( 1- scenario[cenariosused[w]][i] )*prod(scenario[cenariosused[w]][k]*B2[r][i][k]+B2[r][k][i] + 1-B7[r][k] for k in 1:length(DESTC) if k != i) 
            end 
              
            stop=false
          else # end getobj < +.05
            #println("No new scenario on integer solution")
            ##############################################
            #Verify if it can be new incubent
            if result < globalincumbentvalue
              #println("New Incumbent Found")
              globalincumbentvalue = result
              globalincumbentsolutiont = tsol 
              globalincumbentsolutionz = resultz
            end 
            #pos = findfirst(x -> x == current, Queue)
            #deleteat!(Queue,pos)
            #delete!(Queue, current)
            #println("Will break after integer solution")
            stop=true
            break
          end
          global TimeSce += toq()
          tic()
        else #Not a Integer solution
          stop = true
          #If root note, verify scenario to be inserted and continue solving master problem
             #but this version ot inserting scenarios for root node
          # Reduced Cost variable fixing for z
          zdual = getdual(z)
          for i in 1:length(zdual)
            if result + zdual[i] > globalincumbentvalue && !in(i,varzsolutionfixedzero)
              println("node ",current," z solution for route ",routesused[i], " will not improve")
              #read(STDIN,Char)
              #if resultz[i] >= 0.1 && resultz[i] <= 0.9
              #  println("erro z[", i,"]")
              #  read(STDIN,Char)
              #end
              #setlowerbound(z[i],0)
              #setupperbound(z[i],0) 
              JuMP.fix(z[i], 0) 
              push!(varzsolutionfixedzero,routesused[i])
              push!(tree[current].varzfixedsolutionzero,routesused[i])
              #Should still eliminate constraints not needed here  
            end
          end
          #End Section for reduced Cost variable fixing

          #insert rounded capacity cuts if possible (deterministic)
          #println("verify rounded capacity vi")
          for num1 in 1:20    #try at most 20 times
            #create partition
            setpart = sample(collect(1:length(DESTC)), rand(5:length(DESTC)), replace = false)
            num2 = length(setpart)
            for num3 in 1:num2-3
              setpart = setdiff(setpart,[setpart[1]])
              #verify if partition violates x solution
              capacity = 2*ceil(prod((1-PROBOD[i]) for i in setpart)/CAPV)
              sumtemp=0
              for i in setpart
                sumtemp += sum(xsol[i,j] for j in 1:length(DESTC1))
              end
              if sumtemp < capacity
                #insert rounded capacity
                println("vi inserted")
                push!(const11, @constraint(model, sum(sum(x[i,j] for i in setpart ) for j=1:length(DESTC1) if i != j  ) >= capacity ) )
                push!(tree[current].restr,const11[end])
                break
              end
            end 
          end #end 20 times most

          #Insert Gomory custs for variables  "x"
          if  Gomorycutsrounds == 0 && current == 1
            println("...will insert Gomory cut")
            Gomorycutsrounds += 1
            cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
            rowv2 = []

            for i in 1:length(DESTC1)
              for j in 1:length(DESTC1)
                if i != j
                  cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
                  #println("test x[",i," ,",j,"]") 
                  if cbasis[linearindex(x[i,j])] == :Basic && xsol[i,j]-floor(xsol[i,j]) > 0.00001  #is basic variable and is fractional, so insert gomory cut
                    #println("x[",i," ,",j,"] is elegible ", cbasis[linearindex(x[i,j])], ",",xsol[i,j],",",in([i;j],varxsolutionfixedzero)," ,",in([i;j],varxsolutionfixedone)) 
                    #read(STDIN,Char)
                    row = 1.0*zeros(length(cbasis))
                    #identify row where variable is basic
                    for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
                      #print("test x[",i," ,",j,"] in constraint ", k)
                      ci = model.internalModel.inner
                      ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
                      #println("(",linearindex(x[i,j]),",",row[linearindex(x[i,j])],",")
                      #if  row[linearindex(x[i,j])] != 0  read(STDIN,Char) end 
                      if row[linearindex(x[i,j])] == 1
                        #println("found x[",i," ,",j,"] in ", k) 
                        break   
                      end
                    end
                    if row == 1.0*zeros(length(cbasis))
                      println("ERROR GOMORY CUT X SECTION!")
                      read(STDIN,Char)
                    else #insert cut
                      #push!(const11,@constraint(model, sum(min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(x[i,j]) && row[k] != 0 ) >= 1) )
                      for k in 1:length(row)
                        if k != linearindex(x[i,j]) && row[k] != 0
                          if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
                            row[k] =  row[k]/(1- xsol[i,j]+floor(xsol[i,j]) )
                          elseif k <= 2*length(DESTC)+1  && row[k] > 0 
                            row[k] =  row[k]/(xsol[i,j]-floor(xsol[i,j]) )
                          elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= xsol[i,j]-floor(xsol[i,j]) 
                            row[k] =  (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) )
                          elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > xsol[i,j]-floor(xsol[i,j]) 
                            row[k] =  (1-row[k]+floor(row[k]))/(1-xsol[i,j]+floor(xsol[i,j]) )
                          else
                           println("ERROR! NOT COVERED GOMMORY X")
                            read(STDIN,Char)
                          end
                          #row[k] = min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )
                        end
                      end 
                      row[linearindex(x[i,j])]=0
                      #push!(row,xsol[i,j]-floor(xsol[i,j]) )
                      push!(rowv2,row)
                      #println("GOMORY CUT FOR VAR X!")
                      #read(STDIN,Char)
                    end 
                  end
                end
              end
            end
                       
            #Only at the end insert new constraints
            for v = 1:length(rowv2)
              push!(const11,@constraint(model, sum( rowv2[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv2[v][k] != 0 ) >= 1 )  )
              push!(tree[current].restr,const11[end])
            end
          end           #End insert gomory custs        
        
        
          #time for branching
          println("time for branching")
          if maximum(xsol-floor.(xsol)) > 0.00001
            println("1") 
            if length(vartsolutionfixedzero) != 0 || length(vartsolutionfixedone) != 0
              #println("d ja estava fixed e veio xij frac")
              #println(vartsolutionfixedzero)
              #println(vartsolutionfixedone)
              #read(STDIN,Char)
            end 
            println("2") 
            a = maximum(xsol-floor.(xsol))
            f=0
            g=0
            for i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                if xsol[i,j]-floor.(xsol[i,j])== a
                  f = i
                  g  = j
                  break
                end
              end
            end
            println("3") 
            #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
            b =[f;g] 
            push!(tree,TreeNode(current,[],[],[],[],[b],[],[],[],resultz,[],const11))
            push!(tree[current].children,length(tree))
            #push!(Queue,length(tree))
            Queue[length(tree)]=result

             println("4") 
            #push!(tree,TreeNode(current,[],[],[],[b],[],[],[],tsol,resultz,[],const11))
            push!(tree,TreeNode(current,[],[],[],[b],[],[],[],[],resultz,[],const11)) 
            #current = length(tree)
            push!(tree[current].children,length(tree))
            ##push!(Queue,length(tree))
            Queue[length(tree)]=result

            ##pos = findfirst(x -> x == current, Queue)
            ##deleteat!(Queue,pos)
            ##delete!(Queue,current)

            ##current=Queue[end]
          
            ##Gomorycutsrounds = 0
            ##current=length(tree)

            #global NODECOUNT += 1
            #setlowerbound(x[f,g],0)
            #setupperbound(x[f,g],0)
            #push!(varxsolutionfixedzero,b) 
          elseif maximum(tsol-floor.(tsol)) > 0.00001
            println("5") 
            a = maximum(tsol-floor.(tsol))
            f=0
            g=0
            for i in 1:length(cenariosused),j in 1:length(DESTC)
              
                if tsol[i,j]-floor.(tsol[i,j])== a
                  f = i
                  g  = j
                  break
                end
              
            end
            println("6") 
            #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
            b =[f;g]
          #println("Sol d still frac, so split into two nodes using max frac", f)  
          #push!(tree,TreeNode(current,[],[],[],[],[],[],[b],tsol,resultz,[],const11))
          push!(tree,TreeNode(current,[],[],[],[],[],[],[b],[],resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result
          #push!(tree,TreeNode(current,[],[],[],[],[],[b],[],tsol,resultz,[],const11))
          push!(tree,TreeNode(current,[],[],[],[],[],[b],[],[],resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result
          println("7") 
          #pos = findfirst(x -> x == current, Queue)
          #deleteat!(Queue,pos)
          #delete!(Queue,current)



          #current=Queue[end]

         
          #Gomorycutsrounds = 0
          #current=length(tree)

          #global NODECOUNT += 1
          #setlowerbound(d[f],0)
          #setupperbound(d[f],0)
          #push!(vardsolutionfixedzero,f) 
          end #
        end
        println("8")  
        global TimeSce += toq()
        tic()
      end # !mincost
      println("9") 
      stop && break
      println("10")  
    end #end while true 'solve)
    return model, result 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate  Initial Routes, Scenarios, Prices and Incumbent solution
  ###############################################
  #Create vectors for all possible scenarios
  #global scenarioall = Vector{Vector{Int}}() 
  #fillscenario(scenarioall,zeros(length(DESTC)),length(DESTC),1)
  #println("All possible Scenario vectors created of size ", length(scenarioall) )
  
  


  ###############################################
  #Create Vectors for initial scenarios 
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  push!(scenario, vcat(ones(1:length(DESTC)),zeros(1:length(DESTC)) )  )   #Scenario 1 no customer orders. 
  push!(scenario, vcat(zeros(1:length(DESTC)),ones(1:length(DESTC)) )  )   #Scenario 1 all customer orders with ODs.
  ################################################
  #println(scenario)
  #readline()
  

  ################################################
  #Bypass REWOD and define new prices to pay ODs
  global PRICEOD
  PRICEOD = Inf*1.0*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > 1.0*REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = 1.0*REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  ################################################

 
  
  #readline()
  ################################################
  #Create vectors for all possible scenarios
  global scenarioall = Vector{Vector{Int}}() 
  fillscenario(scenarioall,zeros(length(DESTC)),length(DESTC),1)
  println("All possible Scenario vectors created of size ", length(scenarioall) )
  #Create vectors for all possible routings
  global routetemp = Vector{Vector{Int}}() 
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenarioall)
    if  true #sum((1-PROBC[i]) for i in 1:length(DESTC)  ) <= CAPV
      groupe = []
      for c in 1:length(scenarioall[s])
        if scenarioall[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),routetemp)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  #now filter only capacitated routes
  println("Route temp vectors created of size ", length(routetemp))
  ################################################# 
  RDemand = round.(100*(1-PROBC))
  RCAPV=100*CAPV
  for i in 1:length(routetemp)
    if sum(RDemand[routetemp[i][j]] for j in 1:length(routetemp[i])  ) <= RCAPV
       push!(route,routetemp[i])
    end
  end
  println("Route vectors created of size ", length(route))
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = [] #zeros(length(route),length(DESTC))
  for i in 1:length(route)
    push!(B1,1.0*zeros(1:length(DESTC)))
  end
  global B2 = [] #zeros(length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(route)
    push!(B2,[])
    for j in 1:length(DESTC1) 
      push!(B2[i],[])
      for j1 in 1:length(DESTC1)
        push!(B2[i][j],0 )
      end
    end
  end

  global B3 = [] # zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B3,[])
    for j in 1:length(DESTC) 
      push!(B3[i],[])
      for j1 in 1:length(DESTC)
        push!(B3[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B3[i][j][j1],0 )
        end
      end
    end
  end
  global B4 = [] #zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B4,[])
    for j in 1:length(DESTC) 
      push!(B4[i],[])
      for j1 in 1:length(DESTC)
        push!(B4[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B4[i][j][j1],0 )
        end
      end
    end
  end
  global B5 = [] #zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B5,[])
    for j in 1:length(DESTC) 
      push!(B5[i],[])
      for j1 in 1:length(DESTC)
        push!(B5[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B5[i][j][j1],0 )
        end
      end
    end
  end
  global B6 = [] #zeros(length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(route)
    push!(B6,[])
    for j in 1:length(DESTC1) 
      push!(B6[i],[])
      for j1 in 1:length(DESTC1)
        push!(B6[i][j],0)
      end
    end
  end
  global B7 = [] #zeros(length(route),length(DESTC))
  for i in 1:length(route)
    push!(B7,zeros(1:length(DESTC)))
  end 
  global B8 = [] #zeros(length(scenario),length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(scenario)
    push!(B8,[])
    for j in 1:length(route)
      push!(B8[i],[]) 
      for j1 in 1:length(DESTC1) 
        push!(B8[i][j],[])
        for j2 in 1:length(DESTC1)
          push!(B8[i][j][j1],0 )
        end
      end
    end
  end
  
  for r in 1:length(route)
    for i in route[r]
      B7[r][i] = 1.0
      B1[r][i] += 1.0
    end
    B6[r][route[r][end]][length(DESTC1)] = 1.0
    B6[r][length(DESTC1)][route[r][1]] = 1.0
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        B2[r][length(DESTC1)][i] = 1.0
        B2[r][i][length(DESTC1)] = 1.0
        if i != route[r][end]
          B6[r][i][route[r][findfirst(x -> x==i, route[r])+1]] = 1.0
          #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
        end
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r][i][j] = 1.0
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r][i][j][k] = 1.0
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r][i][j][k] = 1.0
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r][i][j][k] = 1.0
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end
  for r in 1:length(route), w in 1:length(scenario),i in 1:length(DESTC),j in 1:length(DESTC)
    B8[w][r][i][j]  = (1- scenario[w][i])*(1 - scenario[w][j])*B2[r][i][j]*prod(scenario[w][k]*B3[r][i][j][k]+B4[r][i][j][k]+B5[r][i][j][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != j && k != i)
  end
  for r in 1:length(route), w in 1:length(scenario),i in 1:length(DESTC)
    B8[w][r][length(DESTC1)][i ]  = (1 - scenario[w][i])*prod(scenario[w][k]*B2[r][k][i]+B2[r][i][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != i)  
    B8[w][r][i][length(DESTC1)] = ( 1- scenario[w][i] )*prod(scenario[w][k]*B2[r][i][k]+B2[r][k][i] + 1-B7[r][k] for k in 1:length(DESTC) if k != i) 
  end 
  #####################################################################
  #read(STDIN,Char)
  global firstpassnode1nonewscenario = false
  println("Start Processing nodes")

  routesused = collect(1:length(route))
  #Initialize branch-and-bound tree and queue
  #Queue = Vector{Int}()
  Queue = PriorityQueue()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,collect(1:length(scenario)),[],[],[],[],[],[],[],[])) 
  #push!(Queue,1)
  Queue[1]= 0
  #Start processing queue in a Deep First mode
  global heursol
  global globalincumbentvalue = +Inf
  global globalincumbentsolutiont = []
  global globalincumbentsolutionz = []

    
  masterx = 0
  master = []
  #always process last element of queue
  tic()
  while length(Queue)>0
    println("11") 
    #current=Queue[end]
    current= dequeue!(Queue)
    (master, masterx) =process(current)
    println("12") 
    if current == 1
      global ROOTSOL=masterx
      global ROOTTIME= TimeMain +TimePric+TimeSce
    end
    println("13") 
    global NODECOUNT += 1
    if length(Queue) > 0 && globalincumbentvalue != +Inf
      a,b=peek(Queue)
      if (globalincumbentvalue-b)/globalincumbentvalue < 0.02
        break 
      end
    end
    println("14") 
    if TimeMain +TimePric+TimeSce >= 36000
      break
    end
  end
  println("15")  
  for i in 1:length(globalincumbentsolutionz)
    if globalincumbentsolutionz[i] > 0.0001 && globalincumbentsolutionz[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  #println("result= ",globalincumbentvalue,"dsol= ", globalincumbentsolutiond, "time = ", TimeMain +TimePric+TimeSce)
  global NROUTES = sum(globalincumbentsolutionz)
  #global PERHIGH = sum(globalincumbentsolutiont)/length(globalincumbentsolutiont)
  #print(route) 
  return globalincumbentvalue,globalincumbentsolutiont
end #end solveSTODROPCEXACTALLCOLUMNS function




#
function solveSTODROPCHEUR(V,A,DESTC,PROBC,REWV,CAPV,PROBOD,REWOD)
  function process(current)
    #Initialize cenarios and solutions to be added

    global firstpassnode1nonewscenario 
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vartsolutionfixedzero = []
    vartsolutionfixedone = []
    varzsolutionfixedzero = []
    restractive = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vartsolutionfixedzero = [vartsolutionfixedzero;tree[node].vartfixedsolutionzero]
      vartsolutionfixedone = [vartsolutionfixedone;tree[node].vartfixedsolutionone]
      varzsolutionfixedzero = [varzsolutionfixedzero;tree[node].varzfixedsolutionzero]
      restractive = [restractive;tree[node].restr]
      node = tree[node].parent
      node == 0 && break
    end
    #println("parent= ", tree[current].parent)
    #println("cenariosused =",cenariosused)
    #read(STDIN,Char)
    #println("routesused =",routesused)
    #read(STDIN,Char)
    #println("varxsolutionfixedzero =",varxsolutionfixedzero)
    #println("varxsolutionfixedone =",varxsolutionfixedone) 
    #println("vartsolutionfixedzero =",vartsolutionfixedzero)
    #println("vartsolutionfixedone =",vartsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    Gomorycutsrounds = 0
    #Now formulate problem (relaxed problem now) 
    z = Vector{Variable}()
    t = Vector{Vector{Variable}}()
    y = Vector{Vector{Vector{Vector{Variable}}}}()
    const1 = Vector{ConstraintRef}(length(cenariosused))
    #const3[w in 1:length(cenariosused),i in 1:length(DESTC)]
    const3= Vector{Vector{ConstraintRef}}()
    #const4[w in 1:length(cenariosused),i in 1:length(DESTC)]
    const4= Vector{Vector{ConstraintRef}}()
    #const5[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j]
    const5 = Vector{Vector{Vector{Vector{ConstraintRef}}}}()
    #const6[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC)]
    const6 = Vector{Vector{Vector{ConstraintRef}}}()
    #const7[w in 1:length(cenariosused), r in 1:length(routesused),i in 1:length(DESTC)]
    const7 = Vector{Vector{Vector{ConstraintRef}}}()
    #const8[w in 1:length(cenariosused),i in 1:length(DESTC); scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0]
    const8 = Vector{Vector{ConstraintRef}}()
    #const9[w in 1:length(cenariosused),i in 1:length(DESTC); scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0]
    const9 = Vector{Vector{ConstraintRef}}()



    const11 = Vector{ConstraintRef}()
    #const12 = Vector{ConstraintRef}(length(scenario))

    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=3600))
    #Reserve first fixed indexes for continuous variables
    @variable(model,  s)
    @variable(model,  u[1:length(DESTC)])
    @variable(model,  v[1:length(DESTC)])
    for i in 1:length(DESTC)
      setupperbound(u[i],100000000)
      setupperbound(v[i],100000000)
    end
      
    
    for r in 1:length(routesused)
      if in(routesused[r],varzsolutionfixedzero)
        push!(z,@variable(model, lowerbound = 0, upperbound =0, basename="z[$r]"))
        JuMP.fix(z[end], 0) 
      else
        push!(z,@variable(model, lowerbound = 0, basename="z[$r]"))
      end
    end

    #for i in 1:length(tree[current].zsolution)
    #   setvalue(z[i],tree[current].zsolution[i])  
    #end
  
    @variable(model,  x[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j]>=0)
    for i in 1:length(DESTC1)
      for j in 1:length(DESTC1)
        if i != j
          setupperbound(x[i,j],1)
        end
      end 
    end
   
    for w in 1:length(cenariosused)
      push!(t,[])
      for i in 1:length(DESTC)
        push!(t[end], @variable(model, lowerbound = 0,  basename="t[$w][$i]"))
      end 
    end

    
    for w in 1:length(cenariosused)
     push!(y,[])
     for r in 1:length(routesused)
       push!(y[end], [] )
       for i in 1:length(DESTC1)
         push!(y[end][end], [] )
         for j in 1:length(DESTC1)
           if i == j
             push!(y[end][end][end], @variable(model, lowerbound = 0, upperbound =0, basename="y[$w][$r][$i][$j]") )
             JuMP.fix(y[end][end][end][j], 0) 
           else
             push!(y[end][end][end], @variable(model, lowerbound = 0, basename="y[$w][$r][$i][$j]") )
           end
         end
       end
     end 
    end

  

    for i in varxsolutionfixedzero
      setlowerbound(x[i[1],i[2]],0)
      setupperbound(x[i[1],i[2]],0) 
    end
    for i in varxsolutionfixedone
      setlowerbound(x[i[1],i[2]],1)
      setupperbound(x[i[1],i[2]],1)
    end

    for i in vartsolutionfixedzero
      JuMP.fix(t[i[1]][i[2]], 0) 
      #setlowerbound(t[i[1]][i[2]],0)
      #setupperbound(t[i[1]][i[2]],0) 
    end
    for i in vartsolutionfixedone
      JuMP.fix(t[i[1]][i[2]], 1.0) 
      #setlowerbound(t[i[1]][i[2]],1)
      #setupperbound(t[i[1]][i[2]],1)
    end

    
    
    @objective(model, Min, s + sum(PROBC[i]*u[i] + PROBOD[i]*v[i] for i in 1:length(DESTC))  )
    
    for w in 1:length(cenariosused)
      const1[w] = @constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] 
+ scenario[cenariosused[w]][length(DESTC)+i]*v[i]  for i in 1:length(DESTC))  -
sum( PRICEOD[i]*t[w][i] for i in 1:length(DESTC)) - sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y[w][r][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= 0)
    end

    @constraint(model, const2[i in 1:length(DESTC)], sum(B1[routesused[r]][i]*z[r] for r in 1:length(routesused)) == 1)
      
    @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r]][i][j]*z[r] for r in 1:length(routesused))==0)

    for w in 1:length(cenariosused)
      push!(const3,[])
      for i in 1:length(DESTC)
        push!(const3[end],@constraint(model,  t[w][i] <= 1-scenario[cenariosused[w]][i]))
      end
    end

    for w in 1:length(cenariosused)
      push!(const4,[])
      for i in 1:length(DESTC)
        push!(const4[end],@constraint(model,  t[w][i] <= scenario[cenariosused[w]][length(DESTC)+i] ))
      end
    end
   
    for w in 1:length(cenariosused)
      push!(const5,[])
      for  r in 1:length(routesused)
        push!(const5[end],[])
        for i in 1:length(DESTC)
          push!(const5[end][end],[])
          for j in 1:length(DESTC)
            if i != j
              push!(const5[end][end][end],@constraint(model, y[w][r][i][j] - z[r] +t[w][i] +t[w][j] - sum( t[w][k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scenario[cenariosused[w]][i] +1 - scenario[cenariosused[w]][j] + B2[routesused[r]][i][j] + sum( scenario[cenariosused[w]][k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 ))
            else
              push!(const5[end][end][end],@constraint(model,0==0))
            end
          end
        end
      end
    end


    for w in 1:length(cenariosused)
      push!(const6,[])
      for  r in 1:length(routesused)
        push!(const6[end],[])
        for i in 1:length(DESTC)
          push!(const6[end][end],@constraint(model, y[w][r][length(DESTC1)][i]- z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const7,[])
      for  r in 1:length(routesused)
        push!(const7[end],[])
        for i in 1:length(DESTC)
          push!(const7[end][end],@constraint(model, y[w][r][i][length(DESTC1)] - z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const8,[])
      for i in 1:length(DESTC)
        if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
          push!(const8[end],@constraint(model, t[w][i] + sum(sum( y[w][r][i][j] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
        else
          push!(const8[end],@constraint(model,0==0))
        end
      end
    end

    for w in 1:length(cenariosused)
      push!(const9,[])
      for i in 1:length(DESTC)
        if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
          push!(const9[end],@constraint(model, t[w][i] + sum(sum( y[w][r][j][i] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
        else
          push!(const9[end],@constraint(model,0==0))
        end
      end
    end
 
    const11 = [const11;restractive] #active cutting planes from parent nodes

    
   

    while true
       
      status = solve(model)
      result = getobjectivevalue(model)
      resultz = getvalue(z)
        
      global TimeMain += toq()
      tic()
      println("Solved Node ",current," problem, result =", result, "; Time is ",TimeMain, " and Queue length is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #println(model)
      #read(STDIN,Char)

      stop = false
      mincost = false
      if status != :Optimal
        #prune   
        #pos = findfirst(x -> x == current, Queue)
        #deleteat!(Queue,pos)
        #delete!(Queue, current) 
        println("Infeasible...prune ", status)
        #read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolutiont
      global globalincumbentsolutionz
      global firstpassnode1nonewscenario
     
      #ysol = getvalue(y)
      xsol1 = getvalue(x)
      xsol = zeros(length(DESTC1),length(DESTC1))

      for i=1:length(DESTC1)
        for j=1:length(DESTC1)
          if i != j    
            if xsol1[i,j]<=0.00001 
              xsol[i,j] = 0
            elseif  xsol1[i,j] >=0.99999
              xsol[i,j]= 1
            else
              xsol[i,j]=xsol1[i,j]
            end
          end
        end
      end
      tsol1 = getvalue(t)
      tsol = zeros(length(cenariosused),length(DESTC))
      for w=1:length(cenariosused)
        for i=1:length(DESTC)
            
            if tsol1[w][i]<=0.00001 
              tsol[w,i] = 0
            elseif  tsol1[w][i] >=0.99999
              tsol[w,i]= 1
            else
              tsol[w,i]=tsol1[w][i]
            end
          
        end
      end
      
      #First execute pricing

      #Collect duals necessary for route reduced cost calculation
      println("current = ",current, "firstpassnode1nonewscenario= ",firstpassnode1nonewscenario)   
      if current == 1  && firstpassnode1nonewscenario #heuristic, only include columns in node 1
      p = getdual(const1) 
      alpha = getdual(const2) 
      neta= getdual(const10)
      gamma2 = getdual(const3)
      gamma1 = getdual(const4)
      beta = getdual(const5)
      beta0i = getdual(const6)
      betai0 = getdual(const7)
      gamma3 = getdual(const8)
      gamma4 = getdual(const9)
      #println(p)
      #println(alpha)
      #println(neta)
      #Calculate parameters K4 and K3 for dynamic prog algorithm
      K3 = 1.0*zeros(length(cenariosused),length(DESTC1),length(DESTC1))
      for w in 1:length(cenariosused)
        for i in 1:length(DESTC1)
          for j in 1:length(DESTC1)
            if i != j
              K3[w,i,j]=1.0*REWV*A[DESTC1[i],DESTC1[j]]*p[w]
              if i != length(DESTC1) && scenario[cenariosused[w]][i]==0 && scenario[cenariosused[w]][length(DESTC)+i]==1
                K3[w,i,j] -= gamma3[w][i]
              end
              if j != length(DESTC1) && scenario[cenariosused[w]][j]==0 && scenario[cenariosused[w]][length(DESTC)+j]==1
               K3[w,i,j] -= gamma4[w][j]
              end
            end
          end
        end
      end
   

      #Calculate K2[w,i]
      K2 = 1.0*zeros(length(cenariosused), length(DESTC)) 
      for w in 1:length(cenariosused)
        for k in 1:length(DESTC)
          K2[w,k] = 1.0*PRICEOD[k]*p[w] - gamma1[w][k] - gamma2[w][k] - sum(sum( beta[w][r][k][j] for j in 1:length(DESTC) if j !=k) for r in 1:length(routesused)) - sum(sum( beta[w][r][j][k] for j in 1:length(DESTC) if j !=k) for r in 1:length(routesused)) - sum( beta0i[w][r][k] for r in 1:length(routesused)) - sum( betai0[w][r][k] for r in 1:length(routesused))
        
          if scenario[cenariosused[w]][k]==0 && scenario[cenariosused[w]][length(DESTC)+k]==1
            K2[w,k] -= gamma3[w][k]
            K2[w,k] -= gamma4[w][k]  
          end
        end
      end
      
      #for w in 1:length(cenariosused)
      #  for k in 1:length(DESTC)
      #    k2iw= 1.0*PRICEOD[k]*p[w] - gamma1[w][k] - gamma2[w][k]
      #    if scenario[cenariosused[w]][k]==0 && scenario[cenariosused[w]][length(DESTC)+k]==1
      #      k2iw -= gamma3[w][k]
      #      k2iw -= gamma4[w][k]  
      #    end
      #    temp = 0
      #    for r in 1:length(routesused)
      #      for j in 1:length(DESTC)
      #        if k != j
      #          temp += beta[w][r][k][j]
      #          temp += beta[w][r][j][k]
      #        end
      #      end
      #      temp += betai0[w][r][k] + betai0[w][r][k]
      #      temp -= sum( 1.0*B2[routesused[r]][k][j]*betai0[w][r][j] for j in 1:length(DESTC) if j != k)
      #      temp -= sum( 1.0*B2[routesused[r]][j][k]*beta0i[w][r][j] for j in 1:length(DESTC) if j != k)
      #      temp -= sum( sum(1.0*B3[routesused[r]][i][j][k]*beta[w][r][i][j] for i in 1:length(DESTC) if i!=j && i !=k) for j in 1:length(DESTC) if j !=k)  
      #    end
      #    println("Dif reducedcost 106 original ",w," ,",k, ", ",temp-k2iw, ", ",getdual(t[w][k]))
          #if temp - k2iw >=0.1
          #   println("Erro dif reducedcost 106 original ",w," ,",k, ", ",temp," ,", k2iw, ", ",tsol[w,k]," ,",getdual(t[w][k]))
          #   #readline()
          #end
      #  end
      #end
     
      #readline()
      #Start dynamic programming algorithm

      #Scale and round up parameters for capacity constraint
      RDemand = round.(100*(1-PROBC))
      RCAPV=100*CAPV
      #println(RDemand, RCAPV)
      minRDemand =- minimum(RDemand)+1
      M = Array{Any}(floor(Int,RCAPV + minRDemand),length(DESTC))
      #Initialize cells
      for i in 1: floor(Int,(RCAPV + minRDemand))
        for j in 1:length(DESTC)
          if i == RDemand[j]+ minRDemand
            reducedcost = alpha[j]
            temp = 1.0*zeros(1:length(cenariosused))
            for w in 1:length(cenariosused)
              #add arcs to new j customer
              #First term of min 
              if scenario[cenariosused[w]][j] == 0 
                 temp[w] += K3[w,j,length(DESTC1)]
                 temp[w] += K3[w,length(DESTC1),j]
              end
              #Second term of min
              mink = Inf
              for k in 1:length(DESTC)
                temp2 = K2[w,k]
                temp2 += temp[w]
                if k == j  && scenario[cenariosused[w]][j] == 0
                  temp2 -= K3[w,k,length(DESTC1)]
                  temp2 -= K3[w,length(DESTC1),k]
                end
                if temp2 < mink
                  mink = temp2
                end
              end
              reducedcost -= min(temp[w], mink)
            end
            M[i,j]= [[j],reducedcost,alpha[j],0,temp]
          else
            M[i,j]= [[],-Inf,0,0,zeros(1:length(cenariosused))]
          end
        end
      end


      #fill up Matrix  label (route pathsofar, costsofar, costalphasofar, costnetasofar; costK3(w)*B8sofar)
      for d in 1: floor(Int,(RCAPV + minRDemand))
        for i in 1:length(DESTC)
          for j in 1:length(DESTC)
            if i != j  
              if length(M[d,i][1]) != 0 && findfirst(x -> x == j, M[d,i][1]) == 0 && M[d,i][1][end] == i && d - minRDemand + RDemand[j] <= RCAPV   
                #calculate reduced cost for extended visit from i to j with demand d
                reducedcost = M[d,i][3] - M[d,i][4]
                reducedcost += alpha[j]
                reducedcost -= neta[i,j]
                temp = 1.0*zeros(1:length(cenariosused))
                for w in 1:length(cenariosused)
                  temp[w] =  M[d,i][5][w]
                  #Eliminate previous links between last available customer towards depot and
                  #add arcs to new j customer 
                  if scenario[cenariosused[w]][j] == 0
                    for k in M[d,i][1]
                      pos = findfirst(x -> x == k, M[d,i][1])
                      if scenario[cenariosused[w]][k] == 0 && (k== M[d,i][1][end] || prod(scenario[cenariosused[w]][M[d,i][1][pos+1:end]] == 1) )
                        temp[w] -= K3[w,k,length(DESTC1)]
                        temp[w] += K3[w,k,j]
                        temp[w] += K3[w,j,length(DESTC1)]
                        break
                      end
                    end
                  end
                  #Second term of min
                  mink = Inf
                  for k in 1:length(DESTC)
                    temp2 = K2[w,k]+temp[w]
                    for k2 in vcat(M[d,i][1],[j])
                      if in(k,vcat(M[d,i][1],[j])) && k2 != k && scenario[cenariosused[w]][k2] == 0 && scenario[cenariosused[w]][k] == 0
                        pos = findfirst(x -> x == k, vcat(M[d,i][1],[j]))
                        pos2 = findfirst(x -> x == k2, vcat(M[d,i][1],[j]))
                        if pos < pos2 
                          if prod(scenario[cenariosused[w]][M[d,i][1][pos+1:pos2-1]]) == 1
                            temp2 -= K3[w,k,k2]
                          end 
                        elseif  pos > pos2
                          if prod(scenario[cenariosused[w]][M[d,i][1][pos2+1:pos-1]]) == 1
                            temp2 -= K3[w,k2,k]
                          end  
                        end
                      end
                    end
                    if temp2 < mink
                      mink = temp2
                    end
                  end 
                  reducedcost -= min(temp[w], mink)
                end
                if reducedcost > M[floor(Int,d+ RDemand[j]),j][2]
                 M[floor(Int,d+ RDemand[j]),j][1]= vcat(M[d,i][1],[j])
                 M[floor(Int,d+ RDemand[j]),j][2]= reducedcost
                 M[floor(Int,d+ RDemand[j]),j][3]= M[d,i][3]+ alpha[j]
                 M[floor(Int,d+ RDemand[j]),j][4]= M[d,i][4]+ neta[i,j]
                 M[floor(Int,d+ RDemand[j]),j][5]= temp
                end 
              end
            end
          end
        end
      end
      
      temproutes=[]
      temp = []
      maxreducedcost = -Inf
      for i in 1:length(DESTC)
        for d =floor(Int,(RCAPV + minRDemand)):-1:1
           #println("d= ",d,"; i= ",i," : ", M[d,i])           
           if M[d,i][2] != -Inf && M[d,i][2]>0.01
             ins=true
             for ind in 1:length(routesused)
               if length(M[d,i][1])==length(route[routesused[ind]]) && setdiff(M[d,i][1],route[routesused[ind]]) == []
                 ins=false
                 break
               end
             end
             for ind in 1:length(temproutes)
               if length(M[d,i][1])==length(temproutes[ind]) && setdiff(M[d,i][1],temproutes[ind]) == []
                 ins=false
                 break
               end
             end
 
             if  ins &&   M[d,i][2] > maxreducedcost
               maxreducedcost = M[d,i][2]
               temp =  M[d,i][1]
               #push!(temproutes, M[d,i][1])
               break
             end
           end
        end
        #readline()
      end
      if temp != []
        push!(temproutes,temp)
      end
      #println(temproutes," ,",maxreducedcost)
      #routes with positive reduced  cost should be inserted
      if length(temproutes) != 0
        mincost = true
        for i = 1:length(temproutes)
          newroute = temproutes[i]
          #println("new routes= ", temproutes[i])
          push!(route,newroute)
          push!(routesused,size(route,1))  
          #Update Current node information on routes used
          push!(tree[current].addedsolutions,size(route,1))  
          global ROUTCOUNT += 1
      
          #Update B parameters
          push!(B1,1.0*zeros(1:length(DESTC)))
          push!(B2,[])
          for j in 1:length(DESTC1) 
            push!(B2[end],[])
            for j1 in 1:length(DESTC1)
              push!(B2[end][j],0 )
            end
          end
          push!(B3,[])
          for j in 1:length(DESTC) 
            push!(B3[end],[])
            for j1 in 1:length(DESTC)
              push!(B3[end][j],[])
              for j2 in 1:length(DESTC)
                push!(B3[end][j][j1],0 )
              end
            end
          end
          push!(B4,[])
          for j in 1:length(DESTC) 
            push!(B4[end],[])
            for j1 in 1:length(DESTC)
              push!(B4[end][j],[])
              for j2 in 1:length(DESTC)
                push!(B4[end][j][j1],0 )
              end
            end
          end
          push!(B5,[])
          for j in 1:length(DESTC) 
            push!(B5[end],[])
            for j1 in 1:length(DESTC)
              push!(B5[end][j],[])
              for j2 in 1:length(DESTC)
                push!(B5[end][j][j1],0 )
              end
            end
          end
          push!(B6,[])
          for j in 1:length(DESTC1) 
            push!(B6[end],[])
            for j1 in 1:length(DESTC1)
              push!(B6[end][j],0)
            end
          end
          push!(B7,zeros(1:length(DESTC)))
          for i in 1:length(scenario)
            push!(B8[i],[]) 
            for j1 in 1:length(DESTC1) 
              push!(B8[i][end],[])
              for j2 in 1:length(DESTC1)
                push!(B8[i][end][j1],0 )
              end
            end
          end
          r = size(route,1)
          for i in route[r]
            B7[r][i] = 1.0
            B1[r][i] += 1.0
          end
          B6[r][route[r][end]][length(DESTC1)] = 1.0
          B6[r][length(DESTC1)][route[r][1]] = 1.0
          for i in 1:length(DESTC)
            if findfirst(x -> x==i, route[r]) != 0
              B2[r][length(DESTC1)][i] = 1.0
              B2[r][i][length(DESTC1)] = 1.0
              if i != route[r][end]
                B6[r][i][route[r][findfirst(x -> x==i, route[r])+1]] = 1.0
                #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
              end
              #println("B1[ ",r,",",i," ]= ", B1[r,i])
              #read(STDIN,Char)
            end 
            for j in 1:length(DESTC)
              if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
                B2[r][i][j] = 1.0
                #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
                #read(STDIN,Char)
              end
              #println()
              for k in 1:length(DESTC)
                if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
                  B3[r][i][j][k] = 1.0
                  #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
                  #read(STDIN,Char)
                end
                #println()
                if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
                  B4[r][i][j][k] = 1.0
                  #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
                  #read(STDIN,Char)
                end
                #println()
                if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
                  B5[r][i][j][k] = 1.0
                  #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
                  #read(STDIN,Char)
                end
                #println()
              end
            end
          end
          for  w in 1:length(scenario),i in 1:length(DESTC),j in 1:length(DESTC)
            B8[w][r][i][j]  = (1- scenario[w][i])*(1 - scenario[w][j])*B2[r][i][j]*prod(scenario[w][k]*B3[r][i][j][k]+B4[r][i][j][k]+B5[r][i][j][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != j && k != i)
          end
          for  w in 1:length(scenario),i in 1:length(DESTC)
            B8[w][r][length(DESTC1)][i ]  = (1 - scenario[w][i])*prod(scenario[w][k]*B2[r][k][i]+B2[r][i][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != i)  
            B8[w][r][i][length(DESTC1)] = ( 1- scenario[w][i] )*prod(scenario[w][k]*B2[r][i][k]+B2[r][k][i] + 1-B7[r][k] for k in 1:length(DESTC) if k != i) 
          end 
          #Define new variables and constraints
          tes = [1.0*B1[r][1]]
          tes1 = [const2[1]]
          for i in 2:length(DESTC)
            push!(tes,1.0*B1[r][i])
            push!(tes1,const2[i])
          end
          for i in 1:length(DESTC1)
            for j in 1:length(DESTC1)
              if i != j
                push!(tes, -1.0*B6[r][i][j])
                push!(tes1, const10[i,j])
              end
            end
          end
          push!(z,@variable(model,  lowerbound = 0,  basename="z[$r]",objective = 0.0,inconstraints = tes1, coefficients = tes ))
      
   
          for w in 1:length(cenariosused)
            push!(y[w], [] )
            for i in 1:length(DESTC1)
              push!(y[w][end], [] )
              for j in 1:length(DESTC1)
                if i == j
                  push!(y[w][end][i], @variable(model, lowerbound = 0, upperbound =0, basename="y[$w][$r][$i][$j]") )
                  JuMP.fix(y[w][end][i][j], 0)
                else
                  tes1 = [const1[w]]
                  tes = [-1.0*REWV*A[DESTC1[i],DESTC1[j]]]
                  if i != length(DESTC1) && scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
                    push!(tes1, const8[w][i])
                    push!(tes, 1.0)
                  end
                  if j != length(DESTC1) && scenario[cenariosused[w]][length(DESTC)+j] == 1 && scenario[cenariosused[w]][j]==0
                    push!(tes1, const9[w][j])
                    push!(tes, 1.0)             
                  end
                  #push!(y[w][end][i], @variable(model, lowerbound = 0, upperbound =1, basename="y[$w][$r][$i][$j]",objective = 0.0, inconstraints = tes1, coefficients = tes) )
                  push!(y[w][end][i], @variable(model, lowerbound = 0, basename="y[$w][$r][$i][$j]",objective = 0.0, inconstraints = tes1, coefficients = tes) )
                end
              end
            end
          end
  
          r = length(routesused)
          for w in 1:length(cenariosused)
            push!(const5[w],[])
            for i in 1:length(DESTC)
              push!(const5[w][end],[])
              for j in 1:length(DESTC)
                if i != j
                  push!(const5[w][end][end],@constraint(model, y[w][r][i][j] - z[r] +t[w][i] +t[w][j] - sum( t[w][k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scenario[cenariosused[w]][i] +1 - scenario[cenariosused[w]][j] + B2[routesused[r]][i][j] + sum( scenario[cenariosused[w]][k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 ))
                else
                  push!(const5[w][end][end],@constraint(model,0==0))
                end
              end
            end
          end

          for w in 1:length(cenariosused)
            push!(const6[w],[])
            for i in 1:length(DESTC)
              push!(const6[w][end],@constraint(model, y[w][r][length(DESTC1)][i]- z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
            end
          end
 

          for w in 1:length(cenariosused)
            push!(const7[w],[])
            for i in 1:length(DESTC)
              push!(const7[w][end],@constraint(model, y[w][r][i][length(DESTC1)] - z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
            end
          end
        end
      end #end temproutes insert
      global TimePric += toq()
      tic()
      end  #end current == 1 heuristic

      #println("mincost =", mincost)
      if !mincost 
        #println("after mincost not true no column inserted")
        if  result > globalincumbentvalue + 0.0001
         #prune  since a lower bound is already worse then incumbent 
         #pos = findfirst(x -> x == current, Queue)
         #deleteat!(Queue,pos)
         #delete!(Queue, current)
         #println("Prune by Bound =", result, ", ", globalincumbentvalue)
         #read(STDIN,Char) 
         stop = true
         break
        end

        
        ssol = getvalue(s)
        usol= getvalue(u)
        vsol= getvalue(v)
        if current == 1  || (maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(tsol-floor.(tsol)) <= 0.00001)  #heuristic inserts all scenarios in node 1
          #new integer solution
          #println("Integer Solution Found")
          ##############################################
          #Formulate Scenario Separation Problem and Solve 
          #Get First Stage results
          if current != 1
            for i in 1:length(resultz) #check for error and display
              if  resultz[i] > 0.0001 && resultz[i] < 0.99
                println("SOLUCAO NAO INTEIRA!! ",i, ",",resultz[i])
                read(STDIN,Char)
              end
            end
          end


          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))
                          
          #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:2*length(DESTC)], Bin)
          @variable(model2, t2[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1)
          @variable(model2, obj1)
          @variable(model2, obj2) 
 
          #@objective(model2, Min, - ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC))  -  sum(scen[length(DESTC)+i]*vsol[i] for i in 1:length(DESTC)) + sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) + sum(PRICEOD[i]*t2[i] for i in 1:length(DESTC)))
          @objective(model2, Min, obj1 + obj2)
          @constraint(model2, obj1a, obj1 >=   ssol + sum(scen[i]*usol[i] for i in 1:length(DESTC))  +  sum(scen[length(DESTC)+i]*vsol[i] for i in 1:length(DESTC)) )
          @constraint(model2, obj2a, obj2 >=   sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) + sum(PRICEOD[i]*t2[i] for i in 1:length(DESTC)))
          @constraint(model2, obj3a, obj1 <= obj2)

          @constraint(model2, const2b1[i in 1:length(DESTC)], t2[i] <= scen[length(DESTC)+i])
          @constraint(model2, const2b2[i in 1:length(DESTC)], t2[i] <= 1-scen[i])
          @constraint(model2, const2b3[i in 1:length(DESTC)], scen[length(DESTC)+i]+scen[i] <=1)

          @constraint(model2, const2a1[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= resultz[r] )
          @constraint(model2, const2a2[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a3[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j] - t2[j])  )
          @constraint(model2, const2a4[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= B2[routesused[r]][i][j] )
          @constraint(model2,  const2a5[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC), k in 1:length(DESTC); i != j && k !=i && k !=j  && resultz[r] > 0.001], y2[r,i,j] <=  (scen[k]+ t2[k])*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1-B7[routesused[r]][k] )
          @constraint(model2,  const2a6[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,i,j] - resultz[r] +t2[i] +t2[j] -
sum( t2[k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scen[i] +1 - scen[j] + B2[routesused[r]][i][j] + sum( scen[k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 )

          @constraint(model2, const2a7[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= resultz[r] )
          @constraint(model2, const2a8[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a9[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= B7[routesused[r]][i]  )
          @constraint(model2,  const2a10[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,length(DESTC1),i] <=  (scen[j]+ t2[j])*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1-B7[routesused[r]][j]) 
          @constraint(model2,  const2a11[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] - resultz[r] +t2[i] - sum( t2[j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scen[i]  + B7[routesused[r]][i] + sum( scen[j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 )

          @constraint(model2, const2a12[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= resultz[r] )
          @constraint(model2, const2a13[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i] - t2[i]) )
          @constraint(model2, const2a14[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= B7[routesused[r]][i]  )
          @constraint(model2,  const2a15[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i != j  && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <=  (scen[j]+ t2[j])*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1-B7[routesused[r]][j] )
          @constraint(model2,  const2a16[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] - resultz[r] +t2[i] - sum( t2[j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scen[i]  + B7[routesused[r]][i] + sum( scen[j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 )

           
          #println("Will solve Separation Problem for integer solution")
          
          solve(model2)
          #println("solve model 2 ",getobjectivevalue(model2))
          #println(getvalue(obj1))
          #println(getvalue(obj2), ", ", result)
          #println(getvalue(obj1)-getvalue(obj2))
          #println(getvalue(scen))
          #read(STDIN,Char)
          scenarionew = round.(getvalue(scen))
          flag=false 
 
          if getvalue(obj1)-getvalue(obj2) <= -0.0001
            println("vou gerar novo cenario")
          else 
            println("Nao vou gerar novo cenario")
          end
          if current == 1 && getvalue(obj1)-getvalue(obj2) <= -0.0001 #&& flag == false#-0.05  <= -0.0001
            println("node 1 e vou gerar novo cenario")
            println("ADD NEW SCENARIO ",getvalue(scen), " get= ", getobjectivevalue(model2))
            #read(STDIN,Char)
            #Add new scenario and continue
  
            push!(scenario,scenarionew)
            #Update Current node information on scenarios used
            push!(tree[current].addedcenarios,size(scenario,1))
            push!(cenariosused,size(scenario,1))
            #Create new variables
            w= length(cenariosused)
            push!(t,[])
            for i in 1:length(DESTC)
              push!(t[end], @variable(model, lowerbound = 0, basename="t[$w][$i]"))
            end 
            push!(y,[])
            for r in 1:length(routesused)
              push!(y[end], [] )
              for i in 1:length(DESTC1)
                push!(y[end][r], [] )
                for j in 1:length(DESTC1)
                  if i == j
                    push!(y[end][r][i], @variable(model, lowerbound = 0, upperbound =0, basename="y[$w][$r][$i][$j]") )
                    JuMP.fix(y[end][r][i][j], 0)
                  else
                    push!(y[end][r][i], @variable(model, lowerbound = 0, basename="y[$w][$r][$i][$j]") )
                  end
                end
              end
            end 
  
            #Create new constraints
            push!(const1,@constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] 
+ scenario[cenariosused[w]][length(DESTC)+i]*v[i]  for i in 1:length(DESTC))  -
sum( PRICEOD[i]*t[w][i] for i in 1:length(DESTC)) - sum(sum(sum(1.0*REWV*A[DESTC1[i],DESTC1[j]]*y[w][r][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= 0) )

    
            push!(const3,[])
            for i in 1:length(DESTC)
              push!(const3[end],@constraint(model,  t[w][i] <= 1-scenario[cenariosused[w]][i]))
            end
        
            push!(const4,[])
            for i in 1:length(DESTC)
              push!(const4[end],@constraint(model,  t[w][i] <= scenario[cenariosused[w]][length(DESTC)+i] ))
            end
    
   
   
            push!(const5,[])
            for  r in 1:length(routesused)
              push!(const5[end],[])
              for i in 1:length(DESTC)
                push!(const5[end][r],[])
                for j in 1:length(DESTC)
                  if i != j
                    push!(const5[end][r][i],@constraint(model, y[w][r][i][j] - z[r] +t[w][i] +t[w][j] - sum( t[w][k]*B3[routesused[r]][i][j][k] for k in 1:length(DESTC) if k !=i && k !=j) >= 1 - scenario[cenariosused[w]][i] +1 - scenario[cenariosused[w]][j] + B2[routesused[r]][i][j] + sum( scenario[cenariosused[w]][k]*B3[routesused[r]][i][j][k] + B4[routesused[r]][i][j][k]+ B5[routesused[r]][i][j][k] +1 - B7[routesused[r]][k] for k in 1:length(DESTC) if k !=i && k !=j) - length(DESTC) -1 ))
                  else
                    push!(const5[end][r][i],@constraint(model,0==0))
                  end
                end
              end
            end

            push!(const6,[])
            for  r in 1:length(routesused)
              push!(const6[end],[])
              for i in 1:length(DESTC)
                push!(const6[end][r],@constraint(model, y[w][r][length(DESTC1)][i]- z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][j][i] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][j][i] + B2[routesused[r]][i][j] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
              end
            end
              
            push!(const7,[])
            for  r in 1:length(routesused)
              push!(const7[end],[])
              for i in 1:length(DESTC)
                push!(const7[end][r],@constraint(model, y[w][r][i][length(DESTC1)] - z[r] +t[w][i] - sum( t[w][j]*B2[routesused[r]][i][j] for j in 1:length(DESTC) if j !=i) >= 1 - scenario[cenariosused[w]][i]  + B7[routesused[r]][i] + sum( scenario[cenariosused[w]][j]*B2[routesused[r]][i][j] + B2[routesused[r]][j][i] +1 - B7[routesused[r]][j] for j in 1:length(DESTC) if j !=i) - length(DESTC) -1 ))
              end
            end
  
            push!(const8,[])
            for i in 1:length(DESTC)
              if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
                push!(const8[end],@constraint(model, t[w][i] + sum(sum( y[w][r][i][j] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
              else
                  push!(const8[end],@constraint(model,0==0))
              end
            end


            push!(const9,[])
            for i in 1:length(DESTC)
              if scenario[cenariosused[w]][length(DESTC)+i] == 1 && scenario[cenariosused[w]][i]==0
                push!(const9[end],@constraint(model, t[w][i] + sum(sum( y[w][r][j][i] for j in 1:length(DESTC1) if j != i) for r in 1:length(routesused)) == 1))
              else
                push!(const9[end],@constraint(model,0==0))
              end
            end
            
            #Calculate additional parameters B8
            global B8    #zeros(length(scenario),length(route),length(DESTC1),length(DESTC1))
            push!(B8,[])
            for j in 1:length(route)
              push!(B8[end],[]) 
              for j1 in 1:length(DESTC1) 
                push!(B8[end][j],[])
                for j2 in 1:length(DESTC1)
                  push!(B8[end][j][j1],0 )
                end
              end
            end
            
            
            for r in 1:length(route),i in 1:length(DESTC),j in 1:length(DESTC)
              B8[cenariosused[w]][r][i][j]  = (1- scenario[cenariosused[w]][i])*(1 - scenario[cenariosused[w]][j])*B2[r][i][j]*prod(scenario[cenariosused[w]][k]*B3[r][i][j][k]+B4[r][i][j][k]+B5[r][i][j][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != j && k != i)
            end
            for r in 1:length(route),i in 1:length(DESTC)
              B8[cenariosused[w]][r][length(DESTC1)][i ]  = (1 - scenario[cenariosused[w]][i])*prod(scenario[cenariosused[w]][k]*B2[r][k][i]+B2[r][i][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != i)  
              B8[cenariosused[w]][r][i][length(DESTC1)] = ( 1- scenario[cenariosused[w]][i] )*prod(scenario[cenariosused[w]][k]*B2[r][i][k]+B2[r][k][i] + 1-B7[r][k] for k in 1:length(DESTC) if k != i) 
            end 
              
            stop=false
          elseif  firstpassnode1nonewscenario && (maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(tsol-floor.(tsol)) <= 0.00001) # end getobj < +.05
            println("ainda nao terminei node 1 inserir cenario e solucao é inteira. Prune a vai pra outro node ")
            #println("No new scenario on integer solution")
            ##############################################
            #Verify if it can be new incubent
            if result < globalincumbentvalue
              #println("New Incumbent Found")
              globalincumbentvalue = result
              globalincumbentsolutiont = tsol 
              globalincumbentsolutionz = resultz
            end 
            #pos = findfirst(x -> x == current, Queue)
            #deleteat!(Queue,pos)
            #delete!(Queue, current)
            #println("Will break after integer solution")
            stop=true
            break
          elseif  current == 1
            println("estou node 1, nao tem cebario, nao é inteiro. Se nao terminou cenario, agora termina")
            if !firstpassnode1nonewscenario
              println("termina cenario")
              global firstpassnode1nonewscenario = true
              stop = false
            else
              println("ja tinha terminado cenario, vou fazer branch")
              stop = true
              #Insert Gomory custs for variables  "x"
              if  Gomorycutsrounds == 0 && current == 1
                println("...will insert Gomory cut")
                Gomorycutsrounds += 1
                cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
                rowv2 = []

                for i in 1:length(DESTC1)
                  for j in 1:length(DESTC1)
                    if i != j
                      cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
                      #println("test x[",i," ,",j,"]") 
                      if cbasis[linearindex(x[i,j])] == :Basic && xsol[i,j]-floor(xsol[i,j]) > 0.00001  #is basic variable and is fractional, so insert gomory cut
                        #println("x[",i," ,",j,"] is elegible ", cbasis[linearindex(x[i,j])], ",",xsol[i,j],",",in([i;j],varxsolutionfixedzero)," ,",in([i;j],varxsolutionfixedone)) 
                        #read(STDIN,Char)
                        row = 1.0*zeros(length(cbasis))
                        #identify row where variable is basic
                        for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
                          #print("test x[",i," ,",j,"] in constraint ", k)
                          ci = model.internalModel.inner
                          ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
                          #println("(",linearindex(x[i,j]),",",row[linearindex(x[i,j])],",")
                          #if  row[linearindex(x[i,j])] != 0  read(STDIN,Char) end 
                          if row[linearindex(x[i,j])] == 1
                            #println("found x[",i," ,",j,"] in ", k) 
                            break   
                          end
                        end
                        if row == 1.0*zeros(length(cbasis))
                          println("ERROR GOMORY CUT X SECTION!")
                          read(STDIN,Char)
                        else #insert cut
                          #push!(const11,@constraint(model, sum(min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(x[i,j]) && row[k] != 0 ) >= 1) )
                          for k in 1:length(row)
                            if k != linearindex(x[i,j]) && row[k] != 0
                              if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
                                row[k] =  row[k]/(1- xsol[i,j]+floor(xsol[i,j]) )
                              elseif k <= 2*length(DESTC)+1  && row[k] > 0 
                                row[k] =  row[k]/(xsol[i,j]-floor(xsol[i,j]) )
                              elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= xsol[i,j]-floor(xsol[i,j]) 
                                row[k] =  (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) )
                              elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > xsol[i,j]-floor(xsol[i,j]) 
                                row[k] =  (1-row[k]+floor(row[k]))/(1-xsol[i,j]+floor(xsol[i,j]) )
                              else
                                println("ERROR! NOT COVERED GOMMORY X")
                                read(STDIN,Char)
                              end
                              #row[k] = min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )
                            end
                          end 
                          row[linearindex(x[i,j])]=0
                          #push!(row,xsol[i,j]-floor(xsol[i,j]) )
                          push!(rowv2,row)
                          #println("GOMORY CUT FOR VAR X!")
                          #read(STDIN,Char)
                        end 
                      end
                    end
                  end
                end

                       
                #Only at the end insert new constraints
                for v = 1:length(rowv2)
                  push!(const11,@constraint(model, sum( rowv2[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv2[v][k] != 0 ) >= 1 )  )
                  push!(tree[current].restr,const11[end])
                end
              end           #End insert gomory custs        
        
        
              #time for branching
              println("time for branching node 1")
              if maximum(xsol-floor.(xsol)) > 0.00001
                 
                if length(vartsolutionfixedzero) != 0 || length(vartsolutionfixedone) != 0
                  #println("d ja estava fixed e veio xij frac")
                  #println(vartsolutionfixedzero)
                  #println(vartsolutionfixedone)
                  #read(STDIN,Char)
                end 
                 
                a = maximum(xsol-floor.(xsol))
                f=0
                g=0
                for i in 1:length(DESTC1),j in 1:length(DESTC1)
                  if i != j
                    if xsol[i,j]-floor.(xsol[i,j])== a
                      f = i
                      g  = j
                      break
                    end
                  end
                end
                 
                #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
                b =[f;g] 
                push!(tree,TreeNode(current,[],[],[],[],[b],[],[],[],resultz,[],const11))
                push!(tree[current].children,length(tree))
                #push!(Queue,length(tree))
                Queue[length(tree)]=result

                
                #push!(tree,TreeNode(current,[],[],[],[b],[],[],[],tsol,resultz,[],const11))
                push!(tree,TreeNode(current,[],[],[],[b],[],[],[],[],resultz,[],const11)) 
                #current = length(tree)
                push!(tree[current].children,length(tree))
                ##push!(Queue,length(tree))
                Queue[length(tree)]=result

                ##pos = findfirst(x -> x == current, Queue)
                ##deleteat!(Queue,pos)
                ##delete!(Queue,current)

                ##current=Queue[end]
          
                ##Gomorycutsrounds = 0
                ##current=length(tree)

                #global NODECOUNT += 1
                #setlowerbound(x[f,g],0)
                #setupperbound(x[f,g],0)
                #push!(varxsolutionfixedzero,b) 
              elseif maximum(tsol-floor.(tsol)) > 0.00001
                println("5a") 
                stop=false
                a = round.(tsol)
                for i in 1:length(cenariosused),j in 1:length(DESTC)
                  if  !in([i;j],vartsolutionfixedzero) && !in([i;j],vartsolutionfixedone)
                    if a[i,j] <0.5
                      JuMP.fix(t[i][j],0)
                    else
                      JuMP.fix(t[i][j],1.0)
                    end
                  end
                end
                for i = 1:length(DESTC1)
                  for j = 1:length(DESTC1)
                    if i != j
                      setlowerbound(x[i,j],xsol[i,j])
                      setupperbound(x[i,j],xsol[i,j])
                    end
                  end
                end
                println("t nao era inteiro, vou arredondar e ficar no mesmo node")
              end #
            end #end firstpassnode1nonewscenario
          end #end current = 1  
        else #Not a Integer solution
          stop = true
          println("nao era sol inteira, vou fazer branching")
          #If root note, verify scenario to be inserted and continue solving master problem
             #but this version ot inserting scenarios for root node
          # Reduced Cost variable fixing for z
          zdual = getdual(z)
          for i in 1:length(zdual)
            if result + zdual[i] > globalincumbentvalue && !in(i,varzsolutionfixedzero)
              println("node ",current," z solution for route ",routesused[i], " will not improve")
              #read(STDIN,Char)
              #if resultz[i] >= 0.1 && resultz[i] <= 0.9
              #  println("erro z[", i,"]")
              #  read(STDIN,Char)
              #end
              #setlowerbound(z[i],0)
              #setupperbound(z[i],0)
              JuMP.fix(z[i], 0) 
              push!(varzsolutionfixedzero,routesused[i])
              push!(tree[current].varzfixedsolutionzero,routesused[i])
              #Should still eliminate constraints not needed here  
            end
          end
          #End Section for reduced Cost variable fixing

          #insert rounded capacity cuts if possible (deterministic)
          #println("verify rounded capacity vi")
          for num1 in 1:20    #try at most 20 times
            #create partition
            setpart = sample(collect(1:length(DESTC)), rand(5:length(DESTC)), replace = false)
            num2 = length(setpart)
            for num3 in 1:num2-3
              setpart = setdiff(setpart,[setpart[1]])
              #verify if partition violates x solution
              capacity = 2*ceil(prod((1-PROBOD[i]) for i in setpart)/CAPV)
              sumtemp=0
              for i in setpart
                sumtemp += sum(xsol[i,j] for j in 1:length(DESTC1))
              end
              if sumtemp < capacity
                #insert rounded capacity
                println("vi inserted")
                push!(const11, @constraint(model, sum(sum(x[i,j] for i in setpart ) for j=1:length(DESTC1) if i != j  ) >= capacity ) )
                push!(tree[current].restr,const11[end])
                break
              end
            end 
          end #end 20 times most

          #Insert Gomory custs for variables  "x"
          if  Gomorycutsrounds == 0 && current == 1
            println("...will insert Gomory cut")
            Gomorycutsrounds += 1
            cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
            rowv2 = []

            for i in 1:length(DESTC1)
              for j in 1:length(DESTC1)
                if i != j
                  cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
                  #println("test x[",i," ,",j,"]") 
                  if cbasis[linearindex(x[i,j])] == :Basic && xsol[i,j]-floor(xsol[i,j]) > 0.00001  #is basic variable and is fractional, so insert gomory cut
                    #println("x[",i," ,",j,"] is elegible ", cbasis[linearindex(x[i,j])], ",",xsol[i,j],",",in([i;j],varxsolutionfixedzero)," ,",in([i;j],varxsolutionfixedone)) 
                    #read(STDIN,Char)
                    row = 1.0*zeros(length(cbasis))
                    #identify row where variable is basic
                    for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
                      #print("test x[",i," ,",j,"] in constraint ", k)
                      ci = model.internalModel.inner
                      ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
                      #println("(",linearindex(x[i,j]),",",row[linearindex(x[i,j])],",")
                      #if  row[linearindex(x[i,j])] != 0  read(STDIN,Char) end 
                      if row[linearindex(x[i,j])] == 1
                        #println("found x[",i," ,",j,"] in ", k) 
                        break   
                      end
                    end
                    if row == 1.0*zeros(length(cbasis))
                      println("ERROR GOMORY CUT X SECTION!")
                      read(STDIN,Char)
                    else #insert cut
                      #push!(const11,@constraint(model, sum(min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(x[i,j]) && row[k] != 0 ) >= 1) )
                      for k in 1:length(row)
                        if k != linearindex(x[i,j]) && row[k] != 0
                          if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
                            row[k] =  row[k]/(1- xsol[i,j]+floor(xsol[i,j]) )
                          elseif k <= 2*length(DESTC)+1  && row[k] > 0 
                            row[k] =  row[k]/(xsol[i,j]-floor(xsol[i,j]) )
                          elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= xsol[i,j]-floor(xsol[i,j]) 
                            row[k] =  (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) )
                          elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > xsol[i,j]-floor(xsol[i,j]) 
                            row[k] =  (1-row[k]+floor(row[k]))/(1-xsol[i,j]+floor(xsol[i,j]) )
                          else
                           println("ERROR! NOT COVERED GOMMORY X")
                            read(STDIN,Char)
                          end
                          #row[k] = min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )
                        end
                      end 
                      row[linearindex(x[i,j])]=0
                      #push!(row,xsol[i,j]-floor(xsol[i,j]) )
                      push!(rowv2,row)
                      #println("GOMORY CUT FOR VAR X!")
                      #read(STDIN,Char)
                    end 
                  end
                end
              end
            end
                       
            #Only at the end insert new constraints
            for v = 1:length(rowv2)
              push!(const11,@constraint(model, sum( rowv2[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv2[v][k] != 0 ) >= 1 )  )
              push!(tree[current].restr,const11[end])
            end
          end           #End insert gomory custs        
        
        
          #time for branching
          println("time for branching")
          if maximum(xsol-floor.(xsol)) > 0.00001
             
            if length(vartsolutionfixedzero) != 0 || length(vartsolutionfixedone) != 0
              #println("d ja estava fixed e veio xij frac")
              #println(vartsolutionfixedzero)
              #println(vartsolutionfixedone)
              #read(STDIN,Char)
            end 
            println("2") 
            a = maximum(xsol-floor.(xsol))
            f=0
            g=0
            for i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                if xsol[i,j]-floor.(xsol[i,j])== a
                  f = i
                  g  = j
                  break
                end
              end
            end
            println("3") 
            #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
            b =[f;g] 
            push!(tree,TreeNode(current,[],[],[],[],[b],[],[],[],resultz,[],const11))
            push!(tree[current].children,length(tree))
            #push!(Queue,length(tree))
            Queue[length(tree)]=result

             println("4") 
            #push!(tree,TreeNode(current,[],[],[],[b],[],[],[],tsol,resultz,[],const11))
            push!(tree,TreeNode(current,[],[],[],[b],[],[],[],[],resultz,[],const11)) 
            #current = length(tree)
            push!(tree[current].children,length(tree))
            ##push!(Queue,length(tree))
            Queue[length(tree)]=result

            ##pos = findfirst(x -> x == current, Queue)
            ##deleteat!(Queue,pos)
            ##delete!(Queue,current)

            ##current=Queue[end]
          
            ##Gomorycutsrounds = 0
            ##current=length(tree)

            #global NODECOUNT += 1
            #setlowerbound(x[f,g],0)
            #setupperbound(x[f,g],0)
            #push!(varxsolutionfixedzero,b) 
          elseif maximum(tsol-floor.(tsol)) > 0.00001
            println("5") 
            #heuristic solution; will round up t variables and fix values
            stop=false
            a = round.(tsol)
            for i in 1:length(cenariosused),j in 1:length(DESTC)
              if  !in([i;j],vartsolutionfixedzero) && !in([i;j],vartsolutionfixedone)
                if a[i,j] <0.5
                  JuMP.fix(t[i][j],0)
                else
                  JuMP.fix(t[i][j],1.0) 
                end
              end
            end
            for i = 1:length(DESTC1)
              for j = 1:length(DESTC1)
                if i != j
                  setlowerbound(x[i,j],xsol[i,j])
                  setupperbound(x[i,j],xsol[i,j])
                end
              end
            end
          end #
        end
          
        global TimeSce += toq()
        tic()
      end # !mincost
      
      stop && break
       
    end #end while true 'solve)
    return model, result 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate  Initial Routes, Scenarios, Prices and Incumbent solution
  ###############################################
  #Create vectors for all possible scenarios
  #global scenarioall = Vector{Vector{Int}}() 
  #fillscenario(scenarioall,zeros(length(DESTC)),length(DESTC),1)
  #println("All possible Scenario vectors created of size ", length(scenarioall) )
  
  #Create vectors for initial routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(DESTC)
    push!(route, [s])
  end
  println("Initial Route vectors created of size ", length(route))
  #################################################
 


  ###############################################
  #Create Vectors for initial scenarios 
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  push!(scenario, vcat(ones(1:length(DESTC)),zeros(1:length(DESTC)) )  )   #Scenario 1 no customer orders. 
  push!(scenario, vcat(zeros(1:length(DESTC)),ones(1:length(DESTC)) )  )   #Scenario 1 all customer orders with ODs.
  ################################################
  #println(scenario)
  #readline()
  

  ################################################
  #Bypass REWOD and define new prices to pay ODs
  global PRICEOD
  PRICEOD = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > 1.0*REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = 1.0*REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  ################################################

   ################################################
  #Initial Incumbent solution
  routetemp = [0]
  for num in DESTC
    minj = 0
    sumj = Inf
    for j in setdiff(collect(1:length(DESTC)), routetemp)
      sumtemp=0
      for i in routetemp
        if i == 0
          sumtemp += 1.0*REWV*A[1,DESTC[j]]*prod(PROBC[routetemp[2:end]])*(1-PROBC[j])
        else
          sumtemp += 1.0*REWV*A[DESTC[i],DESTC[j]]*prod(PROBC[routetemp[findfirst(x -> x==i, routetemp)+1:end]])*(1-PROBC[i])*(1-PROBC[j])
        end 
      end
      if sumtemp<sumj
        minj =j
        sumj = sumtemp
      end 
    end
    push!(routetemp,minj)
  end
  #println(routetemp)
  
  #Break ordering of routetemp by average vehicle capacity to create initial routes
  pos=2
  while true
    avgcap = 0
    routetemp2=[]
    while true
      push!(routetemp2,routetemp[pos])
      #Calculate average capacity so far
      avgcap += 1-PROBC[routetemp[pos]]
      pos +=1
      if avgcap>=CAPV || pos == length(routetemp)+1 break end
    end
    #println(routetemp2)
    push!(route,routetemp2) 
    if pos == length(routetemp)+1 break end
  end
  #println(route)

  #To generate incumbent value for routes above, run first node with routes as integer values and insert all scenarios as row as needed.
  

  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = [] #zeros(length(route),length(DESTC))
  for i in 1:length(route)
    push!(B1,1.0*zeros(1:length(DESTC)))
  end
  global B2 = [] #zeros(length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(route)
    push!(B2,[])
    for j in 1:length(DESTC1) 
      push!(B2[i],[])
      for j1 in 1:length(DESTC1)
        push!(B2[i][j],0 )
      end
    end
  end

  global B3 = [] # zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B3,[])
    for j in 1:length(DESTC) 
      push!(B3[i],[])
      for j1 in 1:length(DESTC)
        push!(B3[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B3[i][j][j1],0 )
        end
      end
    end
  end
  global B4 = [] #zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B4,[])
    for j in 1:length(DESTC) 
      push!(B4[i],[])
      for j1 in 1:length(DESTC)
        push!(B4[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B4[i][j][j1],0 )
        end
      end
    end
  end
  global B5 = [] #zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  for i in 1:length(route)
    push!(B5,[])
    for j in 1:length(DESTC) 
      push!(B5[i],[])
      for j1 in 1:length(DESTC)
        push!(B5[i][j],[])
        for j2 in 1:length(DESTC)
          push!(B5[i][j][j1],0 )
        end
      end
    end
  end
  global B6 = [] #zeros(length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(route)
    push!(B6,[])
    for j in 1:length(DESTC1) 
      push!(B6[i],[])
      for j1 in 1:length(DESTC1)
        push!(B6[i][j],0)
      end
    end
  end
  global B7 = [] #zeros(length(route),length(DESTC))
  for i in 1:length(route)
    push!(B7,zeros(1:length(DESTC)))
  end 
  global B8 = [] #zeros(length(scenario),length(route),length(DESTC1),length(DESTC1))
  for i in 1:length(scenario)
    push!(B8,[])
    for j in 1:length(route)
      push!(B8[i],[]) 
      for j1 in 1:length(DESTC1) 
        push!(B8[i][j],[])
        for j2 in 1:length(DESTC1)
          push!(B8[i][j][j1],0 )
        end
      end
    end
  end
  
  for r in 1:length(route)
    for i in route[r]
      B7[r][i] = 1.0
      B1[r][i] += 1.0
    end
    B6[r][route[r][end]][length(DESTC1)] = 1.0
    B6[r][length(DESTC1)][route[r][1]] = 1.0
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        B2[r][length(DESTC1)][i] = 1.0
        B2[r][i][length(DESTC1)] = 1.0
        if i != route[r][end]
          B6[r][i][route[r][findfirst(x -> x==i, route[r])+1]] = 1.0
          #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
        end
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r][i][j] = 1.0
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r][i][j][k] = 1.0
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r][i][j][k] = 1.0
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r][i][j][k] = 1.0
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end
  for r in 1:length(route), w in 1:length(scenario),i in 1:length(DESTC),j in 1:length(DESTC)
    B8[w][r][i][j]  = (1- scenario[w][i])*(1 - scenario[w][j])*B2[r][i][j]*prod(scenario[w][k]*B3[r][i][j][k]+B4[r][i][j][k]+B5[r][i][j][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != j && k != i)
  end
  for r in 1:length(route), w in 1:length(scenario),i in 1:length(DESTC)
    B8[w][r][length(DESTC1)][i ]  = (1 - scenario[w][i])*prod(scenario[w][k]*B2[r][k][i]+B2[r][i][k] + 1-B7[r][k] for k in 1:length(DESTC) if k != i)  
    B8[w][r][i][length(DESTC1)] = ( 1- scenario[w][i] )*prod(scenario[w][k]*B2[r][i][k]+B2[r][k][i] + 1-B7[r][k] for k in 1:length(DESTC) if k != i) 
  end 
  #####################################################################
  #read(STDIN,Char)
  global firstpassnode1nonewscenario = false
  println("Start Processing nodes")

  routesused = collect(1:length(route))
  #Initialize branch-and-bound tree and queue
  #Queue = Vector{Int}()
  Queue = PriorityQueue()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,collect(1:length(scenario)),[],[],[],[],[],[],[],[])) 
  #push!(Queue,1)
  Queue[1]= 0
  #Start processing queue in a Deep First mode
  global heursol
  global globalincumbentvalue = +Inf
  global globalincumbentsolutiont = []
  global globalincumbentsolutionz = []
    
  masterx = 0
  master = []
  #always process last element of queue
  tic()
  while length(Queue)>0 
    #current=Queue[end]
    current= dequeue!(Queue)
    (master, masterx) =process(current) 
    if current == 1
      global ROOTSOL=masterx
      global ROOTTIME= TimeMain +TimePric+TimeSce
    end
    
    global NODECOUNT += 1
    if length(Queue) > 0 && globalincumbentvalue != +Inf
      a,b=peek(Queue)
      if (globalincumbentvalue-b)/globalincumbentvalue < 0.02
        break 
      end
    end
     
    if TimeMain +TimePric+TimeSce >= 36000
      break
    end
  end
    
  for i in 1:length(globalincumbentsolutionz)
    if globalincumbentsolutionz[i] > 0.0001 && globalincumbentsolutionz[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  #println("result= ",globalincumbentvalue,"dsol= ", globalincumbentsolutiond, "time = ", TimeMain +TimePric+TimeSce)
  global NROUTES = sum(globalincumbentsolutionz)
  #global PERHIGH = sum(globalincumbentsolutiont)/length(globalincumbentsolutiont)
  return globalincumbentvalue,globalincumbentsolutiont
end #end solveSTODROPCHEURfunction

