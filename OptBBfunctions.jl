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
      vardfixedsolutionzero::Vector{Int}
      vardfixedsolutionone::Vector{Int}
      dsolution::Vector{Float64}
      zsolution::Vector{Float64}
      varzfixedsolutionzero::Vector{Int}
      restr::Vector{ConstraintRef}
    end 

#11, from 10, but heuristic approximation: 1) routes continues as is, enumeration at all nodes  2) scenario should be generated at root, but this time will fix to initial set only 3) no branch on d: round variables, fix and run gain to terminate node
#10,from n9, now no variable y(r,w,i,j)!!
# 9, equal to 8 + mixed best first strategy
#DDUEXACTN8, equal to N3, using deterministic rounded capacity cut ; but not using mixed best bound strategy with basis reconstruction  and not using reduced cost variable fixing
#DDUEXACTN3: Same as N2, but now introduce col gen and sce gen to b&b
#N2 was same as N1, but now create our own B&B on variables first xij and then di. Extended formulation. Now decide compensation at master problem, list all routes and scenaros at once, no column generation or cut ...solve relaxed MILP...added B&B)
#correlated marginals, worst case approach, decision dependent probabilistic VRP; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL
#This time vehicle sohould satisfy capacity for all scenarios. Only routes within capacity are generated a priori and symetry eliminated (as in N1)
function solveSTODDUEXACTN11(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vardsolutionfixedzero = []
    vardsolutionfixedone = []
    varzsolutionfixedzero = []
    restractive = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vardsolutionfixedzero = [vardsolutionfixedzero;tree[node].vardfixedsolutionzero]
      vardsolutionfixedone = [vardsolutionfixedone;tree[node].vardfixedsolutionone]
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
    #println("vardsolutionfixedzero =",vardsolutionfixedzero)
    #println("vardsolutionfixedone =",vardsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    Gomorycutsrounds = 0
    #Now formulate problem (relaxed problem now)
    #println("comecar node ",current," com routes ", length(routesused), " e cenarios", length(cenariosused))
    z = Vector{Variable}(length(routesused))
    const4 = Vector{ConstraintRef}(length(cenariosused))
    const11 = Vector{ConstraintRef}()
    const12 = Vector{ConstraintRef}(length(scenarioall))

    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=7200))
    #Reserve first fixed indexes for continuous variables
    @variable(model,  s >= 0)
    @variable(model,  u[1:length(DESTC)]>=0)
   
    #println("varz fixed zero ",varzsolutionfixedzero)
    for r in 1:length(routesused)
      if in(routesused[r],varzsolutionfixedzero)
        z[r] = @variable(model, lowerbound = 0, upperbound =0, basename="z[$r]")
      else
        z[r] = @variable(model, lowerbound = 0, upperbound =1, basename="z[$r]")
      end
    end

    for i in 1:length(tree[current].zsolution)
       setvalue(z[i],tree[current].zsolution[i])  #dual feasible?
    end
  
    @variable(model,  x[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j]>=0)

    if inst_idf[end] != 'A' #probod2 = 0, so no d variabels
      @variable(model,  d[1:length(DESTC)]>=0)
      @variable(model,  t[1:length(DESTC)]>=0)

      for i in vardsolutionfixedzero
        setlowerbound(d[i],0)
        setupperbound(d[i],0) 
      end
      for i in vardsolutionfixedone
        setlowerbound(d[i],1)
        setupperbound(d[i],1)
      end

      for i in 1:length(DESTC)
        #setlowerbound(d[i],0)
        #setupperbound(d[i],0)
        if length(tree[current].dsolution) != 0
          setvalue(d[i],tree[current].dsolution[i])  #dual feasible?
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
  
    if inst_idf[end] != 'A'
      @objective(model, Min, s - sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) - sum(PROBOD2[i]*t[i] for i in 1:length(DESTC)) )

      @constraint(model, const1[i in 1:length(DESTC)], -t[i]  >=  -10^4*d[i]) #be careful with cplex integrality tolerance
      @constraint(model, const2[i in 1:length(DESTC)], -t[i] >= -u[i] )

      for w in 1:length(cenariosused)
        const4[w] = @constraint(model,  s - sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scenario[cenariosused[w]][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],findfirst(x -> x == scenario[cenariosused[w]], scenarioall),i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
      end

      @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)
      @constraint(model, const9[i in 1:length(DESTC)], d[i] <= 1)


      @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
 
      for sce in 1:length(scenarioall) #eliminate scenario all 0 and all 1
        if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
          total = 0
          for k in 1:length(DESTC)
            if scenarioall[sce][k]== 1
              total += 1
            end
          end
          const12[sce]=@constraint(model, sum( sum(x[i,j] for i in 1:length(DESTC)  if i != j ) for j in 1:length(DESTC1))  >= ceil(total/CAPV) )
        else
          const12[sce]=@constraint(model, 0==0)
        end
      end


      const11 = [const11;restractive]
    else
      @objective(model, Min, s - sum(PROBOD[i]*u[i] for i in 1:length(DESTC))  )


      for w in 1:length(cenariosused)
        const4[w] = @constraint(model,  s - sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC))  - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],findfirst(x -> x == scenario[cenariosused[w]], scenarioall),i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
      end

      @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)


      @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
 
      for sce in 1:length(scenarioall) #eliminate scenario all 0 and all 1
        if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
          total = 0
          for k in 1:length(DESTC)
            if scenarioall[sce][k]== 1
              total += 1
            end
          end
          const12[sce]=@constraint(model, sum( sum(x[i,j] for i in 1:length(DESTC)  if i != j ) for j in 1:length(DESTC1))  >= ceil(total/CAPV) )
        else
          const12[sce]=@constraint(model, 0==0)
        end
      end


      const11 = [const11;restractive]
  
    end 

    while true
      status = solve(model)
      result = getobjectivevalue(model)
      resultz = getvalue(z)

      if inst_idf[end] != 'A' && getcategory(d[1]) == :Bin
      global TimeSce += toq()
      tic()
      else
      global TimeMain += toq()
      tic()
      end
      #println("Solved Node ",current," problem ", result, "; Time is ",TimeMain, " and length Queue is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #println(model)
      #read(STDIN,Char)

      stop = false
      mincost = false
      if status != :Optimal
        #prune   
        #pos = findfirst(x -> x == current, Queue)
        #deleteat!(Queue,pos)
        #delete!(Queue, current) 
        #println("Infeasible...prune ", status)
        #read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolutiond
      global globalincumbentsolutionz
     

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
      #for i=1:length(DESTC1)
      #  println(xsol[i,:])
      #end
      #read(STDIN,Char)
      dsol = zeros(1:length(DESTC))
      if inst_idf[end] != 'A'
        dsol = getvalue(d)
        for i=1:length(DESTC)
          if dsol[i]<=0.00001 
            dsol[i] = 0
          elseif  dsol[i] >=0.99999
            dsol[i]= 1
          end
        end
      end

      bestr = 0
      #if inst_idf[end] == 'A' || getcategory(d[1]) != :Bin
      if current == 1 && (inst_idf[end] == 'A' || getcategory(d[1]) != :Bin)
      E1dual = getdual(const5)
      G1dual = getdual(const10)
      D1dual = getdual(const4)
      maxcol = -Inf
      for r in 1:length(route)
        if !in(r,routesused)
          total= sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)) - sum(B6[r,i,j]*G1dual[i,j] for i in 1:length(DESTC1),j in 1:length(DESTC1) if  i != j) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[r,findfirst(x -> x == scenario[cenariosused[w]], scenarioall),i,j]*D1dual[w] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for w in 1:length(cenariosused) ) >= 0.00001  
          if total >= 0.00001
            if maxcol < total
              maxcol = total
              bestr = r
            end
          end
        end #r in route
      end

      if bestr != 0
        #println("time to insert route ",bestr," in node", current)
        r = bestr
        push!(routesused,r)
        push!(tree[current].addedsolutions,r)
        global ROUTCOUNT += 1
        mincost = true
        #read(STDIN,Char)
        tes = [B1[r,i] for i in 1:length(DESTC)]
        tes1 = [const5[i] for i in 1:length(DESTC) ]
        #for w in 1:length(cenariosused)   
        for i in 1:length(DESTC1),j in 1:length(DESTC1)
          if i !=j 
            push!(tes, -B6[r,i,j])
            push!(tes1, const10[i,j])
         end
        end
        for w in 1:length(cenariosused)
          push!(tes, -sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[r,findfirst(x -> x == scenario[cenariosused[w]], scenarioall),i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) )
          push!(tes1, const4[w])
        end
        index = length(tree[current].addedsolutions)
        push!(z,@variable(model,  lowerbound = 0, upperbound =1, basename="z[$index]",objective = 0.0, inconstraints = tes1, coefficients = tes ))
        #setvalue(y[findfirst(x -> x == best_s, solutionsused)],1
      end
      global TimePric += toq()
      tic()
      end

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

        stop2=false
        ssol = getvalue(s)
        usol= getvalue(u)
        if (maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(dsol-floor.(dsol)) <= 0.00001)
          #new integer solution
          #println("Integer Solution Found")
          #read(STDIN,Char)
          ##############################################
          #Formulate Separation Problem and Solve 
          #Get First Stage results
          for i in 1:length(resultz) #check for error and display
            if resultz[i] > 0.0001 && resultz[i] < 0.99
              println("SOLUCAO NAO INTEIRA!! ",i, ",",resultz[i])
              read(STDIN,Char) 
            end
          end

        
          ##############################################
          #Verify if it can be new incubent
          if result < globalincumbentvalue
            #println("New Incumbent Found")
            globalincumbentvalue = result
            globalincumbentsolutiond = dsol 
            globalincumbentsolutionz = resultz
          end 
          #pos = findfirst(x -> x == current, Queue)
          #deleteat!(Queue,pos)
          #delete!(Queue, current)
          #println("Will break after integer solution")
          stop=true
          break
        end
  

        #First insert scenario

        #if length(cenariosused) <= 100 
        if current == 1 && length(cenariosused) <= 50
        if inst_idf[end] != 'A'             
          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0] <= 1) 
 
          @objective(model2, Min, ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scen[i]*dsol[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
          @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

          @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

          @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


          @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,length(DESTC1),i] <= (1-scen[i]) )
          @constraint(model2, const22bb[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,length(DESTC1),i] <= resultz[r] )
          @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] >0], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

          @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
          @constraint(model2, const23bb[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,i,length(DESTC1)] <= resultz[r] )
          @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

          @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= (1-scen[i]) )
          @constraint(model2, const21bb[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= resultz[r] )

          @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

          @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )

          #println("Will solve Separation Problem for integer solution")
        else
          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0] <= 1) 
 
          @objective(model2, Min, ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
          @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

          @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

          @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


          @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,length(DESTC1),i] <= (1-scen[i]) )
          @constraint(model2, const22bb[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,length(DESTC1),i] <= resultz[r] )
          @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] >0], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

          @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
          @constraint(model2, const23bb[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,i,length(DESTC1)] <= resultz[r] )
          @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

          @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= (1-scen[i]) )
          @constraint(model2, const21bb[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= resultz[r] )

          @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

          @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )

        end

           
        solve(model2)
        end  
 
        #println("solve model 2 ",getobjectivevalue(model2))
        #println(getvalue(scen))
        #read(STDIN,Char)
        #if length(cenariosused) <= 100  && getobjectivevalue(model2) <= -0.0001 #-0.05
        if current == 1 && length(cenariosused) <= 50 && getobjectivevalue(model2) <= -0.05
           #println(getvalue(scen), " get= ", getobjectivevalue(model2))
           #read(STDIN,Char)
           #Add new scenario and continue
           scenarionew = round.(getvalue(scen))
           pos = findfirst(x -> x == scenarionew, scenarioall)
           #println("pos= ", pos)
           #read(STDIN,Char)
           #println("Integer Solution not valid. New scenario: ", scenarionew)
           push!(scenario,scenarionew)
           #Update Current node information on scenarios used
           push!(tree[current].addedcenarios,size(scenario,1))
           push!(cenariosused,size(scenario,1))
           #Create new variables
  
           #Create new constraints
           if inst_idf[end] != 'A'
             push!(const4, @constraint(model,  s - sum(scenario[end][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scenario[end][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],pos,i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))
           else
             push!(const4, @constraint(model,  s - sum(scenario[end][i]*u[i] for i in 1:length(DESTC))  - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],pos,i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))

           end

           stop2=true
        end
        
        if !stop2
          zdual = getdual(z)
          for i in 1:length(zdual)
            if result + zdual[i] > globalincumbentvalue && !in(i,varzsolutionfixedzero)
              #println("node ",current," z solution for route ",routesused[i], " will not improve")
              #read(STDIN,Char)
              #if resultz[i] >= 0.1 && resultz[i] <= 0.9
              #  println("erro z[", i,"]")
              #  read(STDIN,Char)
              #end
              setlowerbound(z[i],0)
              setupperbound(z[i],0) 
              push!(varzsolutionfixedzero,routesused[i])
              push!(tree[current].varzfixedsolutionzero,routesused[i])
              #Should still eliminate constraints not needed here  
            end
          end
        end


          #Insert Gomory custs for variables "d" and "x"
        if !stop2 && Gomorycutsrounds == 0 && current == 1
           #stop2 = true
           Gomorycutsrounds += 1
           cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
           #First check on variables d
           rowv = []
           for i in 1:length(DESTC)
             #println("test d[",i,"]", dsol)
             cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
             if  dsol[i]-floor(dsol[i]) > 0.00001  && cbasis[linearindex(d[i])] == :Basic #is basic variable and is fractional, so insert gomory cut
               #println("d[",i,"] is eligible ", d[i] )
               row = 1.0*zeros(length(cbasis))
               #identify row where variable is basic
               for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
                 #print("test d[",i,"] in constraint ",k)
                 ci = model.internalModel.inner
                 ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
                 #print(",",row[linearindex(d[i])],",") 
                 if row[linearindex(d[i])] == 1
                   #println("found d[",i,"] in ", k) 
                   break   
                 end
               end
               #println("now test row")
               if row == 1.0*zeros(length(cbasis))
                 println("ERROR GOMORY CUT D SECTION!")
                 read(STDIN,Char)
               else #insert cut
                 #push!(const11,@constraint(model, sum( min( ( (row[k]-floor(row[k]))/(dsol[i]-floor(dsol[i]) ) ),( (1-row[k]+floor(row[k])) /(1- dsol[i]+floor(dsol[i])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(d[i]) && row[k] != 0 ) >= 1 )  )
                 for k in 1:length(row)
                   if k != linearindex(d[i]) && row[k] != 0
                     if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
                       row[k] =  (dsol[i]-floor(dsol[i]))*row[k]
                     elseif k <= 2*length(DESTC)+1  && row[k] > 0 
                       row[k] =  (1- dsol[i]+floor(dsol[i]) )*row[k]
                     elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= dsol[i]-floor(dsol[i]) 
                       row[k] =  (1- dsol[i]+floor(dsol[i]) )*(row[k]-floor(row[k]))
                     elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > dsol[i]-floor(dsol[i]) 
                       row[k] =  (dsol[i]-floor(dsol[i]))*(1-row[k]+floor(row[k]))
                     else
                       println("ERROR! NOT COVERED GOMMORY D")
                       read(STDIN,Char)
                     end
                     #row[k] = min( ( (row[k]-floor(row[k]))/(dsol[i]-floor(dsol[i]) ) ),( (1-row[k]+floor(row[k])) /(1- dsol[i]+floor(dsol[i])) ) )
                   end
                 end 
                 row[linearindex(d[i])]=0
                 push!(row,(dsol[i]-floor(dsol[i]) )*(1- dsol[i]+floor(dsol[i]) ))
                 push!(rowv,row)
                 #println("GOMORY CUT FOR VAR D!")
                 #read(STDIN,Char)
               end 
             end
           end


           #Now check on variables x
           rowv2 = []

           if inst_idf[end] != 'A'         
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
           end #inst_id 
            
           #Only at the end insert new constraints
           for v = 1:length(rowv) 
             #push!(const11,@constraint(model, sum( rowv[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv[v][k] != 0 ) >= rowv[v][end] )  )
             #stop2 = true 
             #push!(tree[current].restr,const11[end])
             #println(const11[end])
             #println(rowv[v])
             #read(STDIN,Char)
           end
           for v = 1:length(rowv2)
             #push!(const11,@constraint(model, sum( rowv2[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv2[v][k] != 0 ) >= 1 )  )
             #stop2 = true
             #push!(tree[current].restr,const11[end])
           end
        end
        
        
        
        if !stop2  #time for branching
        #println("time for branching")
        if maximum(xsol-floor.(xsol)) > 0.00001
          if length(vardsolutionfixedzero) != 0 || length(vardsolutionfixedone) != 0
            #println("d ja estava fixed e veio xij frac")
            #println(vardsolutionfixedzero)
            #println(vardsolutionfixedone)
            #read(STDIN,Char)
          end 
          a = maximum(xsol-floor.(xsol))
          f=0
          g=0
          for i in 1:length(DESTC1),j in 1:length(DESTC1)
            if i != j
              if xsol[i,j]-floor.(xsol[i,j])== a
                f = i
                g = j
                break
              end
            end
          end
          #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
          b =[f;g] 
          push!(tree,TreeNode(current,[],[],[],[],[b],[],[],dsol,resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result


          push!(tree,TreeNode(current,[],[],[],[b],[],[],[],dsol,resultz,[],const11))
          #current = length(tree)
          push!(tree[current].children,length(tree))
          ##push!(Queue,length(tree))
          Queue[length(tree)]=result

          ##pos = findfirst(x -> x == current, Queue)
          ##deleteat!(Queue,pos)
          ##delete!(Queue,current)

          ##current=Queue[end]
          stop = true
          ##Gomorycutsrounds = 0
          ##current=length(tree)

          #global NODECOUNT += 1
          #setlowerbound(x[f,g],0)
          #setupperbound(x[f,g],0)
          #push!(varxsolutionfixedzero,b) 
        elseif maximum(dsol-floor.(dsol)) > 0.00001
          #println("Will fix d in node ", current, "sol now is ",result)
          #read(STDIN,Char)
          #heuristic fix solution and solve again to force terminate as integer
          for i in 1:length(DESTC)
            setcategory(d[i], :Bin)
        
          end
          for i in 1:length(DESTC1),j in 1:length(DESTC1)
            if i != j
              setlowerbound(x[i,j],xsol[i,j])
              setupperbound(x[i,j],xsol[i,j])
            end
          end

          global TimeMain += toq()
          tic()
        
        else
          println("situation not covered")
          read(STDIN,Char)
        end
        end #stop2
      global TimeSce += toq()
          tic()
      end # !mincost
      stop && break
    end #end while true 'solve)
    return model, result 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate Prices, Routes and Scenarios to be used
  ###############################################
  #Create vectors for all possible scenarios
  global scenarioall = Vector{Vector{Int}}() 
  fillscenario(scenarioall,zeros(length(DESTC)),length(DESTC),1)
  println("All possible Scenario vectors created of size ", length(scenarioall) )
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenarioall)
    if sum(scenarioall[s]) <= CAPV
      groupe = []
      for c in 1:length(scenarioall[s])
        if scenarioall[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #################################################

  ###############################################
  #Re-Create vectors for nly initial scenarios - use col  cut generation afterwards
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  push!(scenario, zeros(1:length(DESTC)))   #Scenario 1 has to be all zeros. All customers present. Scenario 1 is used for calculations later. 
  for i = 1:length(DESTC)
    SAMPLETEST = zeros(1:length(DESTC))
    SAMPLETEST[i]=1
    push!(scenario, SAMPLETEST)
  end
  push!(scenario, ones(1:length(DESTC)))
  ################################################

  #scenario = scenarioall
  scenariosused = collect(1:length(scenario))

  #Bypass REWOD and define new prices to pay ODs
  DESTC1=union(DESTC,[1])
  global PRICEOD
  global PRICEOD2
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance.
  global inst_idf
  if inst_idf[end]== 'A'
    PRICEOD2 = 1.0*zeros(1:length(DESTC))
    PRICEOD = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  elseif inst_idf[end]== 'B'
    PRICEOD = 1.0*zeros(1:length(DESTC))
    PRICEOD2 = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD2[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD2[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD2[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  else  #if inst_idf[end]= "C" or others
    PRICEOD = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
    PRICEOD2 = PRICEOD
  end #end inst_id
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC1),length(DESTC1))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B6 = zeros(length(route),length(DESTC1),length(DESTC1)) 
  global B7 = zeros(length(route),length(scenarioall),length(DESTC1),length(DESTC1))  
  for r in 1:length(route)
    B6[r,route[r][end],length(DESTC1)] = 1
    B6[r,length(DESTC1),route[r][1]] = 1
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      B2[r,length(DESTC1),i] = 1
      B2[r,i,length(DESTC1)] = 1
      if findfirst(x -> x==i, route[r]) != 0
        if i != route[r][end]
          B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]] = 1
          #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
        end
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end
  for r in 1:length(route), w in 1:length(scenarioall),i in 1:length(DESTC),j in 1:length(DESTC)
    if i != j
        B7[r,w,i,j]  = (1- scenarioall[w][i])*(1 - scenarioall[w][j])*B2[r,i,j]*prod(scenarioall[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i)
    end
  end

  for r in 1:length(route), w in 1:length(scenarioall),i in 1:length(DESTC)
    B7[r,w,length(DESTC1),i ]  = (1 - scenarioall[w][i])*prod(scenarioall[w][k]*B2[r,k,i]+B2[r,i,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != i)  
    B7[r,w,i,length(DESTC1)] = ( 1- scenarioall[w][i] )*prod(scenarioall[w][k]*B2[r,i,k]+B2[r,k,i] + 1-B1[r,k] for k in 1:length(DESTC) if k != i) 
  end 
  #####################################################################
  #read(STDIN,Char)
 
  #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  routesused = []
  for r in 1:length(route)
    if length(route[r]) <= 1
       push!(routesused,r)
    end
  end
  #routesused = collect(1:length(route))
  #Initialize branch-and-bound tree and queue
  #Queue = Vector{Int}()
  Queue = PriorityQueue()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,scenariosused,[],[],[],[],[],[],[],[])) 
  #push!(Queue,1)
  Queue[1]= 0
  #Start processing queue in a Deep First mode
  global heursol
  global globalincumbentvalue = +Inf
  global globalincumbentsolutiond = []
  global globalincumbentsolutionz = []
    
  masterx = 0
  master = []
  #always process last element of queue
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
    if TimeMain +TimePric+TimeSce >= 10800
      break
    end
  end
 
  for i in 1:length(globalincumbentsolutionz)
    if globalincumbentsolutionz[i] > 0.0001 && globalincumbentsolutionz[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  println("result= ",globalincumbentvalue," dsol= ", globalincumbentsolutiond, " time= ",TimeMain +TimePric+TimeSce)
  global NROUTES = sum(globalincumbentsolutionz)
  global PERHIGH = sum(globalincumbentsolutiond)/length(globalincumbentsolutiond)
  return globalincumbentvalue,globalincumbentsolutiond
end #end N11 function




#10,from n9, now no variable y(r,w,i,j)!!
# 9, equal to 8 + mixed best first strategy
#DDUEXACTN8, equal to N3, using deterministic rounded capacity cut ; but not using mixed best bound strategy with basis reconstruction  and not using reduced cost variable fixing
#DDUEXACTN3: Same as N2, but now introduce col gen and sce gen to b&b
#N2 was same as N1, but now create our own B&B on variables first xij and then di. Extended formulation. Now decide compensation at master problem, list all routes and scenaros at once, no column generation or cut ...solve relaxed MILP...added B&B)
#correlated marginals, worst case approach, decision dependent probabilistic VRP; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL
#This time vehicle sohould satisfy capacity for all scenarios. Only routes within capacity are generated a priori and symetry eliminated (as in N1)
function solveSTODDUEXACTN10(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vardsolutionfixedzero = []
    vardsolutionfixedone = []
    varzsolutionfixedzero = []
    restractive = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vardsolutionfixedzero = [vardsolutionfixedzero;tree[node].vardfixedsolutionzero]
      vardsolutionfixedone = [vardsolutionfixedone;tree[node].vardfixedsolutionone]
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
    #println("vardsolutionfixedzero =",vardsolutionfixedzero)
    #println("vardsolutionfixedone =",vardsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    Gomorycutsrounds = 0
    #Now formulate problem (relaxed problem now)
    #println("comecar node ",current," com routes ", length(routesused), " e cenarios", length(cenariosused))
    #read(STDIN,Char) 
    z = Vector{Variable}(length(routesused))
    const4 = Vector{ConstraintRef}(length(cenariosused))
    const11 = Vector{ConstraintRef}()
    const12 = Vector{ConstraintRef}(length(scenarioall))

    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=7200))
    #Reserve first fixed indexes for continuous variables
    @variable(model,  s >= 0)
    @variable(model,  u[1:length(DESTC)]>=0)
    
   
    #println("varz fixed zero ",varzsolutionfixedzero)
    for r in 1:length(routesused)
      if in(routesused[r],varzsolutionfixedzero)
        z[r] = @variable(model, lowerbound = 0, upperbound =0, basename="z[$r]")
      else
        z[r] = @variable(model, lowerbound = 0, upperbound =1, basename="z[$r]")
      end
    end

    for i in 1:length(tree[current].zsolution)
       setvalue(z[i],tree[current].zsolution[i])  #dual feasible?
    end
  
    @variable(model,  x[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j]>=0)
   if inst_idf[end] != 'A' #probod2 = 0, so no d variabels
      @variable(model,  d[1:length(DESTC)]>=0)
      @variable(model,  t[1:length(DESTC)]>=0)

      for i in vardsolutionfixedzero
        setlowerbound(d[i],0)
        setupperbound(d[i],0) 
      end
      for i in vardsolutionfixedone
        setlowerbound(d[i],1)
        setupperbound(d[i],1)
      end

      for i in 1:length(DESTC)
        #setlowerbound(d[i],0)
        #setupperbound(d[i],0)
        if length(tree[current].dsolution) != 0
          setvalue(d[i],tree[current].dsolution[i])  #dual feasible?
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

    if inst_idf[end] != 'A'  
 
      @objective(model, Min, s - sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) - sum(PROBOD2[i]*t[i] for i in 1:length(DESTC)) )

      @constraint(model, const1[i in 1:length(DESTC)], -t[i]  >=  -10^4*d[i]) #be careful with cplex integrality tolerance
      @constraint(model, const2[i in 1:length(DESTC)], -t[i] >= -u[i] )

      for w in 1:length(cenariosused)
        const4[w] = @constraint(model,  s - sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scenario[cenariosused[w]][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],findfirst(x -> x == scenario[cenariosused[w]], scenarioall),i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
      end

      @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)
      @constraint(model, const9[i in 1:length(DESTC)], d[i] <= 1)


      @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
 
      for sce in 1:length(scenarioall) #eliminate scenario all 0 and all 1
        if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
          total = 0
          for k in 1:length(DESTC)
            if scenarioall[sce][k]== 1
              total += 1
            end
          end
          const12[sce]=@constraint(model, sum( sum(x[i,j] for i in 1:length(DESTC)  if i != j ) for j in 1:length(DESTC1))  >= ceil(total/CAPV) )
        else
          const12[sce]=@constraint(model, 0==0)
        end
      end


      const11 = [const11;restractive]
   else
      @objective(model, Min, s - sum(PROBOD[i]*u[i] for i in 1:length(DESTC))  )

      for w in 1:length(cenariosused)
        const4[w] = @constraint(model,  s - sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC))  - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],findfirst(x -> x == scenario[cenariosused[w]], scenarioall),i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
      end

      @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)
      


      @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
 
      for sce in 1:length(scenarioall) #eliminate scenario all 0 and all 1
        if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
          total = 0
          for k in 1:length(DESTC)
            if scenarioall[sce][k]== 1
              total += 1
            end
          end
          const12[sce]=@constraint(model, sum( sum(x[i,j] for i in 1:length(DESTC)  if i != j ) for j in 1:length(DESTC1))  >= ceil(total/CAPV) )
        else
          const12[sce]=@constraint(model, 0==0)
        end
      end


      const11 = [const11;restractive]
   end 

    while true
      status = solve(model)
      result = getobjectivevalue(model)
      resultz = getvalue(z)
      global TimeMain += toq()
      tic()
      #println("Solved Node ",current," problem ", result, "; Time is ",TimeMain, " and length Queue is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #println(model)
      #read(STDIN,Char)

      stop = false
      mincost = false
      if status != :Optimal
        #prune   
        #pos = findfirst(x -> x == current, Queue)
        #deleteat!(Queue,pos)
        #delete!(Queue, current) 
        #println("Infeasible...prune ", status)
        #read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolutiond
      global globalincumbentsolutionz
     

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
      #for i=1:length(DESTC1)
      #  println(xsol[i,:])
      #end
      #read(STDIN,Char)
      dsol = zeros(1:length(DESTC))
      if inst_idf[end] != 'A'
        dsol = getvalue(d)
        for i=1:length(DESTC)
          if dsol[i]<=0.00001 
            dsol[i] = 0
          elseif  dsol[i] >=0.99999
            dsol[i]= 1
          end
        end
      end

      bestr = 0
      if inst_idf[end] == 'A' || getcategory(d[1]) != :Bin
      E1dual = getdual(const5)
      G1dual = getdual(const10)
      D1dual = getdual(const4)
      maxcol = -Inf
      for r in 1:length(route)
        if !in(r,routesused)
          total= sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)) - sum(B6[r,i,j]*G1dual[i,j] for i in 1:length(DESTC1),j in 1:length(DESTC1) if  i != j) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[r,findfirst(x -> x == scenario[cenariosused[w]], scenarioall),i,j]*D1dual[w] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for w in 1:length(cenariosused) ) >= 0.00001  
          if total >= 0.00001
            if maxcol < total
              maxcol = total
              bestr = r
            end
          end
        end #r in route
      end
     end

      if bestr != 0
        #println("time to insert route ",bestr," in node", current)
        r = bestr
        push!(routesused,r)
        push!(tree[current].addedsolutions,r)
        global ROUTCOUNT += 1
        mincost = true
        #read(STDIN,Char)
        tes = [B1[r,i] for i in 1:length(DESTC)]
        tes1 = [const5[i] for i in 1:length(DESTC) ]
        #for w in 1:length(cenariosused)   
        for i in 1:length(DESTC1),j in 1:length(DESTC1)
          if i !=j 
            push!(tes, -B6[r,i,j])
            push!(tes1, const10[i,j])
         end
        end
        for w in 1:length(cenariosused)
          push!(tes, -sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[r,findfirst(x -> x == scenario[cenariosused[w]], scenarioall),i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) )
          push!(tes1, const4[w])
        end
        index = length(tree[current].addedsolutions)
        push!(z,@variable(model,  lowerbound = 0, upperbound =1, basename="z[$index]",objective = 0.0, inconstraints = tes1, coefficients = tes ))
        #setvalue(y[findfirst(x -> x == best_s, solutionsused)],1
      end
      global TimePric += toq()
      tic()


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

        stop2=false
        ssol = getvalue(s)
        usol= getvalue(u)
        if maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(dsol-floor.(dsol)) <= 0.00001
          #new integer solution
          #println("Integer Solution Found")
          ##############################################
          #Formulate Separation Problem and Solve 
          #Get First Stage results
          for i in 1:length(resultz) #check for error and display
            if resultz[i] > 0.0001 && resultz[i] < 0.99
              println("SOLUCAO NAO INTEIRA!! ",i, ",",resultz[i])
              read(STDIN,Char) 
            end
          end

          if inst_idf[end] != 'A'      
            model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

            @variable(model2, scen[1:length(DESTC)], Bin)
            @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1) 
 
            @objective(model2, Min, ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scen[i]*dsol[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
            @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

            @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

            @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


            @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i]) )
            @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

            @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
            @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

            @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i]) )

            @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

            @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0.001], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )
          else
            model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

            @variable(model2, scen[1:length(DESTC)], Bin)
            @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1) 
 
            @objective(model2, Min, ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC))  - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
            @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

            @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

            @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


            @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i]) )
            @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

            @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
            @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

            @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i]) )

            @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

            @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0.001], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )

          end




 
          #println("Will solve Separation Problem for integer solution")
           
          solve(model2)
          
 
          #println("solve model 2 ",getobjectivevalue(model2))
          #println(getvalue(scen))
          #read(STDIN,Char)
          if getobjectivevalue(model2) <= -0.0001 #-0.05
            #println("ADD NEW SCENARIO ",getvalue(scen), " get= ", getobjectivevalue(model2))
            #read(STDIN,Char)
            #Add new scenario and continue
            scenarionew = round.(getvalue(scen))
            pos = findfirst(x -> x == scenarionew, scenarioall)
            #println("pos= ", pos)
            #read(STDIN,Char)
            #println("Integer Solution not valid. New scenario: ", scenarionew)
            push!(scenario,scenarionew)
            #Update Current node information on scenarios used
            push!(tree[current].addedcenarios,size(scenario,1))
            push!(cenariosused,size(scenario,1))
            #Create new variables
  
            #Create new constraints
            if inst_idf[end] != 'A' 
              push!(const4, @constraint(model,  s - sum(scenario[end][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scenario[end][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],pos,i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))
            else
              push!(const4, @constraint(model,  s - sum(scenario[end][i]*u[i] for i in 1:length(DESTC))  - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],pos,i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))
            end

            stop2=true
          else # end getobj < +.05
            #println("No new scenario on it solution")
          #end verification of scenario insertion for probable integer solution
          ##############################################
          #Verify if it can be new incubent
            if result < globalincumbentvalue
              #println("New Incumbent Found")
              globalincumbentvalue = result
              globalincumbentsolutiond = dsol 
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
        else

          ########################START SEPARATION PROBLEM IF ROOT NODE
        if length(cenariosused) <= 30  && current == 1 
        if inst_idf[end] != 'A'             
          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0] <= 1) 
 
          @objective(model2, Min, ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scen[i]*dsol[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
          @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

          @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

          @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


          @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,length(DESTC1),i] <= (1-scen[i]) )
          @constraint(model2, const22bb[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,length(DESTC1),i] <= resultz[r] )
          @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] >0], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

          @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
          @constraint(model2, const23bb[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,i,length(DESTC1)] <= resultz[r] )
          @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

          @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= (1-scen[i]) )
          @constraint(model2, const21bb[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= resultz[r] )

          @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

          @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )

          #println("Will solve Separation Problem for integer solution")
        else
          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0] <= 1) 
 
          @objective(model2, Min, ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
          @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

          @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

          @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


          @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,length(DESTC1),i] <= (1-scen[i]) )
          @constraint(model2, const22bb[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,length(DESTC1),i] <= resultz[r] )
          @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] >0], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

          @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
          @constraint(model2, const23bb[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0], y2[r,i,length(DESTC1)] <= resultz[r] )
          @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

          @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= (1-scen[i]) )
          @constraint(model2, const21bb[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= resultz[r] )

          @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

          @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )

        end

           
        solve(model2)
 
        #println("solve model 2 ",getobjectivevalue(model2))
        #println(getvalue(scen))
        #read(STDIN,Char)
        if  getobjectivevalue(model2) <= -0.0001 #-0.05
           #println(getvalue(scen), " get= ", getobjectivevalue(model2))
           #read(STDIN,Char)
           #Add new scenario and continue
           scenarionew = round.(getvalue(scen))
           pos = findfirst(x -> x == scenarionew, scenarioall)
           #println("pos= ", pos)
           #read(STDIN,Char)
           #println("Integer Solution not valid. New scenario: ", scenarionew)
           push!(scenario,scenarionew)
           #Update Current node information on scenarios used
           push!(tree[current].addedcenarios,size(scenario,1))
           push!(cenariosused,size(scenario,1))
           #Create new variables
  
           #Create new constraints
           if inst_idf[end] != 'A'
             push!(const4, @constraint(model,  s - sum(scenario[end][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scenario[end][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],pos,i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))
           else
             push!(const4, @constraint(model,  s - sum(scenario[end][i]*u[i] for i in 1:length(DESTC))  - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*B7[routesused[r],pos,i,j]*z[r] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))

           end

           stop2=true
        end
        end #end current == 1




          ########################END SEPARATION PROBLEM IF ROOT NODE
          ##Test to verify if D cuts works
          #for r = 1:length(routesused)
          #  total = 0 
          #  for i in 1:length(DESTC)
          #     total += B1[routesused[r],i]*resultz[r]*(PROBOD[i] + PROBOD2[i]*dsol[i])
          #  end
          #  if total > CAPV
          #    println("Achei VI para route ",r,"/",routesused[r]," :",total," ; ",sum(B1[routesused[r],i]*PROBOD[i] for i in 1:length(DESTC) )," ; ",sum(B1[routesused[r],i]*PROBOD2[i] for i in 1:length(DESTC) ))
          #    read(STDIN,Char)
          #  else  
          #    println("NAO Achei VI para route ",r,"/",routesused[r]," :",total," ; ",sum(B1[routesused[r],i]*PROBOD[i] for i in 1:length(DESTC) )," ; ",sum(B1[routesused[r],i]*PROBOD2[i] for i in 1:length(DESTC) ))
          #  end
          #end

          #Section for Reduced Cost variable Fixing. 
          #ddual = getdual(d)
          #for i in 1:length(ddual)
          #  if result + ddual[i] > globalincumbentvalue && !in(i,vardsolutionfixedone) && !in(i,vardsolutionfixedzero) #&& dsol[i] <= 0.9
          #    println("node ",current," d solution for customer ",i, " will not improve")
          #    read(STDIN,Char)
          #    if dsol[i] >= 0.1 && dsol[i] <= 0.9
          #      println("erro d[", i,"]")
          #      read(STDIN,Char)
          #    end
          #    setlowerbound(d[i],0)
          #    setupperbound(d[i],0) 
          #    push!(vardsolutionfixedzero,i)
          #    push!(tree[current].vardfixedsolutionzero,i) 
          #  end
          #end
          #xdual = getdual(x)
          #for i in 1:length(DESTC1)
          #  for j in 1:length(DESTC1)
          #    if i != j
          #      if result + xdual[i,j] > globalincumbentvalue && !in([i;j],varxsolutionfixedone) && !in([i;j],varxsolutionfixedzero) #&& dsol[i] <= 0.9
          #        println("node ",current," x solution ",i,",",j, " will not improve")
          #        read(STDIN,Char)
          #        if xsol[i,j] >= 0.1 && xsol[i,j] <= 0.9
          #          println("erro x[",i," ,",j,"]")
          #          read(STDIN,Char)
          #        end
          #        setlowerbound(x[i,j],0)
          #        setupperbound(x[i,j],0) 
          #        push!(varxsolutionfixedzero,[i;j])
          #        push!(tree[current].varxfixedsolutionzero,[i;j]) 
          #      end
          #    end
          #  end
          #end
          zdual = getdual(z)
          for i in 1:length(zdual)
            if result + zdual[i] > globalincumbentvalue && !in(i,varzsolutionfixedzero)
              #println("node ",current," z solution for route ",routesused[i], " will not improve")
              #read(STDIN,Char)
              #if resultz[i] >= 0.1 && resultz[i] <= 0.9
              #  println("erro z[", i,"]")
              #  read(STDIN,Char)
              #end
              setlowerbound(z[i],0)
              setupperbound(z[i],0) 
              push!(varzsolutionfixedzero,routesused[i])
              push!(tree[current].varzfixedsolutionzero,routesused[i])
              #Should still eliminate constraints not needed here  
            end
          end
          #End Section for reduced Cost variable fixing

          ##insert rounded capacity cuts if possible (deterministic)
          ##if Gomorycutsrounds == 0
          #for sce in 1:length(scenarioall) #eliminate scenario all 0 and all 1
          #  if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
          #    sumtemp=0
          #    #println("test scenario= ", scenarioall[s])
          #    #read(STDIN,Char)
          #    for i=1:length(scenarioall[sce])
          #      if scenarioall[sce][i]== 1
          #        sumtemp += xsol[i,length(DESTC1)]
          #        #sumtemp += xsol[length(DESTC1),i]
          #      end
          #      for j=1:length(scenarioall[sce])
          #        if scenarioall[sce][i]== 1 && scenarioall[sce][j]== 0
          #           sumtemp += xsol[i,j]
          #        end
          #      end
          #    end
          #    #println("resumo ",sumtemp + 0.00001,",", ceil(sum(scenarioall[s])/CAPV))
          #    #read(STDIN,Char)
          #    if sumtemp + 0.00001 < ceil(sum(scenarioall[sce])/CAPV)
          #      #insert
          #      push!(const11, @constraint(model, sum(sum(x[i,j] for i=1:length(DESTC) if scenarioall[sce][i]== 1 && (j==length(DESTC1) || scenarioall[sce][j]== 0) ) for j=1:length(DESTC1)  ) >= ceil(sum(scenarioall[sce])/CAPV) ) )
          #      push!(tree[current].restr,const11[end])
          #      println("insert RCC and out on node ", current)
          #      println(const11[end])
          #      #read(STDIN,Char)
          #      stop2=true
          #      break
          #    end
          #  end
          #end
          ##end #Gomorycutsrounds==0

          #Insert Gomory custs for variables "d" and "x"
          if !stop2 && Gomorycutsrounds == 0 && current == 1
            #stop2 = true
            Gomorycutsrounds += 1
            cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
            #First check on variables d
            rowv = []
            for i in 1:length(DESTC)
              #println("test d[",i,"]", dsol)
              cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
              if  dsol[i]-floor(dsol[i]) > 0.00001 && cbasis[linearindex(d[i])] == :Basic #is basic variable and is fractional, so insert gomory cut
                #println("d[",i,"] is eligible ", d[i] )
                row = 1.0*zeros(length(cbasis))
                #identify row where variable is basic
                for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
                  #print("test d[",i,"] in constraint ",k)
                  ci = model.internalModel.inner
                  ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
                  #print(",",row[linearindex(d[i])],",") 
                  if row[linearindex(d[i])] == 1
                    #println("found d[",i,"] in ", k) 
                    break   
                  end
                end
                #println("now test row")
                if row == 1.0*zeros(length(cbasis))
                  println("ERROR GOMORY CUT D SECTION!")
                  read(STDIN,Char)
                else #insert cut
                  #push!(const11,@constraint(model, sum( min( ( (row[k]-floor(row[k]))/(dsol[i]-floor(dsol[i]) ) ),( (1-row[k]+floor(row[k])) /(1- dsol[i]+floor(dsol[i])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(d[i]) && row[k] != 0 ) >= 1 )  )
                  for k in 1:length(row)
                    if k != linearindex(d[i]) && row[k] != 0
                      if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
                        row[k] =  (dsol[i]-floor(dsol[i]))*row[k]
                      elseif k <= 2*length(DESTC)+1  && row[k] > 0 
                        row[k] =  (1- dsol[i]+floor(dsol[i]) )*row[k]
                      elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= dsol[i]-floor(dsol[i]) 
                        row[k] =  (1- dsol[i]+floor(dsol[i]) )*(row[k]-floor(row[k]))
                      elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > dsol[i]-floor(dsol[i]) 
                        row[k] =  (dsol[i]-floor(dsol[i]))*(1-row[k]+floor(row[k]))
                      else
                        println("ERROR! NOT COVERED GOMMORY D")
                        read(STDIN,Char)
                      end
                      #row[k] = min( ( (row[k]-floor(row[k]))/(dsol[i]-floor(dsol[i]) ) ),( (1-row[k]+floor(row[k])) /(1- dsol[i]+floor(dsol[i])) ) )
                    end
                  end 
                  row[linearindex(d[i])]=0
                  push!(row,(dsol[i]-floor(dsol[i]) )*(1- dsol[i]+floor(dsol[i]) ))
                  push!(rowv,row)
                  #println("GOMORY CUT FOR VAR D!")
                  #read(STDIN,Char)
                end 
              end
            end


            #Now check on variables x
            rowv2 = []

            if inst_idf[end] != 'A'        
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
            end #inst_idf
            
            #Only at the end insert new constraints
            for v = 1:length(rowv) 
              #push!(const11,@constraint(model, sum( rowv[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv[v][k] != 0 ) >= rowv[v][end] )  )
              #stop2 = true 
              #push!(tree[current].restr,const11[end])
              #println(const11[end])
              #println(rowv[v])
              #read(STDIN,Char)
            end
            for v = 1:length(rowv2)
              #push!(const11,@constraint(model, sum( rowv2[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv2[v][k] != 0 ) >= 1 )  )
              #stop2 = true
              #push!(tree[current].restr,const11[end])
            end
          end
          #End insert gomory custs        
        end
        
        
        if !stop2  #time for branching
        #println("time for branching")
        if maximum(xsol-floor.(xsol)) > 0.00001
          if length(vardsolutionfixedzero) != 0 || length(vardsolutionfixedone) != 0
            #println("d ja estava fixed e veio xij frac")
            #println(vardsolutionfixedzero)
            #println(vardsolutionfixedone)
            #read(STDIN,Char)
          end 
          a = maximum(xsol-floor.(xsol))
          f=0
          g=0
          for i in 1:length(DESTC1),j in 1:length(DESTC1)
            if i != j
              if xsol[i,j]-floor.(xsol[i,j])== a
                f = i
                g = j
                break
              end
            end
          end
          #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
          b =[f;g] 
          push!(tree,TreeNode(current,[],[],[],[],[b],[],[],dsol,resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result


          push!(tree,TreeNode(current,[],[],[],[b],[],[],[],dsol,resultz,[],const11))
          #current = length(tree)
          push!(tree[current].children,length(tree))
          ##push!(Queue,length(tree))
          Queue[length(tree)]=result

          ##pos = findfirst(x -> x == current, Queue)
          ##deleteat!(Queue,pos)
          ##delete!(Queue,current)

          ##current=Queue[end]
          stop = true
          ##Gomorycutsrounds = 0
          ##current=length(tree)

          #global NODECOUNT += 1
          #setlowerbound(x[f,g],0)
          #setupperbound(x[f,g],0)
          #push!(varxsolutionfixedzero,b) 
        elseif maximum(dsol-floor.(dsol)) > 0.00001
          a = maximum(dsol-floor.(dsol))
          f=0
          for i in 1:length(DESTC)
            if dsol[i]-floor.(dsol[i])== a
              f = i
              break
            end
          end
          #println("Sol d still frac, so split into two nodes using max frac", f)  
          push!(tree,TreeNode(current,[],[],[],[],[],[],[f],dsol,resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result
          push!(tree,TreeNode(current,[],[],[],[],[],[f],[],dsol,resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result

          #pos = findfirst(x -> x == current, Queue)
          #deleteat!(Queue,pos)
          #delete!(Queue,current)

          #current=Queue[end]

          stop = true
          #Gomorycutsrounds = 0
          #current=length(tree)

          #global NODECOUNT += 1
          #setlowerbound(d[f],0)
          #setupperbound(d[f],0)
          #push!(vardsolutionfixedzero,f) 
        else
          println("situation not covered")
          read(STDIN,Char)
        end
        end #stop2
      global TimeSce += toq()
          tic()
      end # !mincost
      stop && break
    end #end while true 'solve)
    return model, result 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate Prices, Routes and Scenarios to be used
  ###############################################
  #Create vectors for all possible scenarios
  global scenarioall = Vector{Vector{Int}}() 
  fillscenario(scenarioall,zeros(length(DESTC)),length(DESTC),1)
  println("All possible Scenario vectors created of size ", length(scenarioall) )
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenarioall)
    if sum(scenarioall[s]) <= CAPV
      groupe = []
      for c in 1:length(scenarioall[s])
        if scenarioall[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #################################################

  ###############################################
  #Re-Create vectors for nly initial scenarios - use col  cut generation afterwards
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  push!(scenario, zeros(1:length(DESTC)))   #Scenario 1 has to be all zeros. All customers present. Scenario 1 is used for calculations later. 
  for i = 1:length(DESTC)
    SAMPLETEST = zeros(1:length(DESTC))
    SAMPLETEST[i]=1
    push!(scenario, SAMPLETEST)
  end
  push!(scenario, ones(1:length(DESTC)))
  ################################################

  #scenario = scenarioall
  scenariosused = collect(1:length(scenario))

  #Bypass REWOD and define new prices to pay ODs
  DESTC1=union(DESTC,[1])
  global PRICEOD
  global PRICEOD2
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance.
  global inst_idf
  if inst_idf[end]== 'A'
    PRICEOD2 = 1.0*zeros(1:length(DESTC))
    PRICEOD = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  elseif inst_idf[end]== 'B'
    PRICEOD = 1.0*zeros(1:length(DESTC))
    PRICEOD2 = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD2[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD2[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD2[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  else  #if inst_id[end]= "C" or others
    PRICEOD = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
    PRICEOD2 = PRICEOD
  end #end inst_id
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC1),length(DESTC1))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B6 = zeros(length(route),length(DESTC1),length(DESTC1)) 
  global B7 = zeros(length(route),length(scenarioall),length(DESTC1),length(DESTC1))  
  for r in 1:length(route)
    B6[r,route[r][end],length(DESTC1)] = 1
    B6[r,length(DESTC1),route[r][1]] = 1
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      B2[r,length(DESTC1),i] = 1
      B2[r,i,length(DESTC1)] = 1
      if findfirst(x -> x==i, route[r]) != 0
        if i != route[r][end]
          B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]] = 1
          #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
        end
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end
  for r in 1:length(route), w in 1:length(scenarioall),i in 1:length(DESTC),j in 1:length(DESTC)
    if i != j
        B7[r,w,i,j]  = (1- scenarioall[w][i])*(1 - scenarioall[w][j])*B2[r,i,j]*prod(scenarioall[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i)
    end
  end

  for r in 1:length(route), w in 1:length(scenarioall),i in 1:length(DESTC)
    B7[r,w,length(DESTC1),i ]  = (1 - scenarioall[w][i])*prod(scenarioall[w][k]*B2[r,k,i]+B2[r,i,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != i)  
    B7[r,w,i,length(DESTC1)] = ( 1- scenarioall[w][i] )*prod(scenarioall[w][k]*B2[r,i,k]+B2[r,k,i] + 1-B1[r,k] for k in 1:length(DESTC) if k != i) 
  end 
  #####################################################################
  #read(STDIN,Char)
 
  #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  routesused = []
  for r in 1:length(route)
    if length(route[r]) <= 1
       push!(routesused,r)
    end
  end
  #routesused = collect(1:length(route))
  #Initialize branch-and-bound tree and queue
  #Queue = Vector{Int}()
  Queue = PriorityQueue()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,scenariosused,[],[],[],[],[],[],[],[])) 
  #push!(Queue,1)
  Queue[1]= 0
  #Start processing queue in a Deep First mode
  global heursol
  global globalincumbentvalue = +Inf
  global globalincumbentsolutiond = []
  global globalincumbentsolutionz = []
    
  masterx = 0
  master = []
  #always process last element of queue
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
    if TimeMain +TimePric+TimeSce >= 10800+3600
      break
    end
  end
 
  for i in 1:length(globalincumbentsolutionz)
    if globalincumbentsolutionz[i] > 0.0001 && globalincumbentsolutionz[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  println("result= ",globalincumbentvalue,"dsol= ", globalincumbentsolutiond, "time = ", TimeMain +TimePric+TimeSce)
  global NROUTES = sum(globalincumbentsolutionz)
  global PERHIGH = sum(globalincumbentsolutiond)/length(globalincumbentsolutiond)
  return globalincumbentvalue,globalincumbentsolutiond
end #end N10 function





# 9, equal to 8 + mixed best first strategy
#DDUEXACTN8, equal to N3, using deterministic rounded capacity cut ; but not using mixed best bound strategy with basis reconstruction  and not using reduced cost variable fixing
#DDUEXACTN3: Same as N2, but now introduce col gen and sce gen to b&b
#N2 was same as N1, but now create our own B&B on variables first xij and then di. Extended formulation. Now decide compensation at master problem, list all routes and scenaros at once, no column generation or cut ...solve relaxed MILP...added B&B)
#correlated marginals, worst case approach, decision dependent probabilistic VRP; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL
#This time vehicle sohould satisfy capacity for all scenarios. Only routes within capacity are generated a priori and symetry eliminated (as in N1)
function solveSTODDUEXACTN9(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vardsolutionfixedzero = []
    vardsolutionfixedone = []
    varzsolutionfixedzero = []
    restractive = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vardsolutionfixedzero = [vardsolutionfixedzero;tree[node].vardfixedsolutionzero]
      vardsolutionfixedone = [vardsolutionfixedone;tree[node].vardfixedsolutionone]
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
    #println("vardsolutionfixedzero =",vardsolutionfixedzero)
    #println("vardsolutionfixedone =",vardsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    Gomorycutsrounds = 0
    #Now formulate problem (relaxed problem now)
    #println("comecar node ",current," com routes ", length(routesused), " e cenarios", length(cenariosused))
    z = Vector{Variable}(length(routesused))
    y = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const4 = Vector{ConstraintRef}(length(cenariosused))
    const6 = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const7 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const8 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    #const10 = [[Vector{ConstraintRef}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)]
    const11 = Vector{ConstraintRef}()
    const12 = [Vector{ConstraintRef}(length(cenariosused))  for _ in 1:length(scenarioall) ]
    

    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=7200))
    #Reserve first fixed indexes for continuous variables
    @variable(model,  s >= 0)
    @variable(model,  u[1:length(DESTC)]>=0)
    @variable(model,  t[1:length(DESTC)]>=0)
    #@variable(model,  z[r in routesused]>=0)
    #@variable(model,  y[r in routesused,w in cenariosused,i in 1:length(DESTC1),j in 1:length(DESTC1); i != j] >=0)
   
    #println("varz fixed zero ",varzsolutionfixedzero)
    for r in 1:length(routesused)
      if in(routesused[r],varzsolutionfixedzero)
        z[r] = @variable(model, lowerbound = 0, upperbound =0, basename="z[$r]")
      else
        z[r] = @variable(model, lowerbound = 0, upperbound =1, basename="z[$r]")
      end
    end

    for i in 1:length(tree[current].zsolution)
       setvalue(z[i],tree[current].zsolution[i])  #dual feasible?
    end
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
      if i != j && !in(routesused[r],varzsolutionfixedzero) &&  (i==length(DESTC1) || scenario[cenariosused[w]][i] == 0) && (j==length(DESTC1) || scenario[cenariosused[w]][j] == 0) &&(B2[routesused[r],i,j] == 1) 
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1, basename="y[$r][$w][$i][$j]")
      else
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0, basename="y[$r][$w][$i][$j]")
      end
    end

     for sce in 1:length(scenarioall) #eliminate scenario all 0 and all 1
       if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
         for w in 1:length(cenariosused)
           total = 0
           for k in 1:length(DESTC)
             if scenarioall[sce][k]== 1 && scenario[cenariosused[w]][k]==0
               total += 1
             end
           end
           const12[sce][w]=@constraint(model, sum( sum( sum(y[r][w][i][j] for i in 1:length(DESTC)  if i != j && scenario[cenariosused[w]][i]==0 && scenarioall[sce][i]== 1 && (j==length(DESTC1) || scenarioall[sce][j]== 0) && (j==length(DESTC1) || scenario[cenariosused[w]][j]==0) ) for j in 1:length(DESTC1)) for r in 1:length(routesused) ) >= ceil(total/CAPV) )
         end
       else
         for w in 1:length(cenariosused)
           const12[sce][w]=@constraint(model, 0==0)
         end
       end
     end
 








    @variable(model,  x[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j]>=0)
    @variable(model,  d[1:length(DESTC)]>=0)

    for i in vardsolutionfixedzero
      setlowerbound(d[i],0)
      setupperbound(d[i],0) 
    end
    for i in vardsolutionfixedone
      setlowerbound(d[i],1)
      setupperbound(d[i],1)
    end

    for i in 1:length(DESTC)
      setlowerbound(d[i],0)
      setupperbound(d[i],0)
      if length(tree[current].dsolution) != 0
        #setvalue(d[i],tree[current].dsolution[i])  #dual feasible?
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
 
    @objective(model, Min, s - sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) - sum(PROBOD2[i]*t[i] for i in 1:length(DESTC)) )

    @constraint(model, const1[i in 1:length(DESTC)], -t[i]  >=  -10^4*d[i]) #be careful with cplex integrality tolerance
    @constraint(model, const2[i in 1:length(DESTC)], -t[i] >= -u[i] )

    for w in 1:length(cenariosused)
      const4[w] = @constraint(model,  s - sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[cenariosused[w]][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][w][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
    end

    @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)
    #@constraint(model, const5[i in 1:length(DESTC)], sum(x[i,j] for j in 1:length(DESTC1) if i != j) == 1)
    @constraint(model, const9[i in 1:length(DESTC)], d[i] <= 1)


    @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
    #println("varx fixed = ",varxsolutionfixedzero, " : ", varxsolutionfixedone)
    #for w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
    #  if i != j && (i==length(DESTC1) || scenario[cenariosused[w]][i]==0) && (j==length(DESTC1) || scenario[cenariosused[w]][j]==0) 
    #    if  in([i;j],varxsolutionfixedzero)
    #      const10[w][i][j]=@constraint(model, sum(y[r][w][i][j] for r in 1:length(routesused)) ==0)
    #      #const10[w][i][j]=@constraint(model, sum(y[r][w][i][j] for r in 1:length(routesused)) - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
    #    elseif  in([i;j],varxsolutionfixedone)
    #      const10[w][i][j]=@constraint(model, sum(y[r][w][i][j] for r in 1:length(routesused)) ==1)
    #      #const10[w][i][j]=@constraint(model, sum(y[r][w][i][j] for r in 1:length(routesused)) - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
    #    else
    #      const10[w][i][j]=@constraint(model, sum(y[r][w][i][j] for r in 1:length(routesused)) - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
    #    end
    #  else
    #    const10[w][i][j]= @constraint(model, 0 ==0)
    #  end 
    #end 
    
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
      if i != j && !in(routesused[r],varzsolutionfixedzero) &&   scenario[cenariosused[w]][i] == 0 &&  scenario[cenariosused[w]][j] == 0 &&(B2[routesused[r],i,j] == 1) 
        const6[r][w][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[cenariosused[w]][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
      else
        const6[r][w][i][j] = @constraint(model, 0 == 0)
      end
    end

    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC)
      if !in(routesused[r],varzsolutionfixedzero)  &&  scenario[cenariosused[w]][i] == 0  
        const7[r][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
        const8[r][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
      else
        const7[r][w][i] = @constraint(model, 0 == 0)
        const8[r][w][i] = @constraint(model, 0 == 0)
      end
    end
      #@constraint(model, const22, sum(Variable(model,i+1) for i in 1:length(routesused)) >= 20 )
    #Initial cuts constraints
    #println(const11)
    const11 = [const11;restractive] 

    while true
      status = solve(model)
      global TimeMain += toq()
      tic()
      println("Solved Node ",current," problem ", getobjectivevalue(model), "; Time is ",TimeMain, " and length Queue is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #println(model)
      #cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
      #println("d= ",cbasis[length(cbasis)-length(DESTC)+1:length(cbasis)])
      #read(STDIN,Char)
      
      #a= length(cbasis)-length(DESTC1)*(length(DESTC1)-1) +1
      #b=  length(cbasis)-length(DESTC)
      #println("x= ",cbasis[a:b])
      #read(STDIN,Char)
      #println("t e u= ",cbasis[(length(cbasis)-length(DESTC1)*(length(DESTC1)-1) -2*length(DESTC) +1):(length(cbasis)-length(DESTC1)*(length(DESTC1)-1))])
      #read(STDIN,Char)
      #mpb = model.internalModel  #MathprogBase Model
      #cpx = mpb.inner  #CPLEX Model
      #println(linearindex(d[6]))
      #println(Variable(model,2))
      #for k in 1:length(model.linconstr)
      #  #row = CPLEX.get_tableau_row(cpx,k-1)
      #  ci = model.internalModel.inner
      #  row = zeros(length(cbasis))
      #  ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k,row)
      #  read(STDIN,Char)
      #  if row[length(cbasis)] == 1
      #    println(row)
      #    println(k,",",model.linconstr[k])
      #    read(STDIN,Char)
      #  end
      #end
      #print(";",NODECOUNT,":",TimeMain,"...")
      result = getobjectivevalue(model)
      resultz = getvalue(z)
      #read(STDIN,Char)
      #println("d= ",getvalue(d))
      #read(STDIN,Char)
      #for r in 1:length(routesused)
      #  for w in 1:length(cenariosused)
      #    for i in 1:length(DESTC1)
      #      for j in 1:length(DESTC1)
      #        if i != j
      #          if getvalue(y[r][w][i][j])>0.0001
      #            #print("y [",r,",",w,",",i,",",j," ]:",getvalue(y[r,w,i,j]))
      #            #read(STDIN,Char)
      #          end
      #        end
      #      end
      #    end
      #  end 
      #end
      #println()
      #for r in 1:length(routesused)
      #  #println(r," ;",route[r])
      #  if getvalue(z[r]) >= 0.00001
      #    println(route[r],": z[ ",r," ]= ",getvalue(z[r])," ;") 
      #    read(STDIN,Char)
      #  end
      #end  #r 
      
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
      global globalincumbentsolution
     

      xsol1 = getvalue(x)
      #ysol1 = zeros(length(routesused),1,length(DESTC1),length(DESTC1))
      #for r in 1:length(routesused), i in 1:length(DESTC1), j in 1:length(DESTC1)
      #  if i != j
      #    ysol1[r,1,i,j] = getvalue(y[r][1][i][j])
      #  end
      #end
      xsol = zeros(length(DESTC1),length(DESTC1))

      for i=1:length(DESTC1)
        for j=1:length(DESTC1)
          if i != j 
            #xsol[i,j]= sum( ysol1[r,1,i,j] for r in 1:length(routesused))        
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
      #for i=1:length(DESTC1)
      #  for j=1:length(DESTC1)
      #    if i != j         
      #      if xsol1[i,j]<=0.00001 
      #        xsol[i,j] = 0
      #      elseif  xsol1[i,j] >=0.99999
      #        xsol[i,j]= 1
      #      else
      #        xsol[i,j]= xsol1[i,j]
      #        #print("xsol[",i," ,",j," ]",xsol[i,j])
      #      end
      #    end
      #  end
      #end
      #for i=1:length(DESTC1)
      #  println(xsol[i,:])
      #end
      #read(STDIN,Char)
      dsol = getvalue(d)
      for i=1:length(DESTC)
        if dsol[i]<=0.00001 
          dsol[i] = 0
        elseif  dsol[i] >=0.99999
          dsol[i]= 1
        end
      end

      
      E1dual = getdual(const5)
      G1dual = getdual(const10)
      D1dual = getdual(const4)
      #maxcol = 0
      #for r in 1:length(routesused)    
      #  if sum( B1[routesused[r],i]*E1dual[i] for i in 1:length(DESTC)) - sum(B2[routesused[r],i,j]*G1dual[i,j] for i in 1:length(DESTC),j in 1:length(DESTC) if  i != j)< -0.1
      #     setlowerbound(z[r],0)
      #     setupperbound(z[r],0) 
      #  end 
      #end

      maxcol = -Inf
      bestr = 0
      for r in 1:length(route)
        if !in(r,routesused)
          total= sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)) - sum(B6[r,i,j]*G1dual[i,j] for i in 1:length(DESTC1),j in 1:length(DESTC1) if  i != j)>= 0.00001
          #total = 0
          #for i in 1:length(DESTC)
          #  total += B1[r,i]*E1dual[i]
          #end
          #for w in 1:length(cenariosused), i in 1:length(DESTC1),j in 1:length(DESTC1)
          #  if  i != j && (i==length(DESTC1) || scenario[cenariosused[w]][i]==0) && (j==length(DESTC1) || scenario[cenariosused[w]][j]==0) && !in([i;j],varxsolutionfixedzero) && !in([i;j],varxsolutionfixedone)  
          #    total += -B6[r,i,j]*G1dual[w][i][j]
          #  end  
          #end

           #total = sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)) + sum(sum(B6[r,i,j]*G1dual[w][i][j] for i in 1:length(DESTC1),j in 1:length(DESTC1) if  i != j && (i==length(DESTC1) || scenario[cenariosused[w]][i]==0) && (j==length(DESTC1) || scenario[cenariosused[w]][j]==0) && !in([i;j],varxsolutionfixedzero) && !in([i;j],varxsolutionfixedone) ) for w in 1:length(cenariosused))
           if total >= 0.00001
             if maxcol < total
               maxcol = total
               bestr = r
             end
          end
        end #r in route
      end

      if bestr != 0
        #println("time to insert route ",bestr," in node", current)
        r = bestr
         push!(routesused,r)
            push!(tree[current].addedsolutions,r)
            global ROUTCOUNT += 1
            mincost = true
            #read(STDIN,Char)
            tes = [B1[r,i] for i in 1:length(DESTC)]
            tes1 = [const5[i] for i in 1:length(DESTC) ]
            #tes1 = []
            #tes = []
            #for w in 1:length(cenariosused)   
            #for i in 1:length(DESTC1),j in 1:length(DESTC1)
            #  if i !=j && (i==length(DESTC1) || scenario[cenariosused[w]][i]==0) && (j==length(DESTC1) || scenario[cenariosused[w]][j]==0) && !in([i;j],varxsolutionfixedzero) && !in([i;j],varxsolutionfixedone)
            #    push!(tes, -B6[r,i,j])
            #    push!(tes1, const10[w][i][j])
            #  end
            #end
            #end
            #println(tes)
            #println(tes1)
            index = length(tree[current].addedsolutions)
            push!(z,@variable(model,  lowerbound = 0, upperbound =1, basename="z[$index]",objective = 0.0, inconstraints = tes1, coefficients = tes ))
            #setvalue(y[findfirst(x -> x == best_s, solutionsused)],1)

            ynew = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:1]
 
            for  w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j &&  (i==length(DESTC1) || scenario[cenariosused[w]][i] == 0) && (j==length(DESTC1) || scenario[cenariosused[w]][j] == 0) &&(B2[r,i,j] == 1) 
                tes1 =[const4[w]]
                tes = [-REWV*A[DESTC1[i],DESTC1[j]]]
                for sce in 1:length(scenarioall) 
                  if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce])) && i!=length(DESTC1) 
                    if scenario[cenariosused[w]][i]==0 &&  scenarioall[sce][i]== 1 && (j==length(DESTC1) || scenarioall[sce][j]== 0) && (j==length(DESTC1) || scenario[cenariosused[w]][j]==0)
                      #push!(tes1,const12[sce][w])
                      #push!(tes,1.0)
                    end
                  end
                end
             
                #push!(tes1,const10[w][i][j])
                #push!(tes,1.0)
               
                ynew[1][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1, basename="y[$index][$w][$i][$j]", objective = 0.0, inconstraints = tes1, coefficients = tes)
              else
                ynew[1][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0, basename="y[$index][$w][$i][$j]")
              end
            end
            y = [y;ynew]




            const6new = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:1]
            const7new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:1]
            const8new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:1]
            for  w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
              if i != j && !in(r,varzsolutionfixedzero) &&   scenario[cenariosused[w]][i] == 0 &&  scenario[cenariosused[w]][j] == 0 &&(B2[r,i,j] == 1) 
                const6new[1][w][i][j] = @constraint(model, -(B2[r,i,j] + sum(scenario[cenariosused[w]][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*z[end] + y[end][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
              else
                const6new[1][w][i][j] = @constraint(model, 0 == 0)
              end
            end
            for w in 1:length(cenariosused),i in 1:length(DESTC)
              if !in(r,varzsolutionfixedzero)  &&  scenario[cenariosused[w]][i] == 0 
              const7new[1][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[r,k,i]+0+B2[r,i,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[end] + y[end][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
              const8new[1][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[end] + y[end][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
              else
                const7new[1][w][i] = @constraint(model, 0==0)
                const8new[1][w][i] = @constraint(model, 0==0)
              end
            end
            
            const6 = [const6;const6new]
            const7 = [const7;const7new]
            const8 = [const8;const8new]

      end
          global TimePric += toq()
            tic()


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

        stop2=false
        ssol = getvalue(s)
        usol= getvalue(u)
        if maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(dsol-floor.(dsol)) <= 0.00001
          #new integer solution
          #println("Integer Solution Found")
          ##############################################
          #Formulate Separation Problem and Solve 
          #Get First Stage results
          for i in 1:length(resultz) #check for error and display
            if resultz[i] > 0.0001 && resultz[i] < 0.99
              println("SOLUCAO NAO INTEIRA!! ",i, ",",resultz[i])
              read(STDIN,Char) 
            end
          end
      
          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1) 
 
          @objective(model2, Min, ssol - sum(scen[i]*usol[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scen[i]*dsol[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
          @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

          @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

          @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


          @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i]) )
          @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

          @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
          @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

          @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i]) )

          @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

          @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0.001], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )




 
          #println("Will solve Separation Problem for integer solution")
           
          solve(model2)

 
          #println("solve model 2",getobjectivevalue(model2))
          if getobjectivevalue(model2) <= -0.0001 #-0.05

            ##free all reduced cost fixed variables
            #for i in 1:length(DESTC)
            #  if !in(i,vardsolutionfixedzero) && !in(i,vardsolutionfixedone)
            #    setlowerbound(d[i],0)
            #    setupperbound(d[i],1)  
            #  end
            #end

            #println(getvalue(scen), " get= ", getobjectivevalue(model2))
            #read(STDIN,Char)
            #Add new scenario and continue
            scenarionew = round.(getvalue(scen))
            #println("Integer Solution not valid. New scenario: ", scenarionew)
            push!(scenario,scenarionew)
            #Update Current node information on scenarios used
            push!(tree[current].addedcenarios,size(scenario,1))
            push!(cenariosused,size(scenario,1))
            index = length(tree[current].addedcenarios)
            #Create new variables
            ynew = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:1] for _ in 1:length(routesused)]

            for  r in 1:length(routesused),i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j && !in(routesused[r],varzsolutionfixedzero) &&(i==length(DESTC1) || scenario[cenariosused[end]][i] == 0) && (j==length(DESTC1) || scenario[cenariosused[end]][j] == 0) &&(B2[routesused[r],i,j] == 1) 
                ynew[r][1][i][j] = @variable(model,  lowerbound = 0, upperbound =1, basename="y[$r][$index][$i][$j]")
              else
                ynew[r][1][i][j] = @variable(model,  lowerbound = 0, upperbound =0, basename="y[$r][$index][$i][$j]")
              end
            end
            for r in 1:length(routesused)
             y[r] = [y[r];ynew[r]]
            end 
            #Create new constraints
            push!(const4, @constraint(model,  s - sum(scenario[end][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[end][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][end][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))

            for sce = 1:length(scenarioall)
              if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
              total = 0
              for k in 1:length(DESTC)
                if scenarioall[sce][k]== 1 && scenario[cenariosused[end]][k]==0
                  total += 1
                end
              end            
              #push!(const12[sce],@constraint(model, sum( sum( sum(y[r][end][i][j] for i in 1:length(DESTC)  if scenario[cenariosused[end]][i]==0 && scenarioall[sce][i]== 1 && (j==length(DESTC1) || (scenarioall[sce][j]== 0 && scenario[cenariosused[end]][j]==0))) for j in 1:length(DESTC1)) for r in 1:length(routesused) ) >= ceil(total/CAPV) ))
              end
            end
 
            #const10new = [[Vector{ConstraintRef}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:1]
            #for i in 1:length(DESTC1),j in 1:length(DESTC1)
            #  if i != j
            #    if (i==length(DESTC1) || scenario[end][i]==0) && (j==length(DESTC1) || scenario[end][j]==0) && in([i;j],varxsolutionfixedzero)
            #      const10new[1][i][j]=@constraint(model, sum(y[r][end][i][j] for r in 1:length(routesused)) ==0)
            #    elseif (i==length(DESTC1) || scenario[end][i]==0) && (j==length(DESTC1) || scenario[end][j]==0) && in([i;j],varxsolutionfixedone)
            #      const10new[1][i][j]=@constraint(model, sum(y[r][end][i][j] for r in 1:length(routesused)) ==1)
            #    elseif (i==length(DESTC1) || scenario[end][i]==0) && (j==length(DESTC1) || scenario[end][j]==0) && !in([i;j],varxsolutionfixedzero) && !in([i;j],varxsolutionfixedone)
            #      const10new[1][i][j]=@constraint(model, sum(y[r][end][i][j] for r in 1:length(routesused)) - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
            #    else
            #      const10new[1][i][j]= @constraint(model, 0 ==0)
            #    end
            #  else
            #    const10new[1][i][j]= @constraint(model, 0 ==0)
            #  end 
            #end 
            #const10=[const10;const10new]

            
            const6new = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:1] for _ in 1:length(routesused)]
            const7new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:1] for _ in 1:length(routesused)]
            const8new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:1] for _ in 1:length(routesused)]
            for  r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC)
              if i != j && !in(routesused[r],varzsolutionfixedzero) &&scenario[end][i] == 0 && scenario[end][j] == 0 &&(B2[routesused[r],i,j] == 1) 
                const6new[r][1][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[end][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][end][i][j] >= 1- scenario[end][i] + 1 - scenario[end][j] - length(DESTC))
              else
                const6new[r][1][i][j] = @constraint(model, 0 == 0)
              end
            end
            for r in 1:length(routesused),i in 1:length(DESTC)
              if !in(routesused[r],varzsolutionfixedzero)  &&  scenario[end][i] == 0 
              const7new[r][1][i] = @constraint(model,  -(0 + sum(scenario[end][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][end][length(DESTC1)][i] >=  1 - scenario[end][i] - length(DESTC) + 1)
              const8new[r][1][i]= @constraint(model, -(0 + sum(scenario[end][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][end][i][length(DESTC1)] >= 1- scenario[end][i]  - length(DESTC) + 1)
              else
                const7new[r][1][i] = @constraint(model, 0==0)
                const8new[r][1][i] = @constraint(model, 0==0)
              end
            end
            for r in 1:length(routesused)
              const6[r]=[const6[r];const6new[r]]
              const7[r]=[const7[r];const7new[r]]
              const8[r]=[const8[r];const8new[r]]
            end
            stop2=true
          else # end getobj < +.05
          #end verification of scenario insertion for probable integer solution
          ##############################################
          #Verify if it can be new incubent
            if result < globalincumbentvalue
              #println("New Incumbent Found")
              globalincumbentvalue = result
              globalincumbentsolution = dsol #resultz
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
        else
          #Section for Reduced Cost variable Fixing. 
          #ddual = getdual(d)
          #for i in 1:length(ddual)
          #  if result + ddual[i] > globalincumbentvalue && !in(i,vardsolutionfixedone) && !in(i,vardsolutionfixedzero) #&& dsol[i] <= 0.9
          #    println("node ",current," d solution for customer ",i, " will not improve")
          #    read(STDIN,Char)
          #    if dsol[i] >= 0.1 && dsol[i] <= 0.9
          #      println("erro d[", i,"]")
          #      read(STDIN,Char)
          #    end
          #    setlowerbound(d[i],0)
          #    setupperbound(d[i],0) 
          #    push!(vardsolutionfixedzero,i)
          #    push!(tree[current].vardfixedsolutionzero,i) 
          #  end
          #end
          #xdual = getdual(x)
          #for i in 1:length(DESTC1)
          #  for j in 1:length(DESTC1)
          #    if i != j
          #      if result + xdual[i,j] > globalincumbentvalue && !in([i;j],varxsolutionfixedone) && !in([i;j],varxsolutionfixedzero) #&& dsol[i] <= 0.9
          #        println("node ",current," x solution ",i,",",j, " will not improve")
          #        read(STDIN,Char)
          #        if xsol[i,j] >= 0.1 && xsol[i,j] <= 0.9
          #          println("erro x[",i," ,",j,"]")
          #          read(STDIN,Char)
          #        end
          #        setlowerbound(x[i,j],0)
          #        setupperbound(x[i,j],0) 
          #        push!(varxsolutionfixedzero,[i;j])
          #        push!(tree[current].varxfixedsolutionzero,[i;j]) 
          #      end
          #    end
          #  end
          #end
          zdual = getdual(z)
          for i in 1:length(zdual)
            if result + zdual[i] > globalincumbentvalue && !in(i,varzsolutionfixedzero)
              println("node ",current," z solution for route ",routesused[i], " will not improve")
              #read(STDIN,Char)
              #if resultz[i] >= 0.1 && resultz[i] <= 0.9
              #  println("erro z[", i,"]")
              #  read(STDIN,Char)
              #end
              setlowerbound(z[i],0)
              setupperbound(z[i],0) 
              push!(varzsolutionfixedzero,routesused[i])
              push!(tree[current].varzfixedsolutionzero,routesused[i])
              for w in 1:length(cenariosused)
                for k in 1:length(DESTC1)
                  for j in 1:length(DESTC1)
                    if k != j
                      setlowerbound(y[i][w][k][j],0)
                      setupperbound(y[i][w][k][j],0)
                    end
                  end
                end
              end 
              #Should still eliminate constraints not needed here  
            end
          end
          #End Section for reduced Cost variable fixing

          ##insert rounded capacity cuts if possible (deterministic)
          ##if Gomorycutsrounds == 0
          #for sce in 1:length(scenarioall) #eliminate scenario all 0 and all 1
          #  if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
          #    sumtemp=0
          #    #println("test scenario= ", scenarioall[s])
          #    #read(STDIN,Char)
          #    for i=1:length(scenarioall[sce])
          #      if scenarioall[sce][i]== 1
          #        sumtemp += xsol[i,length(DESTC1)]
          #        #sumtemp += xsol[length(DESTC1),i]
          #      end
          #      for j=1:length(scenarioall[sce])
          #        if scenarioall[sce][i]== 1 && scenarioall[sce][j]== 0
          #           sumtemp += xsol[i,j]
          #        end
          #      end
          #    end
          #    #println("resumo ",sumtemp + 0.00001,",", ceil(sum(scenarioall[s])/CAPV))
          #    #read(STDIN,Char)
          #    if sumtemp + 0.00001 < ceil(sum(scenarioall[sce])/CAPV)
          #      #insert
          #      push!(const11, @constraint(model, sum(sum(x[i,j] for i=1:length(DESTC) if scenarioall[sce][i]== 1 && (j==length(DESTC1) || scenarioall[sce][j]== 0) ) for j=1:length(DESTC1)  ) >= ceil(sum(scenarioall[sce])/CAPV) ) )
          #      push!(tree[current].restr,const11[end])
          #      println("insert RCC and out on node ", current)
          #      println(const11[end])
          #      #read(STDIN,Char)
          #      stop2=true
          #      break
          #    end
          #  end
          #end
          ##end #Gomorycutsrounds==0

          #Insert Gomory custs for variables "d" and "x"
          if !stop2 && Gomorycutsrounds == 0
            #stop2 = true
            Gomorycutsrounds += 1
            cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
            #First check on variables d
            rowv = []
            for i in 1:length(DESTC)
              #println("test d[",i,"]", dsol)
              cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
              if cbasis[linearindex(d[i])] == :Basic && dsol[i]-floor(dsol[i]) > 0.00001  #is basic variable and is fractional, so insert gomory cut
                #println("d[",i,"] is eligible ", d[i] )
                row = 1.0*zeros(length(cbasis))
                #identify row where variable is basic
                for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
                  #print("test d[",i,"] in constraint ",k)
                  ci = model.internalModel.inner
                  ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
                  #print(",",row[linearindex(d[i])],",") 
                  if row[linearindex(d[i])] == 1
                    #println("found d[",i,"] in ", k) 
                    break   
                  end
                end
                #println("now test row")
                if row == 1.0*zeros(length(cbasis))
                  println("ERROR GOMORY CUT D SECTION!")
                  read(STDIN,Char)
                else #insert cut
                  #push!(const11,@constraint(model, sum( min( ( (row[k]-floor(row[k]))/(dsol[i]-floor(dsol[i]) ) ),( (1-row[k]+floor(row[k])) /(1- dsol[i]+floor(dsol[i])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(d[i]) && row[k] != 0 ) >= 1 )  )
                  for k in 1:length(row)
                    if k != linearindex(d[i]) && row[k] != 0
                      if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
                        row[k] =  (dsol[i]-floor(dsol[i]))*row[k]
                      elseif k <= 2*length(DESTC)+1  && row[k] > 0 
                        row[k] =  (1- dsol[i]+floor(dsol[i]) )*row[k]
                      elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= dsol[i]-floor(dsol[i]) 
                        row[k] =  (1- dsol[i]+floor(dsol[i]) )*(row[k]-floor(row[k]))
                      elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > dsol[i]-floor(dsol[i]) 
                        row[k] =  (dsol[i]-floor(dsol[i]))*(1-row[k]+floor(row[k]))
                      else
                        println("ERROR! NOT COVERED GOMMORY D")
                        read(STDIN,Char)
                      end
                      #row[k] = min( ( (row[k]-floor(row[k]))/(dsol[i]-floor(dsol[i]) ) ),( (1-row[k]+floor(row[k])) /(1- dsol[i]+floor(dsol[i])) ) )
                    end
                  end 
                  row[linearindex(d[i])]=0
                  push!(row,(dsol[i]-floor(dsol[i]) )*(1- dsol[i]+floor(dsol[i]) ))
                  push!(rowv,row)
                  println("GOMORY CUT FOR VAR D!")
                  #read(STDIN,Char)
                end 
              end
            end


            ##Now check on variables x
            #rowv2 = []
            ##if current==1
            #for i in 1:length(DESTC1)
            #  for j in 1:length(DESTC1)
            #    if i != j
            #      cbasis, rbasis = MathProgBase.getbasis(internalmodel(model))
            #      #println("test x[",i," ,",j,"]") 
            #      if cbasis[linearindex(x[i,j])] == :Basic && xsol[i,j]-floor(xsol[i,j]) > 0.00001  #is basic variable and is fractional, so insert gomory cut
            #        #println("x[",i," ,",j,"] is elegible ", cbasis[linearindex(x[i,j])], ",",xsol[i,j],",",in([i;j],varxsolutionfixedzero)," ,",in([i;j],varxsolutionfixedone)) 
            #        #read(STDIN,Char)
            #        row = 1.0*zeros(length(cbasis))
            #        #identify row where variable is basic
            #        for k in 1:length(rbasis) #length(model.linconstr) #try all constraints
            #          #print("test x[",i," ,",j,"] in constraint ", k)
            #          ci = model.internalModel.inner
            #          ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, k-1,row) #use CPLEX C function library
            #          #println("(",linearindex(x[i,j]),",",row[linearindex(x[i,j])],",")
            #          #if  row[linearindex(x[i,j])] != 0  read(STDIN,Char) end 
            #          if row[linearindex(x[i,j])] == 1
            #            #println("found x[",i," ,",j,"] in ", k) 
            #            break   
            #          end
            #        end
            #        if row == 1.0*zeros(length(cbasis))
            #          println("ERROR GOMORY CUT X SECTION!")
            #          read(STDIN,Char)
            #        else #insert cut
            #          #push!(const11,@constraint(model, sum(min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )*Variable(model,k) for k=1:length(cbasis) if k != linearindex(x[i,j]) && row[k] != 0 ) >= 1) )
            #          for k in 1:length(row)
            #            if k != linearindex(x[i,j]) && row[k] != 0
            #              if k <= 2*length(DESTC)+1  && row[k] < 0                   #First index positions to continous variables
            #                row[k] =  row[k]/(1- xsol[i,j]+floor(xsol[i,j]) )
            #              elseif k <= 2*length(DESTC)+1  && row[k] > 0 
            #                row[k] =  row[k]/(xsol[i,j]-floor(xsol[i,j]) )
            #              elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) <= xsol[i,j]-floor(xsol[i,j]) 
            #                row[k] =  (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) )
            #              elseif k > 2*length(DESTC)+1  && row[k]-floor(row[k]) > xsol[i,j]-floor(xsol[i,j]) 
            #                row[k] =  (1-row[k]+floor(row[k]))/(1-xsol[i,j]+floor(xsol[i,j]) )
            #              else
            #               println("ERROR! NOT COVERED GOMMORY X")
            #                read(STDIN,Char)
            #              end
            #              #row[k] = min( ( (row[k]-floor(row[k]))/(xsol[i,j]-floor(xsol[i,j]) ) ),( (1-row[k]+floor(row[k])) /(1- xsol[i,j]+floor(xsol[i,j])) ) )
            #            end
            #          end 
            #          row[linearindex(x[i,j])]=0
            #          #push!(row,xsol[i,j]-floor(xsol[i,j]) )
            #          push!(rowv2,row)
            #          println("GOMORY CUT FOR VAR X!")
            #          #read(STDIN,Char)
            #        end 
            #      end
            #    end
            #  end
            #end
            #end #current=1
            #Only at the end insert new constraints
            for v = 1:length(rowv) 
              push!(const11,@constraint(model, sum( rowv[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv[v][k] != 0 ) >= rowv[v][end] )  )
              if current == 1
                stop2 = true 
                push!(tree[current].restr,const11[end])
              end
              #println(const11[end])
              #println(rowv[v])
              #read(STDIN,Char)
            end
            #for v = 1:length(rowv2)
            #  push!(const11,@constraint(model, sum( rowv2[v][k]*Variable(model,k) for k=1:length(cbasis) if  rowv2[v][k] != 0 ) >= 1 )  )
            #  if current == 1
            #    push!(tree[current].restr,const11[end])
            #  end
            #end
          end
          #End insert gomory custs        
        end
        
        
        if !stop2  #time for branching
         #println("time for branching")
        #if maximum(xsol-floor.(xsol)) <= 0.00001
        #  #Use solution to define binary d
        #  for i in 1:length(DESTC)
        #    if -PROBOD2[i]*usol[i] > sum( D1dual[w]*PRICEOD[i]*scenario[cenariosused[w]][i] for w in 1:length(cenariosused))
        #       if in(i,vardsolutionfixedzero)
        #          println("d[ ",i,"] deveria ser 1 mas esta fixo em zero")
        #          read(STDIN,Char)
        #       else
        #         setlowerbound(d[i],1)
        #         setupperbound(d[i],1)
        #       end
        #     else
        #       if in(i,vardsolutionfixedone)
        #          println("d[ ",i,"] deveria ser 0 mas esta fixo em um")
        #          read(STDIN,Char)
        #       else
        #         setlowerbound(d[i],0)
        #         setupperbound(d[i],0)
        #       end
        #
        #     end
        #  end

        if maximum(xsol-floor.(xsol)) > 0.00001
          if length(vardsolutionfixedzero) != 0 || length(vardsolutionfixedone) != 0
            println("d ja estava fixed e veio xij frac")
            println(vardsolutionfixedzero)
            println(vardsolutionfixedone)
            #read(STDIN,Char)
          end 
          a = maximum(xsol-floor.(xsol))
          f=0
          g=0
          for i in 1:length(DESTC1),j in 1:length(DESTC1)
            if i != j
              if xsol[i,j]-floor.(xsol[i,j])== a
                f = i
                g = j
                break
              end
            end
          end
          #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
          b =[f;g] 
          push!(tree,TreeNode(current,[],[],[],[],[b],[],[],dsol,resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result


          push!(tree,TreeNode(current,[],[],[],[b],[],[],[],dsol,resultz,[],const11))
          #current = length(tree)
          push!(tree[current].children,length(tree))
          ##push!(Queue,length(tree))
          Queue[length(tree)]=result

          ##pos = findfirst(x -> x == current, Queue)
          ##deleteat!(Queue,pos)
          ##delete!(Queue,current)

          ##current=Queue[end]
          stop = true
          ##Gomorycutsrounds = 0
          ##current=length(tree)

          #global NODECOUNT += 1
          #setlowerbound(x[f,g],0)
          #setupperbound(x[f,g],0)
          #push!(varxsolutionfixedzero,b) 
        elseif maximum(dsol-floor.(dsol)) > 0.00001
          a = maximum(dsol-floor.(dsol))
          f=0
          for i in 1:length(DESTC)
            if dsol[i]-floor.(dsol[i])== a
              f = i
              break
            end
          end
          #println("Sol d still frac, so split into two nodes using max frac", f)  
          push!(tree,TreeNode(current,[],[],[],[],[],[],[f],dsol,resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result
          push!(tree,TreeNode(current,[],[],[],[],[],[f],[],dsol,resultz,[],const11))
          push!(tree[current].children,length(tree))
          #push!(Queue,length(tree))
          Queue[length(tree)]=result

          #pos = findfirst(x -> x == current, Queue)
          #deleteat!(Queue,pos)
          #delete!(Queue,current)

          #current=Queue[end]

          stop = true
          #Gomorycutsrounds = 0
          #current=length(tree)

          #global NODECOUNT += 1
          #setlowerbound(d[f],0)
          #setupperbound(d[f],0)
          #push!(vardsolutionfixedzero,f) 
        else
          println("situation not covered")
          read(STDIN,Char)
        end
        end #stop2
      global TimeSce += toq()
          tic()
      end # !mincost
      stop && break
    end #end while true 'solve)
    return model, result 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate Prices, Routes and Scenarios to be used
  ###############################################
  #Create vectors for all possible scenarios
  global scenarioall = Vector{Vector{Int}}() 
  fillscenario(scenarioall,zeros(length(DESTC)),length(DESTC),1)
  println("All possible Scenario vectors created of size ", length(scenarioall) )
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenarioall)
    if sum(scenarioall[s]) <= CAPV
      groupe = []
      for c in 1:length(scenarioall[s])
        if scenarioall[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #################################################

  ###############################################
  #Re-Create vectors for nly initial scenarios - use col  cut generation afterwards
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  push!(scenario, zeros(1:length(DESTC)))   #Scenario 1 has to be all zeros. All customers present. Scenario 1 is used for calculations later. 
  for i = 1:length(DESTC)
    SAMPLETEST = zeros(1:length(DESTC))
    SAMPLETEST[i]=1
    push!(scenario, SAMPLETEST)
  end
  push!(scenario, ones(1:length(DESTC)))
  ################################################
  #scenario = scenarioall

  #Bypass REWOD and define new prices to pay ODs
  DESTC1=union(DESTC,[1])
  global PRICEOD
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance. 
  global PRICEOD = Inf*ones(1:length(DESTC))
  for i in DESTC
    for j in DESTC1
      for r in DESTC1
        if (j !=r && j != i && r != i) || (j==1 && r==1)
          if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
          end
        end
      end
    end
  end
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC1),length(DESTC1))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B6 = zeros(length(route),length(DESTC1),length(DESTC1)) 
  for r in 1:length(route)
    B6[r,route[r][end],length(DESTC1)] = 1
    B6[r,length(DESTC1),route[r][1]] = 1
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      B2[r,length(DESTC1),i] = 1
      B2[r,i,length(DESTC1)] = 1
      if findfirst(x -> x==i, route[r]) != 0
        if i != route[r][end]
          B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]] = 1
          #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
        end
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end 
  #####################################################################
  #read(STDIN,Char)
 
  #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  routesused = []
  for r in 1:length(route)
    if length(route[r]) <= 1
       push!(routesused,r)
    end
  end
  #routesused = collect(1:length(route))
  scenariosused = collect(1:length(scenario))
  #Initialize branch-and-bound tree and queue
  #Queue = Vector{Int}()
  Queue = PriorityQueue()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,scenariosused,[],[],[],[],[],[],[],[])) 
  #push!(Queue,1)
  Queue[1]= 0
  #Start processing queue in a Deep First mode
  global heursol
  global globalincumbentvalue = +Inf
  global globalincumbentsolution = []
    
  masterx = 0
  master = []
  #always process last element of queue
  while length(Queue)>0
    #current=Queue[end]
    current= dequeue!(Queue)
    (master, masterx) =process(current)
    global NODECOUNT += 1
    if length(Queue) > 0 && globalincumbentvalue != +Inf
      a,b=peek(Queue)
      if (globalincumbentvalue-b)/globalincumbentvalue < 0.02
        break 
      end
    end
    #if TimeMain +TimePric+TimeSce >= 10800
    #  break
    #end
  end
 
  for i in 1:length(globalincumbentsolution)
    if globalincumbentsolution[i] > 0.0001 && globalincumbentsolution[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  println("dsol= ", globalincumbentsolution)
  return globalincumbentvalue,globalincumbentsolution
end #end N9 function






#DDUEXACTN8, equal to N3, using deterministic rounded capacity cut ; but not using mixed best bound strategy with basis reconstruction  and not using reduced cost variable fixing
#DDUEXACTN3: Same as N2, but now introduce col gen and sce gen to b&b
#N2 was same as N1, but now create our own B&B on variables first xij and then di. Extended formulation. Now decide compensation at master problem, list all routes and scenaros at once, no column generation or cut ...solve relaxed MILP...added B&B)
#correlated marginals, worst case approach, decision dependent probabilistic VRP; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL
#This time vehicle sohould satisfy capacity for all scenarios. Only routes within capacity are generated a priori and symetry eliminated (as in N1)
function solveSTODDUEXACTN8(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vardsolutionfixedzero = []
    vardsolutionfixedone = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vardsolutionfixedzero = [vardsolutionfixedzero;tree[node].vardfixedsolutionzero]
      vardsolutionfixedone = [vardsolutionfixedone;tree[node].vardfixedsolutionone]
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
    #println("vardsolutionfixedzero =",vardsolutionfixedzero)
    #println("vardsolutionfixedone =",vardsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    #Now formulate problem (relaxed problem now)

    z = Vector{Variable}(length(routesused))
    y = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const4 = Vector{ConstraintRef}(length(cenariosused))
    const6 = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const7 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const8 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const11 = Vector{ConstraintRef}()

    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=7200))
    @variable(model,  s >= 0)
    #@variable(model,  z[r in routesused]>=0)
    #@variable(model,  y[r in routesused,w in cenariosused,i in 1:length(DESTC1),j in 1:length(DESTC1); i != j] >=0)

    for r in 1:length(routesused)
        z[r] = @variable(model, lowerbound = 0, upperbound =1)
    end
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
      if i != j
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1)
      else
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
      end
    end
    @variable(model,  x[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j]>=0)
    @variable(model,  u[1:length(DESTC)]<=0)
    @variable(model,  t[1:length(DESTC)]<=0)
    @variable(model,  d[1:length(DESTC)]>=0)

    for i in vardsolutionfixedzero
      setlowerbound(d[i],0)
      setupperbound(d[i],0) 
    end
    for i in vardsolutionfixedone
      setlowerbound(d[i],1)
      setupperbound(d[i],1)
    end

    for i in varxsolutionfixedzero
      setlowerbound(x[i[1],i[2]],0)
      setupperbound(x[i[1],i[2]],0) 
    end
    for i in varxsolutionfixedone
      setlowerbound(x[i[1],i[2]],1)
      setupperbound(x[i[1],i[2]],1)
    end  
 
    @objective(model, Min, s + sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) + sum(PROBOD2[i]*t[i] for i in 1:length(DESTC)) )

    @constraint(model, const1[i in 1:length(DESTC)], t[i]  >=  -10^4*d[i]) #be careful with cplex integrality tolerance
    @constraint(model, const2[i in 1:length(DESTC)], t[i] >= u[i] )

    for w in 1:length(cenariosused)
      const4[w] = @constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[cenariosused[w]][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][w][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
    end

    @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)
    @constraint(model, const9[i in 1:length(DESTC)], d[i] <= 1)
    @constraint(model,  const10[i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], x[i,j] - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
    
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
      if i != j
        const6[r][w][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[cenariosused[w]][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
      else
        const6[r][w][i][j] = @constraint(model, 0 == 0)
      end
    end

    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC)
      const7[r][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
      const8[r][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
    end

 

    while true
      status = solve(model)
      global TimeMain += toq()
      tic()
      println("Solved Node ",current," problem ", getobjectivevalue(model), "; Time is ",TimeMain, " and length Queue is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #read(STDIN,Char)
      #print(";",NODECOUNT,":",TimeMain,"...")
      result = getobjectivevalue(model)
      resultz = getvalue(z)
      #read(STDIN,Char)
      #println("d= ",getvalue(d))
      #read(STDIN,Char)
      #for r in 1:length(routesused)
      #  for w in 1:length(cenariosused)
      #    for i in 1:length(DESTC1)
      #      for j in 1:length(DESTC1)
      #        if i != j
      #          if getvalue(y[r][w][i][j])>0.0001
      #            #print("y [",r,",",w,",",i,",",j," ]:",getvalue(y[r,w,i,j]))
      #            #read(STDIN,Char)
      #          end
      #        end
      #      end
      #    end
      #  end 
      #end
      #println()
      #for r in 1:length(routesused)
      #  #println(r," ;",route[r])
      #  if getvalue(z[r]) >= 0.00001
      #    println(route[r],": z[ ",r," ]= ",getvalue(z[r])," ;") 
      #    read(STDIN,Char)
      #  end
      #end  #r 
      
      stop = false
      mincost = false
      if status != :Optimal
        #prune   
        pos = findfirst(x -> x == current, Queue)
        deleteat!(Queue,pos)
        println("Infeasible...prune ", status)
        read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolution
     

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
              xsol[i,j]= xsol1[i,j]
              #print("xsol[",i," ,",j," ]",xsol[i,j])
            end
          end
        end
      end
      #for i=1:length(DESTC1)
      #  println(xsol[i,:])
      #end
      #read(STDIN,Char)
      dsol = getvalue(d)
      for i=1:length(DESTC)
        if dsol[i]<=0.00001 
          dsol[i] = 0
        elseif  dsol[i] >=0.99999
          dsol[i]= 1
        end
      end

      
      E1dual = getdual(const5)
      G1dual = getdual(const10)
      maxcol = 0
      #for r in 1:length(routesused)    
      #  if sum( B1[routesused[r],i]*E1dual[i] for i in 1:length(DESTC)) - sum(B2[routesused[r],i,j]*G1dual[i,j] for i in 1:length(DESTC),j in 1:length(DESTC) if  i != j)< -0.1
      #     setlowerbound(z[r],0)
      #     setupperbound(z[r],0) 
      #  end 
      #end

      #Section for Reduced Cost variable Fixing. Var need to be zero and use primal heuristic of this node as upper bound
      #Need still to verify why can be fixed for the whole subtree
      #ddual = getdual(d)
      #for i in 1:length(ddual)
      #  if result + ddual[i] > globalincumbentvalue && !in(i,vardsolutionfixedone) && !in(i,vardsolutionfixedzero) #&& dsol[i] <= 0.9
      #    #println("node ",current," d solution for customer ",i, " will not improve")
      #    #read(STDIN,Char)
      #   if dsol[i] >= 0.1 && dsol[i] <= 0.9
      #      println("valeu a pena")
      #      read(STDIN,Char)
      #    end
      #    setlowerbound(d[i],0)
      #    setupperbound(d[i],0)  
      #  end
      #end
      for r in 1:length(route)
        if !in(r,routesused)
          if sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)) + sum(B6[r,i,j]*G1dual[i,j] for i in 1:length(DESTC1),j in 1:length(DESTC1) if  i != j)>= 0.00001
            #println("Node =",current," ,vou inserir route ", r, " valor constraint= ", sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)))
            push!(routesused,r)
            push!(tree[current].addedsolutions,r)
            global ROUTCOUNT += 1
            mincost = true
            #read(STDIN,Char)
            tes = [B1[r,i] for i in 1:length(DESTC)]
            tes1 = [const5[i] for i in 1:length(DESTC) ]
            for i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i !=j
                push!(tes, -B6[r,i,j])
                push!(tes1, const10[i,j])
              end
            end
            #println(tes)
            #println(tes1)
            push!(z,@variable(model,  lowerbound = 0, upperbound =1, objective = 0.0, inconstraints = tes1, coefficients = tes ))
            #setvalue(y[findfirst(x -> x == best_s, solutionsused)],1)

            ynew = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:1]
            for  w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                ynew[1][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1, objective = 0.0, inconstraints = [const4[w] ], coefficients = [-REWV*A[DESTC1[i],DESTC1[j]]])
              else
                ynew[1][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
              end
            end
            y = [y;ynew]


            const6new = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:1]
            const7new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:1]
            const8new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:1]
            for  w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
              if i != j
                const6new[1][w][i][j] = @constraint(model, -(B2[r,i,j] + sum(scenario[cenariosused[w]][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*z[end] + y[end][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
              else
                const6new[1][w][i][j] = @constraint(model, 0 == 0)
              end
            end
            for w in 1:length(cenariosused),i in 1:length(DESTC)
              const7new[1][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[r,k,i]+0+B2[r,i,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[end] + y[end][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
              const8new[1][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[end] + y[end][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
            end
            
            const6 = [const6;const6new]
            const7 = [const7;const7new]
            const8 = [const8;const8new]
            global TimePric += toq()
            tic()
            if maxcol == 5
             break
            end
            maxcol +=1
          end
        end #r in route
      end
  
      if !mincost

        if  result > globalincumbentvalue + 0.0001
         #prune  since a lower bound is already worse then incumbent 
         pos = findfirst(x -> x == current, Queue)
         deleteat!(Queue,pos)
         #println("Prune by Bound =", result, ", ", globalincumbentvalue)
         #read(STDIN,Char) 
         stop = true
         break
        end

        stop2=false
        if maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(dsol-floor.(dsol)) <= 0.00001
          #new integer solution
          #println("Integer Solution Found")

          ##############################################
          #Formulate Separation Problem and Solve 
          #Get First Stage results
          ssol = getvalue(s)
          usol= getvalue(u)
          for i in 1:length(resultz) #check for error and display
            if resultz[i] > 0.0001 && resultz[i] < 0.99
              println("SOLUCAO NAO INTEIRA!! ",i)
              read(STDIN,Char) 
            end
          end
      
          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1) 
 
          @objective(model2, Min, ssol + sum(scen[i]*usol[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scen[i]*dsol[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
          @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

          @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

          @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


          @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i]) )
          @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

          @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
          @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

          @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i]) )

          @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

          @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0.001], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )




 
          #println("Will solve Separation Problem for integer solution")
           
          solve(model2)

 
          #println("solve model 2",getobjectivevalue(model2))
          if getobjectivevalue(model2) <= -0.0001 #-0.05

            ##free all reduced cost fixed variables
            #for i in 1:length(DESTC)
            #  if !in(i,vardsolutionfixedzero) && !in(i,vardsolutionfixedone)
            #    setlowerbound(d[i],0)
            #    setupperbound(d[i],1)  
            #  end
            #end

            #println(getvalue(scen), " get= ", getobjectivevalue(model2))
            #read(STDIN,Char)
            #Add new scenario and continue
            scenarionew = round.(getvalue(scen))
            #println("Integer Solution not valid. New scenario: ", scenarionew)
            push!(scenario,scenarionew)
            #Update Current node information on scenarios used
            push!(tree[current].addedcenarios,size(scenario,1))
            push!(cenariosused,size(scenario,1))
            #Create new variables
            ynew = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:1] for _ in 1:length(routesused)]
            for  r in 1:length(routesused),i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                ynew[r][1][i][j] = @variable(model,  lowerbound = 0, upperbound =1)
              else
                ynew[r][1][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
              end
            end
            for r in 1:length(routesused)
             y[r] = [y[r];ynew[r]]
            end 
            #Create new constraints
            push!(const4, @constraint(model,  s + sum(scenario[end][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[end][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][end][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))


            
            const6new = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:1] for _ in 1:length(routesused)]
            const7new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:1] for _ in 1:length(routesused)]
            const8new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:1] for _ in 1:length(routesused)]
            for  r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC)
              if i != j
                const6new[r][1][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[end][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][end][i][j] >= 1- scenario[end][i] + 1 - scenario[end][j] - length(DESTC))
              else
                const6new[r][1][i][j] = @constraint(model, 0 == 0)
              end
            end
            for r in 1:length(routesused),i in 1:length(DESTC)
              const7new[r][1][i] = @constraint(model,  -(0 + sum(scenario[end][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][end][length(DESTC1)][i] >=  1 - scenario[end][i] - length(DESTC) + 1)
              const8new[r][1][i]= @constraint(model, -(0 + sum(scenario[end][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][end][i][length(DESTC1)] >= 1- scenario[end][i]  - length(DESTC) + 1)
            end
            for r in 1:length(routesused)
              const6[r]=[const6[r];const6new[r]]
              const7[r]=[const7[r];const7new[r]]
              const8[r]=[const8[r];const8new[r]]
            end
            stop2=true
          else # end getobj < +.05
          #end verification of scenario insertion for probable integer solution
          ##############################################
          #Verify if it can be new incubent
            if result < globalincumbentvalue
              #println("New Incumbent Found")
              globalincumbentvalue = result
              globalincumbentsolution = dsol #resultz
            end 
            pos = findfirst(x -> x == current, Queue)
            deleteat!(Queue,pos)
            stop=true
            break
          end
          global TimeSce += toq()
          tic()
        else
          #insert rounded capacity cuts if possible (deterministic)
          for sce in 1:length(scenarioall) #eliminate scenario all 0 and all 1
            if scenarioall[sce] != zeros(1:length(scenarioall[sce])) && scenarioall[sce] != ones(1:length(scenarioall[sce]))
              sumtemp=0
              #println("test scenario= ", scenarioall[s])
              #read(STDIN,Char)
              for i=1:length(scenarioall[sce])
                if scenarioall[sce][i]== 1
                  sumtemp += xsol[i,length(DESTC1)]
                  #sumtemp += xsol[length(DESTC1),i]
                end
                for j=1:length(scenarioall[sce])
                  if scenarioall[sce][i]== 1 && scenarioall[sce][j]== 0
                     sumtemp += xsol[i,j]
                  end
                end
              end
              #println("resumo ",sumtemp + 0.00001,",", ceil(sum(scenarioall[s])/CAPV))
              #read(STDIN,Char)
              if sumtemp + 0.00001 < ceil(sum(scenarioall[sce])/CAPV)
                #insert
                push!(const11, @constraint(model, sum(sum(x[i,j] for i=1:length(DESTC) if scenarioall[sce][i]== 1 && (j==length(DESTC1) || scenarioall[sce][j]== 0) ) for j=1:length(DESTC1)  ) >= ceil(sum(scenarioall[sce])/CAPV) ) )
                println("insert RCC and out on node ", current)
                #println(const11[end])
                #read(STDIN,Char)
                stop2=true
                break
              end
            end
         end
        end
        
        
        if !stop2  #time for branching
        if maximum(xsol-floor.(xsol)) > 0.00001
          if length(vardsolutionfixedzero) != 0 || length(vardsolutionfixedone) != 0
            println("d ja estava fixed e veio xij frac")
            println(vardsolutionfixedzero)
            println(vardsolutionfixedone)
            read(STDIN,Char)
          end 
          a = maximum(xsol-floor.(xsol))
          f=0
          g=0
          for i in 1:length(DESTC1),j in 1:length(DESTC1)
            if i != j
              if xsol[i,j]-floor.(xsol[i,j])== a
                f = i
                g = j
                break
              end
            end
          end
          #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
          b =[f;g] 
          push!(tree,TreeNode(current,[],[],[],[],[b],[],[]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))
          push!(tree,TreeNode(current,[],[],[],[b],[],[],[]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))

          pos = findfirst(x -> x == current, Queue)
          deleteat!(Queue,pos)

          current=Queue[end]
          global NODECOUNT += 1
          setlowerbound(x[f,g],0)
          setupperbound(x[f,g],0)
          push!(varxsolutionfixedzero,b) 
        elseif maximum(dsol-floor.(dsol)) > 0.00001
          a = maximum(dsol-floor.(dsol))
          f=0
          for i in 1:length(DESTC)
            if dsol[i]-floor.(dsol[i])== a
              f = i
              break
            end
          end
          #println("Sol d still frac, so split into two nodes using max frac", f)  
          push!(tree,TreeNode(current,[],[],[],[],[],[],[f]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))
          push!(tree,TreeNode(current,[],[],[],[],[],[f],[]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))

          pos = findfirst(x -> x == current, Queue)
          deleteat!(Queue,pos)

          current=Queue[end]
          global NODECOUNT += 1
          setlowerbound(d[f],0)
          setupperbound(d[f],0)
          push!(vardsolutionfixedzero,f) 
        else
          println("situation not covered")
          read(STDIN,Char)
        end
        end #stop2
      end # !mincost
      stop && break
    end #end while true 'solve)
    return model, resultz 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate Prices, Routes and Scenarios to be used
  ###############################################
  #Create vectors for all possible scenarios
  global scenarioall = Vector{Vector{Int}}() 
  fillscenario(scenarioall,zeros(length(DESTC)),length(DESTC),1)
  println("All possible Scenario vectors created of size ", length(scenarioall) )
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenarioall)
    if sum(scenarioall[s]) <= CAPV
      groupe = []
      for c in 1:length(scenarioall[s])
        if scenarioall[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #################################################

  ###############################################
  #Re-Create vectors for nly initial scenarios - use col  cut generation afterwards
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  for i = 1:length(DESTC)
    SAMPLETEST = zeros(1:length(DESTC))
    SAMPLETEST[i]=1
    push!(scenario, SAMPLETEST)
  end
  push!(scenario, ones(1:length(DESTC)))
  push!(scenario, zeros(1:length(DESTC)))
  ################################################

  #Bypass REWOD and define new prices to pay ODs
  DESTC1=union(DESTC,[1])
  global PRICEOD
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance. 
  global PRICEOD = Inf*ones(1:length(DESTC))
  for i in DESTC
    for j in DESTC1
      for r in DESTC1
        if (j !=r && j != i && r != i) || (j==1 && r==1)
          if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
          end
        end
      end
    end
  end
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC),length(DESTC))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B6 = zeros(length(route),length(DESTC1),length(DESTC1)) 
  for r in 1:length(route)
    B6[r,route[r][end],length(DESTC1)] = 1
    B6[r,length(DESTC1),route[r][1]] = 1
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        if i != route[r][end]
          B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]] = 1
          #print("B6[ ",r," ,",i," ,",route[r][findfirst(x -> x==i, route[r])+1]," ]= ",B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]], "; ",route[r])
        end
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end 
  #####################################################################
  #read(STDIN,Char)
 
  #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  routesused = []
  for r in 1:length(route)
    if length(route[r]) <= 1
       push!(routesused,r)
    end
  end
  #routesused = collect(1:length(route))
  scenariosused = collect(1:length(scenario))
  #Initialize branch-and-bound tree and queue
  Queue = Vector{Int}()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,scenariosused,[],[],[],[])) 
  push!(Queue,1)
  #Start processing queue in a Deep First mode

  global globalincumbentvalue = +Inf
  global globalincumbentsolution = []
    
  masterx = []
  master = []
  #always process last element of queue
  while length(Queue)>0
    current=Queue[end]
    (master, masterx) =process(current)
    global NODECOUNT += 1
    #if TimeMain +TimePric+TimeSce >= 10800
    #  break
    #end
  end
 
  for i in 1:length(globalincumbentsolution)
    if globalincumbentsolution[i] > 0.0001 && globalincumbentsolution[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  println("dsol= ", globalincumbentsolution)
  return globalincumbentvalue,globalincumbentsolution
end #end N8 function


#DDUEXACTN7: same as N3, but no compensation decision, only first option for price
#DDUEXACTN3: Same as N2, but now introduce col gen and sce gen to b&b
#N2 was same as N1, but now create our own B&B on variables first xij and then di. Extended formulation. Now decide compensation at master problem, list all routes and scenaros at once, no column generation or cut ...solve relaxed MILP...added B&B)
#correlated marginals, worst case approach, decision dependent probabilistic VRP; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL
#This time vehicle sohould satisfy capacity for all scenarios. Only routes within capacity are generated a priori and symetry eliminated (as in N1)
function solveSTODDUEXACTN7(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vardsolutionfixedzero = []
    vardsolutionfixedone = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vardsolutionfixedzero = [vardsolutionfixedzero;tree[node].vardfixedsolutionzero]
      vardsolutionfixedone = [vardsolutionfixedone;tree[node].vardfixedsolutionone]
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
    #println("vardsolutionfixedzero =",vardsolutionfixedzero)
    #println("vardsolutionfixedone =",vardsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    #Now formulate problem (relaxed problem now)

    z = Vector{Variable}(length(routesused))
    y = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const4 = Vector{ConstraintRef}(length(cenariosused))
    const6 = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const7 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const8 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]


    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=7200))
    @variable(model,  s >= 0)
    #@variable(model,  z[r in routesused]>=0)
    #@variable(model,  y[r in routesused,w in cenariosused,i in 1:length(DESTC1),j in 1:length(DESTC1); i != j] >=0)

    for r in 1:length(routesused)
        z[r] = @variable(model, lowerbound = 0, upperbound =1)
    end
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
      if i != j
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1)
      else
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
      end
    end
    @variable(model,  x[i in 1:length(DESTC),j in 1:length(DESTC); i != j]>=0)
    @variable(model,  u[1:length(DESTC)]<=0)


    for i in varxsolutionfixedzero
      setlowerbound(x[i[1],i[2]],0)
      setupperbound(x[i[1],i[2]],0) 
    end
    for i in varxsolutionfixedone
      setlowerbound(x[i[1],i[2]],1)
      setupperbound(x[i[1],i[2]],1)
    end  
 
    @objective(model, Min, s + sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) )


    for w in 1:length(cenariosused)
      const4[w] = @constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC))  - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][w][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
    end

    @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)
    @constraint(model,  const10[i in 1:length(DESTC),j in 1:length(DESTC); i != j], x[i,j] - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
    
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
      if i != j
        const6[r][w][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[cenariosused[w]][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
      else
        const6[r][w][i][j] = @constraint(model, 0 == 0)
      end
    end

    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC)
      const7[r][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
      const8[r][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
    end

 

    while true
      status = solve(model)
      global TimeMain += toq()
      tic()
      println("Solved Node ",current," problem ", getobjectivevalue(model), "; Time is ",TimeMain, " and length Queue is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #read(STDIN,Char)
      #print(";",NODECOUNT,":",TimeMain,"...")
      result = getobjectivevalue(model)
      resultz = getvalue(z)
      #read(STDIN,Char)
      #println("d= ",getvalue(d))
      #read(STDIN,Char)
      #for r in 1:length(routesused)
      #  for w in 1:length(cenariosused)
      #    for i in 1:length(DESTC1)
      #      for j in 1:length(DESTC1)
      #        if i != j
      #          if getvalue(y[r][w][i][j])>0.0001
      #            #print("y [",r,",",w,",",i,",",j," ]:",getvalue(y[r,w,i,j]))
      #            #read(STDIN,Char)
      #          end
      #        end
      #      end
      #    end
      #  end 
      #end
      #println()
      #for r in 1:length(routesused)
      #  #println(r," ;",route[r])
      #  if getvalue(z[r]) >= 0.00001
      #    println(route[r],": z[ ",r," ]= ",getvalue(z[r])," ;") 
      #    read(STDIN,Char)
      #  end
      #end  #r 
      
      stop = false
      mincost = false
      if status != :Optimal
        #prune   
        pos = findfirst(x -> x == current, Queue)
        deleteat!(Queue,pos)
        println("Infeasible...prune ", status)
        read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolution
     

      xsol1 = getvalue(x)
      xsol = zeros(length(DESTC),length(DESTC))
      for i=1:length(DESTC)
        for j=1:length(DESTC)
          if i != j
            if xsol1[i,j]<=0.00001 
              xsol[i,j] = 0
            elseif  xsol1[i,j] >=0.99999
              xsol[i,j]= 1
            else
              xsol[i,j]= xsol1[i,j]
              #print("xsol[",i," ,",j," ]",xsol[i,j])
            end
          end
        end
      end
      #println()
      #read(STDIN,Char)
      

      
      E1dual = getdual(const5)
      G1dual = getdual(const10)
      maxcol = 0
      #for r in 1:length(routesused)    
      #  if sum( B1[routesused[r],i]*E1dual[i] for i in 1:length(DESTC)) - sum(B2[routesused[r],i,j]*G1dual[i,j] for i in 1:length(DESTC),j in 1:length(DESTC) if  i != j)< -0.1
      #     setlowerbound(z[r],0)
      #     setupperbound(z[r],0) 
      #  end 
      #end
      for r in 1:length(route)
        if !in(r,routesused)
          if sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)) + sum(B6[r,i,j]*G1dual[i,j] for i in 1:length(DESTC),j in 1:length(DESTC) if  i != j)>= 0.00001
            #println("Node =",current," ,vou inserir route ", r, " valor constraint= ", sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)))
            push!(routesused,r)
            push!(tree[current].addedsolutions,r)
            global ROUTCOUNT += 1
            mincost = true
            #read(STDIN,Char)
            tes = [B1[r,i] for i in 1:length(DESTC)]
            tes1 = [const5[i] for i in 1:length(DESTC) ]
            for i in 1:length(DESTC),j in 1:length(DESTC)
              if i !=j
                push!(tes, -B6[r,i,j])
                push!(tes1, const10[i,j])
              end
            end
            #println(tes)
            #println(tes1)
            push!(z,@variable(model,  lowerbound = 0, upperbound =1, objective = 0.0, inconstraints = tes1, coefficients = tes ))
            #setvalue(y[findfirst(x -> x == best_s, solutionsused)],1)

            ynew = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:1]
            for  w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                ynew[1][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1, objective = 0.0, inconstraints = [const4[w] ], coefficients = [-REWV*A[DESTC1[i],DESTC1[j]]])
              else
                ynew[1][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
              end
            end
            y = [y;ynew]


            const6new = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:1]
            const7new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:1]
            const8new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:1]
            for  w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
              if i != j
                const6new[1][w][i][j] = @constraint(model, -(B2[r,i,j] + sum(scenario[cenariosused[w]][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*z[end] + y[end][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
              else
                const6new[1][w][i][j] = @constraint(model, 0 == 0)
              end
            end
            for w in 1:length(cenariosused),i in 1:length(DESTC)
              const7new[1][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[r,k,i]+0+B2[r,i,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[end] + y[end][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
              const8new[1][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[end] + y[end][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
            end
            
            const6 = [const6;const6new]
            const7 = [const7;const7new]
            const8 = [const8;const8new]
            global TimePric += toq()
            tic()
            if maxcol == 5
             break
            end
            maxcol +=1
          end
        end #r in route
      end
  
      if !mincost
        if  result > globalincumbentvalue + 0.0001
         #prune  since a lower bound is already worse then incumbent 
         pos = findfirst(x -> x == current, Queue)
         deleteat!(Queue,pos)
         #println("Prune by Bound =", result, ", ", globalincumbentvalue)
         #read(STDIN,Char) 
         stop = true
         break
        end 
        if maximum(xsol-floor.(xsol)) <= 0.00001
          #new integer solution
          #println("Integer Solution Found")

          ##############################################
          #Formulate Separation Problem and Solve 
          #Get First Stage results
          ssol = getvalue(s)
          usol= getvalue(u)
          for i in 1:length(resultz) #check for error and display
            if resultz[i] > 0.0001 && resultz[i] < 0.99
              println("SOLUCAO NAO INTEIRA!! ",i)
              read(STDIN,Char) 
            end
          end
      
          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1) 
 
          @objective(model2, Min, ssol + sum(scen[i]*usol[i] for i in 1:length(DESTC))  - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
          @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

          @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

          @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


          @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i]) )
          @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

          @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
          @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

          @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i]) )

          @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

          @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0.001], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )




 
          #println("Will solve Separation Problem for integer solution")
           
          solve(model2)

 
          #println("solve model 2",getobjectivevalue(model2))
          if getobjectivevalue(model2) <= -0.0001 #-0.05
            #println(getvalue(scen), " get= ", getobjectivevalue(model2))
            #read(STDIN,Char)
            #Add new scenario and continue
            scenarionew = round.(getvalue(scen))
            #println("Integer Solution not valid. New scenario: ", scenarionew)
            push!(scenario,scenarionew)
            #Update Current node information on scenarios used
            push!(tree[current].addedcenarios,size(scenario,1))
            push!(cenariosused,size(scenario,1))
            #Create new variables
            ynew = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:1] for _ in 1:length(routesused)]
            for  r in 1:length(routesused),i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                ynew[r][1][i][j] = @variable(model,  lowerbound = 0, upperbound =1)
              else
                ynew[r][1][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
              end
            end
            for r in 1:length(routesused)
             y[r] = [y[r];ynew[r]]
            end 
            #Create new constraints
            push!(const4, @constraint(model,  s + sum(scenario[end][i]*u[i] for i in 1:length(DESTC))- sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][end][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))


            
            const6new = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:1] for _ in 1:length(routesused)]
            const7new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:1] for _ in 1:length(routesused)]
            const8new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:1] for _ in 1:length(routesused)]
            for  r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC)
              if i != j
                const6new[r][1][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[end][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][end][i][j] >= 1- scenario[end][i] + 1 - scenario[end][j] - length(DESTC))
              else
                const6new[r][1][i][j] = @constraint(model, 0 == 0)
              end
            end
            for r in 1:length(routesused),i in 1:length(DESTC)
              const7new[r][1][i] = @constraint(model,  -(0 + sum(scenario[end][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][end][length(DESTC1)][i] >=  1 - scenario[end][i] - length(DESTC) + 1)
              const8new[r][1][i]= @constraint(model, -(0 + sum(scenario[end][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][end][i][length(DESTC1)] >= 1- scenario[end][i]  - length(DESTC) + 1)
            end
            for r in 1:length(routesused)
              const6[r]=[const6[r];const6new[r]]
              const7[r]=[const7[r];const7new[r]]
              const8[r]=[const8[r];const8new[r]]
            end
          else # end getobj < +.05
          #end verification of scenario insertion for probable integer solution
          ##############################################
          #Verify if it can be new incubent
            if result < globalincumbentvalue
              #println("New Incumbent Found")
              globalincumbentvalue = result
              globalincumbentsolution = resultz
            end 
            pos = findfirst(x -> x == current, Queue)
            deleteat!(Queue,pos)
            stop=true
            break
          end
          global TimeSce += toq()
          tic()
        elseif maximum(xsol-floor.(xsol)) > 0.00001
          if length(vardsolutionfixedzero) != 0 || length(vardsolutionfixedone) != 0
            println("d ja estava fixed e veio xij frac")
            println(vardsolutionfixedzero)
            println(vardsolutionfixedone)
            read(STDIN,Char)
          end 
          a = maximum(xsol-floor.(xsol))
          f=0
          g=0
          for i in 1:length(DESTC),j in 1:length(DESTC)
            if i != j
              if xsol[i,j]-floor.(xsol[i,j])== a
                f = i
                g = j
                break
              end
            end
          end
          #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
          b =[f;g] 
          push!(tree,TreeNode(current,[],[],[],[],[b],[],[]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))
          push!(tree,TreeNode(current,[],[],[],[b],[],[],[]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))

          pos = findfirst(x -> x == current, Queue)
          deleteat!(Queue,pos)

          current=Queue[end]
          global NODECOUNT += 1
          setlowerbound(x[f,g],0)
          setupperbound(x[f,g],0) 
        else
          println("situation not covered")
          read(STDIN,Char)
        end
      end # !mincost
      stop && break
    end #end while true 'solve)
    return model, resultz 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate Prices, Routes and Scenarios to be used
  ###############################################
  #Create vectors for all possible scenarios
  global scenario = Vector{Vector{Int}}() 
  fillscenario(scenario,zeros(length(DESTC)),length(DESTC),1)
  println("All possible Scenario vectors created of size ", length(scenario) )
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenario)
    if sum(scenario[s]) <= CAPV
      groupe = []
      for c in 1:length(scenario[s])
        if scenario[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #################################################

  ###############################################
  #Re-Create vectors for nly initial scenarios - use col  cut generation afterwards
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  for i = 1:length(DESTC)
    SAMPLETEST = zeros(1:length(DESTC))
    SAMPLETEST[i]=1
    push!(scenario, SAMPLETEST)
  end
  push!(scenario, ones(1:length(DESTC)))
  push!(scenario, zeros(1:length(DESTC)))
  ################################################

  #Bypass REWOD and define new prices to pay ODs
  DESTC1=union(DESTC,[1])
  global PRICEOD
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance. 
  global PRICEOD = Inf*ones(1:length(DESTC))
  for i in DESTC
    for j in DESTC1
      for r in DESTC1
        if (j !=r && j != i && r != i) || (j==1 && r==1)
          if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
          end
        end
      end
    end
  end
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC),length(DESTC))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC)) 
  global B6 = zeros(length(route),length(DESTC),length(DESTC))
  for r in 1:length(route)
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        if i != route[r][end]
          B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]] = 1
        end
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end 
  #####################################################################
 
  #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  routesused = []
  for r in 1:length(route)
    if length(route[r]) <= 1
       push!(routesused,r)
    end
  end
  #routesused = collect(1:length(route))
  scenariosused = collect(1:length(scenario))
  #Initialize branch-and-bound tree and queue
  Queue = Vector{Int}()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,scenariosused,[],[],[],[])) 
  push!(Queue,1)
  #Start processing queue in a Deep First mode

  global globalincumbentvalue = +Inf
  global globalincumbentsolution = []
    
  x = []
  master = []
  #always process last element of queue
  while length(Queue)>0
    current=Queue[end]
    (master, x) =process(current)
    global NODECOUNT += 1
    #if TimeMain +TimePric+TimeSce >= 10800
    #  break
    #end
  end
 
  for i in 1:length(globalincumbentsolution)
    if globalincumbentsolution[i] > 0.0001 && globalincumbentsolution[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  return globalincumbentvalue,globalincumbentsolution
end #end N7 function

#DDUEXACTN3: Same as N2, but now introduce col gen and sce gen to b&b
#N2 was same as N1, but now create our own B&B on variables first xij and then di. Extended formulation. Now decide compensation at master problem, list all routes and scenaros at once, no column generation or cut ...solve relaxed MILP...added B&B)
#correlated marginals, worst case approach, decision dependent probabilistic VRP; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL
#This time vehicle sohould satisfy capacity for all scenarios. Only routes within capacity are generated a priori and symetry eliminated (as in N1)
function solveSTODDUEXACTN3(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vardsolutionfixedzero = []
    vardsolutionfixedone = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vardsolutionfixedzero = [vardsolutionfixedzero;tree[node].vardfixedsolutionzero]
      vardsolutionfixedone = [vardsolutionfixedone;tree[node].vardfixedsolutionone]
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
    #println("vardsolutionfixedzero =",vardsolutionfixedzero)
    #println("vardsolutionfixedone =",vardsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    #Now formulate problem (relaxed problem now)

    z = Vector{Variable}(length(routesused))
    y = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const4 = Vector{ConstraintRef}(length(cenariosused))
    const6 = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const7 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const8 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]


    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=7200))
    @variable(model,  s >= 0)
    #@variable(model,  z[r in routesused]>=0)
    #@variable(model,  y[r in routesused,w in cenariosused,i in 1:length(DESTC1),j in 1:length(DESTC1); i != j] >=0)

    for r in 1:length(routesused)
        z[r] = @variable(model, lowerbound = 0, upperbound =1)
    end
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
      if i != j
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1)
      else
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
      end
    end
    @variable(model,  x[i in 1:length(DESTC),j in 1:length(DESTC); i != j]>=0)
    @variable(model,  u[1:length(DESTC)]<=0)
    @variable(model,  t[1:length(DESTC)]<=0)
    @variable(model,  d[1:length(DESTC)]>=0)

    for i in vardsolutionfixedzero
      setlowerbound(d[i],0)
      setupperbound(d[i],0) 
    end
    for i in vardsolutionfixedone
      setlowerbound(d[i],1)
      setupperbound(d[i],1)
    end

    for i in varxsolutionfixedzero
      setlowerbound(x[i[1],i[2]],0)
      setupperbound(x[i[1],i[2]],0) 
    end
    for i in varxsolutionfixedone
      setlowerbound(x[i[1],i[2]],1)
      setupperbound(x[i[1],i[2]],1)
    end  
 
    @objective(model, Min, s + sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) + sum(PROBOD2[i]*t[i] for i in 1:length(DESTC)) )

    @constraint(model, const1[i in 1:length(DESTC)], t[i]  >=  -10^4*d[i]) #be careful with cplex integrality tolerance
    @constraint(model, const2[i in 1:length(DESTC)], t[i] >= u[i] )

    for w in 1:length(cenariosused)
      const4[w] = @constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[cenariosused[w]][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][w][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
    end

    @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)
    @constraint(model, const9[i in 1:length(DESTC)], d[i] <= 1)
    @constraint(model,  const10[i in 1:length(DESTC),j in 1:length(DESTC); i != j], x[i,j] - sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused))==0)
    
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
      if i != j
        const6[r][w][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[cenariosused[w]][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
      else
        const6[r][w][i][j] = @constraint(model, 0 == 0)
      end
    end

    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC)
      const7[r][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
      const8[r][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
    end

 

    while true
      status = solve(model)
      global TimeMain += toq()
      tic()
      println("Solved Node ",current," problem ", getobjectivevalue(model), "; Time is ",TimeMain, " and length Queue is ", length(Queue),". Incumbent is ", globalincumbentvalue,". Number of Routes/Scenario: ",length(routesused)," / ",length(cenariosused))
      #read(STDIN,Char)
      #print(";",NODECOUNT,":",TimeMain,"...")
      result = getobjectivevalue(model)
      resultz = getvalue(z)
      #read(STDIN,Char)
      #println("d= ",getvalue(d))
      #read(STDIN,Char)
      #for r in 1:length(routesused)
      #  for w in 1:length(cenariosused)
      #    for i in 1:length(DESTC1)
      #      for j in 1:length(DESTC1)
      #        if i != j
      #          if getvalue(y[r][w][i][j])>0.0001
      #            #print("y [",r,",",w,",",i,",",j," ]:",getvalue(y[r,w,i,j]))
      #            #read(STDIN,Char)
      #          end
      #        end
      #      end
      #    end
      #  end 
      #end
      #println()
      #for r in 1:length(routesused)
      #  #println(r," ;",route[r])
      #  if getvalue(z[r]) >= 0.00001
      #    println(route[r],": z[ ",r," ]= ",getvalue(z[r])," ;") 
      #    read(STDIN,Char)
      #  end
      #end  #r 
      
      stop = false
      mincost = false
      if status != :Optimal
        #prune   
        pos = findfirst(x -> x == current, Queue)
        deleteat!(Queue,pos)
        println("Infeasible...prune ", status)
        read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolution
     

      xsol1 = getvalue(x)
      xsol = zeros(length(DESTC),length(DESTC))
      for i=1:length(DESTC)
        for j=1:length(DESTC)
          if i != j
            if xsol1[i,j]<=0.00001 
              xsol[i,j] = 0
            elseif  xsol1[i,j] >=0.99999
              xsol[i,j]= 1
            else
              xsol[i,j]= xsol1[i,j]
              #print("xsol[",i," ,",j," ]",xsol[i,j])
            end
          end
        end
      end
      #println()
      #read(STDIN,Char)
      dsol = getvalue(d)
      for i=1:length(DESTC)
        if dsol[i]<=0.00001 
          dsol[i] = 0
        elseif  dsol[i] >=0.99999
          dsol[i]= 1
        end
      end

      
      E1dual = getdual(const5)
      G1dual = getdual(const10)
      maxcol = 0
      #for r in 1:length(routesused)    
      #  if sum( B1[routesused[r],i]*E1dual[i] for i in 1:length(DESTC)) - sum(B2[routesused[r],i,j]*G1dual[i,j] for i in 1:length(DESTC),j in 1:length(DESTC) if  i != j)< -0.1
      #     setlowerbound(z[r],0)
      #     setupperbound(z[r],0) 
      #  end 
      #end
      for r in 1:length(route)
        if !in(r,routesused)
          if sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)) - sum(B2[r,i,j]*G1dual[i,j] for i in 1:length(DESTC),j in 1:length(DESTC) if  i != j)>= 0.00001
            #println("Node =",current," ,vou inserir route ", r, " valor constraint= ", sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)))
            push!(routesused,r)
            push!(tree[current].addedsolutions,r)
            global ROUTCOUNT += 1
            mincost = true
            #read(STDIN,Char)
            tes = [B1[r,i] for i in 1:length(DESTC)]
            tes1 = [const5[i] for i in 1:length(DESTC) ]
            for i in 1:length(DESTC),j in 1:length(DESTC)
              if i !=j
                push!(tes, -B2[r,i,j])
                push!(tes1, const10[i,j])
              end
            end
            #println(tes)
            #println(tes1)
            push!(z,@variable(model,  lowerbound = 0, upperbound =1, objective = 0.0, inconstraints = tes1, coefficients = tes ))
            #setvalue(y[findfirst(x -> x == best_s, solutionsused)],1)

            ynew = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:1]
            for  w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                ynew[1][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1, objective = 0.0, inconstraints = [const4[w] ], coefficients = [-REWV*A[DESTC1[i],DESTC1[j]]])
              else
                ynew[1][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
              end
            end
            y = [y;ynew]


            const6new = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:1]
            const7new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:1]
            const8new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:1]
            for  w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
              if i != j
                const6new[1][w][i][j] = @constraint(model, -(B2[r,i,j] + sum(scenario[cenariosused[w]][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*z[end] + y[end][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
              else
                const6new[1][w][i][j] = @constraint(model, 0 == 0)
              end
            end
            for w in 1:length(cenariosused),i in 1:length(DESTC)
              const7new[1][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[r,k,i]+0+B2[r,i,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[end] + y[end][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
              const8new[1][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[end] + y[end][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
            end
            
            const6 = [const6;const6new]
            const7 = [const7;const7new]
            const8 = [const8;const8new]
            global TimePric += toq()
            tic()
            if maxcol == 5
             break
            end
            maxcol +=1
          end
        end #r in route
      end
  
      if !mincost
        if  result > globalincumbentvalue + 0.0001
         #prune  since a lower bound is already worse then incumbent 
         pos = findfirst(x -> x == current, Queue)
         deleteat!(Queue,pos)
         #println("Prune by Bound =", result, ", ", globalincumbentvalue)
         #read(STDIN,Char) 
         stop = true
         break
        end 
        if maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(dsol-floor.(dsol)) <= 0.00001
          #new integer solution
          #println("Integer Solution Found")

          ##############################################
          #Formulate Separation Problem and Solve 
          #Get First Stage results
          ssol = getvalue(s)
          usol= getvalue(u)
          for i in 1:length(resultz) #check for error and display
            if resultz[i] > 0.0001 && resultz[i] < 0.99
              println("SOLUCAO NAO INTEIRA!! ",i)
              read(STDIN,Char) 
            end
          end
      
          model2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000)) #CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=1000))

          @variable(model2, scen[1:length(DESTC)], Bin)
          @variable(model2, 0 <= y2[r in 1:length(routesused),i in 1:length(DESTC1), j in 1:length(DESTC1); i != j && resultz[r]> 0.001] <= 1) 
 
          @objective(model2, Min, ssol + sum(scen[i]*usol[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scen[i]*dsol[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y2[r,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused) if resultz[r]> 0.001 ) - sum(PRICEOD[i]*scen[i] for i in 1:length(DESTC)))

      
          @constraint(model2, const21[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], -(B2[routesused[r],i,j] + sum(scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*resultz[r] + y2[r,i,j] >= 1- scen[i] + 1 - scen[j] - length(DESTC))

          @constraint(model2, const22[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,length(DESTC1),i] >=  1 - scen[i] - length(DESTC) + 1)

          @constraint(model2, const23[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], -(0 + sum(scen[k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*resultz[r] + y2[r,i,length(DESTC1)] >= 1- scen[i]  - length(DESTC) + 1) 


          @constraint(model2, const22b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,length(DESTC1),i] <= (1-scen[i]) )
          @constraint(model2, const22c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,length(DESTC1),i] <= scen[j]*B2[routesused[r],j,i]+0+B2[routesused[r],i,j] + 1-B1[routesused[r],j] )

          @constraint(model2, const23b[r in 1:length(routesused),i in 1:length(DESTC); resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= (1-scen[i]) )
          @constraint(model2, const23c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC); i!=j && resultz[r] > 0.001], y2[r,i,length(DESTC1)] <= scen[j]*B2[routesused[r],i,j]+B2[routesused[r],j,i]+0 + 1-B1[routesused[r],j] )

          @constraint(model2, const21b[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[i]) )

          @constraint(model2, const21c[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC);i != j && resultz[r] > 0.001], y2[r,i,j] <= (1-scen[j])*B2[routesused[r],i,j] )

          @constraint(model2, const21d[r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC),k in 1:length(DESTC);i != j && k!=i && k!=j && resultz[r] > 0.001], y2[r,i,j] <= scen[k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k]  )




 
          #println("Will solve Separation Problem for integer solution")
           
          solve(model2)

 
          #println("solve model 2",getobjectivevalue(model2))
          if getobjectivevalue(model2) <= -0.0001 #-0.05
            #println(getvalue(scen), " get= ", getobjectivevalue(model2))
            #read(STDIN,Char)
            #Add new scenario and continue
            scenarionew = round.(getvalue(scen))
            #println("Integer Solution not valid. New scenario: ", scenarionew)
            push!(scenario,scenarionew)
            #Update Current node information on scenarios used
            push!(tree[current].addedcenarios,size(scenario,1))
            push!(cenariosused,size(scenario,1))
            #Create new variables
            ynew = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:1] for _ in 1:length(routesused)]
            for  r in 1:length(routesused),i in 1:length(DESTC1),j in 1:length(DESTC1)
              if i != j
                ynew[r][1][i][j] = @variable(model,  lowerbound = 0, upperbound =1)
              else
                ynew[r][1][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
              end
            end
            for r in 1:length(routesused)
             y[r] = [y[r];ynew[r]]
            end 
            #Create new constraints
            push!(const4, @constraint(model,  s + sum(scenario[end][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[end][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][end][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[end][i] for i in 1:length(DESTC)) ))


            
            const6new = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:1] for _ in 1:length(routesused)]
            const7new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:1] for _ in 1:length(routesused)]
            const8new = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:1] for _ in 1:length(routesused)]
            for  r in 1:length(routesused),i in 1:length(DESTC),j in 1:length(DESTC)
              if i != j
                const6new[r][1][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[end][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][end][i][j] >= 1- scenario[end][i] + 1 - scenario[end][j] - length(DESTC))
              else
                const6new[r][1][i][j] = @constraint(model, 0 == 0)
              end
            end
            for r in 1:length(routesused),i in 1:length(DESTC)
              const7new[r][1][i] = @constraint(model,  -(0 + sum(scenario[end][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][end][length(DESTC1)][i] >=  1 - scenario[end][i] - length(DESTC) + 1)
              const8new[r][1][i]= @constraint(model, -(0 + sum(scenario[end][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][end][i][length(DESTC1)] >= 1- scenario[end][i]  - length(DESTC) + 1)
            end
            for r in 1:length(routesused)
              const6[r]=[const6[r];const6new[r]]
              const7[r]=[const7[r];const7new[r]]
              const8[r]=[const8[r];const8new[r]]
            end
          else # end getobj < +.05
          #end verification of scenario insertion for probable integer solution
          ##############################################
          #Verify if it can be new incubent
            if result < globalincumbentvalue
              #println("New Incumbent Found")
              globalincumbentvalue = result
              globalincumbentsolution = resultz
            end 
            pos = findfirst(x -> x == current, Queue)
            deleteat!(Queue,pos)
            stop=true
            break
          end
          global TimeSce += toq()
          tic()
        elseif maximum(xsol-floor.(xsol)) > 0.00001
          if length(vardsolutionfixedzero) != 0 || length(vardsolutionfixedone) != 0
            println("d ja estava fixed e veio xij frac")
            println(vardsolutionfixedzero)
            println(vardsolutionfixedone)
            read(STDIN,Char)
          end 
          a = maximum(xsol-floor.(xsol))
          f=0
          g=0
          for i in 1:length(DESTC),j in 1:length(DESTC)
            if i != j
              if xsol[i,j]-floor.(xsol[i,j])== a
                f = i
                g = j
                break
              end
            end
          end
          #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
          b =[f;g] 
          push!(tree,TreeNode(current,[],[],[],[],[b],[],[]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))
          push!(tree,TreeNode(current,[],[],[],[b],[],[],[]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))

          pos = findfirst(x -> x == current, Queue)
          deleteat!(Queue,pos)

          current=Queue[end]
          global NODECOUNT += 1
          setlowerbound(x[f,g],0)
          setupperbound(x[f,g],0) 
        elseif maximum(dsol-floor.(dsol)) > 0.00001
          a = maximum(dsol-floor.(dsol))
          f=0
          for i in 1:length(DESTC)
            if dsol[i]-floor.(dsol[i])== a
              f = i
              break
            end
          end
          #println("Sol d still frac, so split into two nodes using max frac", f)  
          push!(tree,TreeNode(current,[],[],[],[],[],[],[f]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))
          push!(tree,TreeNode(current,[],[],[],[],[],[f],[]))
          push!(tree[current].children,length(tree))
          push!(Queue,length(tree))

          pos = findfirst(x -> x == current, Queue)
          deleteat!(Queue,pos)

          current=Queue[end]
          global NODECOUNT += 1
          setlowerbound(d[f],0)
          setupperbound(d[f],0) 
        else
          println("situation not covered")
          read(STDIN,Char)
        end
      end # !mincost
      stop && break
    end #end while true 'solve)
    return model, resultz 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate Prices, Routes and Scenarios to be used
  ###############################################
  #Create vectors for all possible scenarios
  global scenario = Vector{Vector{Int}}() 
  fillscenario(scenario,zeros(length(DESTC)),length(DESTC),1)
  println("All possible Scenario vectors created of size ", length(scenario) )
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenario)
    if sum(scenario[s]) <= CAPV
      groupe = []
      for c in 1:length(scenario[s])
        if scenario[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #################################################

  ###############################################
  #Re-Create vectors for nly initial scenarios - use col  cut generation afterwards
  global scenario = Vector{Vector{Int}}() 
  #Create  initial scenarios
  for i = 1:length(DESTC)
    SAMPLETEST = zeros(1:length(DESTC))
    SAMPLETEST[i]=1
    push!(scenario, SAMPLETEST)
  end
  push!(scenario, ones(1:length(DESTC)))
  push!(scenario, zeros(1:length(DESTC)))
  ################################################

  #Bypass REWOD and define new prices to pay ODs
  DESTC1=union(DESTC,[1])
  global PRICEOD
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance. 
  global PRICEOD = Inf*ones(1:length(DESTC))
  for i in DESTC
    for j in DESTC1
      for r in DESTC1
        if (j !=r && j != i && r != i) || (j==1 && r==1)
          if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
          end
        end
      end
    end
  end
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC),length(DESTC))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC)) 
  for r in 1:length(route)
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end 
  #####################################################################
 
  #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  routesused = []
  for r in 1:length(route)
    if length(route[r]) <= 1
       push!(routesused,r)
    end
  end
  #routesused = collect(1:length(route))
  scenariosused = collect(1:length(scenario))
  #Initialize branch-and-bound tree and queue
  Queue = Vector{Int}()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,scenariosused,[],[],[],[])) 
  push!(Queue,1)
  #Start processing queue in a Deep First mode

  global globalincumbentvalue = +Inf
  global globalincumbentsolution = []
    
  x = []
  master = []
  #always process last element of queue
  while length(Queue)>0
    current=Queue[end]
    (master, x) =process(current)
    global NODECOUNT += 1
    #if TimeMain +TimePric+TimeSce >= 10800
    #  break
    #end
  end
 
  for i in 1:length(globalincumbentsolution)
    if globalincumbentsolution[i] > 0.0001 && globalincumbentsolution[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  return globalincumbentvalue,globalincumbentsolution
end #end N3 function


#DDUEXACTN2: Same as N1, but now create our own B&B on variables first xij and then di. Extended formulation. Now decide compensation at master problem, list all routes and scenaros at once, no column generation or cut ...solve relaxed MILP...added B&B)
#correlated marginals, worst case approach, decision dependent probabilistic VRP; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL
#This time vehicle sohould satisfy capacity for all scenarios. Only routes within capacity are generated a priori and symetry eliminated (as in N1)
function solveSTODDUEXACTN2(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
  function process(current)
    #Initialize cenarios and solutions to be added
    cenariosused = []
    routesused = []
    varxsolutionfixedzero = []
    varxsolutionfixedone = [] 
    vardsolutionfixedzero = []
    vardsolutionfixedone = []  
    node = current
    while true
      cenariosused = [cenariosused;tree[node].addedcenarios]
      routesused = [routesused;tree[node].addedsolutions]
      varxsolutionfixedzero = [varxsolutionfixedzero;tree[node].varxfixedsolutionzero]
      varxsolutionfixedone = [varxsolutionfixedone;tree[node].varxfixedsolutionone]
      vardsolutionfixedzero = [vardsolutionfixedzero;tree[node].vardfixedsolutionzero]
      vardsolutionfixedone = [vardsolutionfixedone;tree[node].vardfixedsolutionone]
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
    #println("vardsolutionfixedzero =",vardsolutionfixedzero)
    #println("vardsolutionfixedone =",vardsolutionfixedone) 
    #read(STDIN,Char)

    result = 0
    resultz= []
    #Now formulate problem (relaxed problem now)

    z = Vector{Variable}(length(routesused))
    y = [[[Vector{Variable}(length(DESTC1)) for _ in 1:length(DESTC1)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const4 = Vector{ConstraintRef}(length(cenariosused))
    const6 = [[[Vector{ConstraintRef}(length(DESTC)) for _ in 1:length(DESTC)]  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const7 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]
    const8 = [[Vector{ConstraintRef}(length(DESTC))  for _ in 1:length(cenariosused)] for _ in 1:length(routesused)]


    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=7200))
    @variable(model,  s >= 0)
    #@variable(model,  z[r in routesused]>=0)
    #@variable(model,  y[r in routesused,w in cenariosused,i in 1:length(DESTC1),j in 1:length(DESTC1); i != j] >=0)

    for r in 1:length(routesused)
        z[r] = @variable(model, lowerbound = 0, upperbound =1)
    end
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC1),j in 1:length(DESTC1)
      if i != j
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =1)
      else
        y[r][w][i][j] = @variable(model,  lowerbound = 0, upperbound =0)
      end
    end
    @variable(model,  x[i in 1:length(DESTC),j in 1:length(DESTC); i != j]>=0)
    @variable(model,  u[1:length(DESTC)]<=0)
    @variable(model,  t[1:length(DESTC)]<=0)
    @variable(model,  d[1:length(DESTC)]>=0)

    for i in vardsolutionfixedzero
      setlowerbound(d[i],0)
      setupperbound(d[i],0) 
    end
    for i in vardsolutionfixedone
      setlowerbound(d[i],1)
      setupperbound(d[i],1)
    end

    for i in varxsolutionfixedzero
      setlowerbound(x[i[1],i[2]],0)
      setupperbound(x[i[1],i[2]],0) 
    end
    for i in varxsolutionfixedone
      setlowerbound(x[i[1],i[2]],1)
      setupperbound(x[i[1],i[2]],1)
    end  
 
    @objective(model, Min, s + sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) + sum(PROBOD2[i]*t[i] for i in 1:length(DESTC)) )

    @constraint(model, const1[i in 1:length(DESTC)], t[i]  >=  -10^4*d[i]) #be careful with cplex integrality tolerance
    @constraint(model, const2[i in 1:length(DESTC)], t[i] >= u[i] )

    for w in 1:length(cenariosused)
      const4[w] = @constraint(model,  s + sum(scenario[cenariosused[w]][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[cenariosused[w]][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r][w][i][j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in 1:length(routesused)) >= sum(PRICEOD[i]*scenario[cenariosused[w]][i] for i in 1:length(DESTC)))
    end

    @constraint(model, const5[i in 1:length(DESTC)], sum(B1[routesused[r],i]*z[r] for r in 1:length(routesused)) == 1)
    @constraint(model, const9[i in 1:length(DESTC)], d[i] <= 1)
    @constraint(model,  const10[i in 1:length(DESTC),j in 1:length(DESTC); i != j], x[i,j] == sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused)))
    
    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC),j in 1:length(DESTC)
      if i != j
        const6[r][w][i][j] = @constraint(model, -(B2[routesused[r],i,j] + sum(scenario[cenariosused[w]][k]*B3[routesused[r],i,j,k]+B4[routesused[r],i,j,k]+B5[routesused[r],i,j,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r][w][i][j] >= 1- scenario[cenariosused[w]][i] + 1 - scenario[cenariosused[w]][j] - length(DESTC))
      else
        const6[r][w][i][j] = @constraint(model, 0 == 0)
      end
    end

    for r in 1:length(routesused), w in 1:length(cenariosused),i in 1:length(DESTC)
      const7[r][w][i] = @constraint(model,  -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],k,i]+0+B2[routesused[r],i,k] + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][length(DESTC1)][i] >=  1 - scenario[cenariosused[w]][i] - length(DESTC) + 1)
      const8[r][w][i]= @constraint(model, -(0 + sum(scenario[cenariosused[w]][k]*B2[routesused[r],i,k]+B2[routesused[r],k,i]+0 + 1-B1[routesused[r],k] for k in 1:length(DESTC) if k != i))*z[r] + y[r][w][i][length(DESTC1)] >= 1- scenario[cenariosused[w]][i]  - length(DESTC) + 1)
    end

    while true
      status = solve(model)
      global TimeMain += toq()
      tic()
      #println("Solved Main problem ", getobjectivevalue(model), "; Time is ",TimeMain,"; Node is = ",current,"of ", NODECOUNT, " and length Queue is ", length(Queue))
      #read(STDIN,Char)
      #print(";",NODECOUNT,":",TimeMain,"...")
      result = getobjectivevalue(model)
      resultz = getvalue(z)
      #read(STDIN,Char)
      #println("d= ",getvalue(d))
      #read(STDIN,Char)
      #for r in 1:length(routesused)
      #  for w in 1:length(cenariosused)
      #    for i in 1:length(DESTC1)
      #      for j in 1:length(DESTC1)
      #        if i != j
      #          if getvalue(y[r][w][i][j])>0.0001
      #            #print("y [",r,",",w,",",i,",",j," ]:",getvalue(y[r,w,i,j]))
      #            #read(STDIN,Char)
      #          end
      #        end
      #      end
      #    end
      #  end 
      #end
      #println()
      #for r in 1:length(routesused)
      #  #println(r," ;",route[r])
      #  if getvalue(z[r]) >= 0.00001
      #    println(route[r],": z[ ",r," ]= ",getvalue(z[r])," ;") 
      #    read(STDIN,Char)
      #  end
      #end  #r 
      
      stop = false
      if status != :Optimal
        #prune   
        pos = findfirst(x -> x == current, Queue)
        deleteat!(Queue,pos)
        #println("Infeasible...prune ", status)
        #read(STDIN,Char) 
        stop = true
        break 
      end

      global globalincumbentvalue
      global globalincumbentsolution
      if  result > globalincumbentvalue + 0.0001
        #prune  since a lower bound is already worse then incumbent 
        pos = findfirst(x -> x == current, Queue)
        deleteat!(Queue,pos)
        #println("Prune by Bound =", result, ", ", globalincumbentvalue)
        #read(STDIN,Char) 
        stop = true
        break
      end 

      xsol1 = getvalue(x)
      xsol = zeros(length(DESTC),length(DESTC))
      for i=1:length(DESTC)
        for j=1:length(DESTC)
          if i != j
            if xsol1[i,j]<=0.00001 
              xsol[i,j] = 0
            elseif  xsol1[i,j] >=0.99999
              xsol[i,j]= 1
            else
              xsol[i,j]= xsol1[i,j]
              #print("xsol[",i," ,",j," ]",xsol[i,j])
            end
          end
        end
      end
      #println()
      #read(STDIN,Char)
      dsol = getvalue(d)
      for i=1:length(DESTC)
        if dsol[i]<=0.00001 
          dsol[i] = 0
        elseif  dsol[i] >=0.99999
          dsol[i]= 1
        end
      end

      if maximum(xsol-floor.(xsol)) <= 0.00001 && maximum(dsol-floor.(dsol)) <= 0.00001
        #new integer solution
        #println("Integer Solution Found")
        #Verify if it can be new incubent
        if result < globalincumbentvalue
          #println("New Incumbent Found")
          globalincumbentvalue = result
          globalincumbentsolution = resultz
        end 
        pos = findfirst(x -> x == current, Queue)
        deleteat!(Queue,pos)
        stop=true
      elseif maximum(xsol-floor.(xsol)) > 0.00001
        if length(vardsolutionfixedzero) != 0 || length(vardsolutionfixedone) != 0
            println("d ja estava fixed e veio xij frac")
            println(vardsolutionfixedzero)
            println(vardsolutionfixedone)
            read(STDIN,Char)
        end 
        a = maximum(xsol-floor.(xsol))
        f=0
        g=0
        for i in 1:length(DESTC),j in 1:length(DESTC)
          if i != j
            if xsol[i,j]-floor.(xsol[i,j])== a
             f = i
             g = j
             break
            end
          end
        end
        #println("Sol x still frac, so split into two nodes using max frac", f," ,",g) 
        b =[f;g] 
        push!(tree,TreeNode(current,[],[],[],[],[b],[],[]))
        push!(tree[current].children,length(tree))
        push!(Queue,length(tree))
        push!(tree,TreeNode(current,[],[],[],[b],[],[],[]))
        push!(tree[current].children,length(tree))
        push!(Queue,length(tree))

        pos = findfirst(x -> x == current, Queue)
        deleteat!(Queue,pos)

        current=Queue[end]
        global NODECOUNT += 1
        setlowerbound(x[f,g],0)
        setupperbound(x[f,g],0) 
      elseif maximum(dsol-floor.(dsol)) > 0.00001
        a = maximum(dsol-floor.(dsol))
        f=0
        for i in 1:length(DESTC)
          if dsol[i]-floor.(dsol[i])== a
            f = i
            break
          end
        end
        #println("Sol d still frac, so split into two nodes using max frac", f)  
        push!(tree,TreeNode(current,[],[],[],[],[],[],[f]))
        push!(tree[current].children,length(tree))
        push!(Queue,length(tree))
        push!(tree,TreeNode(current,[],[],[],[],[],[f],[]))
        push!(tree[current].children,length(tree))
        push!(Queue,length(tree))

        pos = findfirst(x -> x == current, Queue)
        deleteat!(Queue,pos)

        current=Queue[end]
        global NODECOUNT += 1
        setlowerbound(d[f],0)
        setupperbound(d[f],0) 


      
      else
        println("situation not covered")
        read(STDIN,Char)
      end
      stop && break
    end #end while true 'solve)
    return model, resultz 
  end #end Process function

  
  #########################################################################
  # MAIN PROGRAM 
  #Generate Prices, Routes and Scenarios to be used
  ###############################################
  #Create vectors for all possible scenarios
  global scenario = Vector{Vector{Int}}() 
  fillscenario(scenario,zeros(length(DESTC)),length(DESTC),1)
  println("Scenario vectors created of size ", length(scenario))
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenario)
    if sum(scenario[s]) <= CAPV
      groupe = []
      for c in 1:length(scenario[s])
        if scenario[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #################################################
  #Bypass REWOD and define new prices to pay ODs
  DESTC1=union(DESTC,[1])
  global PRICEOD
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance. 
  global PRICEOD = Inf*ones(1:length(DESTC))
  for i in DESTC
    for j in DESTC1
      for r in DESTC1
        if (j !=r && j != i && r != i) || (j==1 && r==1)
          if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
          end
        end
      end
    end
  end
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC),length(DESTC))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC)) 
  for r in 1:length(route)
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end 
  #####################################################################
 
  #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  #routesused = []
  #for r in 1:length(route)
  #  if length(route[r]) <= 1
  #     push!(routesused,r)
  #  end
  #end
  routesused = collect(1:length(route))
  scenariosused = collect(1:length(scenario))
  #Initialize branch-and-bound tree and queue
  Queue = Vector{Int}()
  tree = Vector{TreeNode}()
  push!(tree,TreeNode(0,[],routesused,scenariosused,[],[],[],[])) 
  push!(Queue,1)
  #Start processing queue in a Deep First mode

  global globalincumbentvalue = +Inf
  global globalincumbentsolution = []
    
  x = []
  master = []
  #always process last element of queue
  while length(Queue)>0
    current=Queue[end]
    (master, x) =process(current)
    global NODECOUNT += 1
    if TimeMain +TimePric+TimeSce >= 10800
      break
    end
  end
  
  for i in 1:length(globalincumbentsolution)
    if globalincumbentsolution[i] > 0.0001 && globalincumbentsolution[i] < 0.99
      println("SOLUCAO NAO INTEIRA!! ",i)
      read(STDIN,Char) 
    end
  end
  
  return globalincumbentvalue,globalincumbentsolution
end #end N2 function


#DDUEXACTN1: Extended formulation. Now decide compensation at master problem, list all routes and scenaros at once, no column generation or cut ...solve MILP...no added B&B)
#correlated marginals, worst case approach, decision dependent probabilistic VRP; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL
#This time vehicle sohould satisfy capacity for all scenarios. Only routes within capacity are generated a priori
function solveSTODDUEXACTN1(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
 
  
  #########################################################################
  # MAIN PROGRAM 
  #Generate Prices, Routes and Scenarios to be used
  ###############################################
  #Create vectors for all possible scenarios
  global scenario = Vector{Vector{Int}}() 
  fillscenario(scenario,zeros(length(DESTC)),length(DESTC),1)
  println("Scenario vectors created of size ", length(scenario))
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenario)
    if sum(scenario[s]) <= CAPV
      groupe = []
      for c in 1:length(scenario[s])
        if scenario[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #################################################
  #Bypass REWOD and define new prices to pay ODs
  DESTC1=union(DESTC,[1])
  global PRICEOD
  global PRICEOD2
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance.
  global inst_idf
  #println(inst_idf[end])
  #read(STDIN,Char)
  if inst_idf[end]== 'A'
    println("Instance A for PRICE")
    PRICEOD2 = 1.0*zeros(1:length(DESTC))
    PRICEOD = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  elseif inst_idf[end]== 'B'
    println("Instance B for PRICE")
    PRICEOD = 1.0*zeros(1:length(DESTC))
    PRICEOD2 = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD2[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD2[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD2[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
  else  #if inst_id[end]= "C" or others
    println("Instance C for PRICE")
    PRICEOD = Inf*ones(1:length(DESTC))
    for i in DESTC
      for j in DESTC1
        for r in DESTC1
          if (j !=r && j != i && r != i) || (j==1 && r==1)
            if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
              PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
            end
          end
        end
      end
    end
    PRICEOD2 = PRICEOD
  end #end inst_id
  #read(STDIN,Char)
  ################################################
  #Calculate needed parameters for routes and store:
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC),length(DESTC))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B6 = zeros(length(route),length(DESTC),length(DESTC)) 
  for r in 1:length(route)
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        if i != route[r][end]
          B6[r,i,route[r][findfirst(x -> x==i, route[r])+1]] = 1
        end
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end 
  #####################################################################
 
    #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  #routesused = []
  #for r in 1:length(route)
  #  if length(route[r]) <= 1
  #     push!(routesused,r)
  #  end
  #end
 
  routesused = collect(1:length(route))

  result = 0
  resultz= []
  #Now formulate problem (relaxed problem now)
   
  model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=7200,CPX_PARAM_MIPDISPLAY=0)) #CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=3600
  @variable(model,  s >= 0)
  @variable(model,  z[r in routesused], Bin)
  @variable(model,  y[r in routesused,1:length(scenario),i in 1:length(DESTC1),j in 1:length(DESTC1); i != j]>=0)
  @variable(model,  u[1:length(DESTC)]<=0)
  @variable(model,  t[1:length(DESTC)]<=0)
  @variable(model,  d[1:length(DESTC)], Bin) 
  #@variable(model,  x[i in 1:length(DESTC),j in 1:length(DESTC); i != j],Bin)

  if inst_idf[end]== 'A' #probod2 = 0, so no d variabels
      for i in 1:length(DESTC)
        setlowerbound(d[i],0)
        setupperbound(d[i],0)
      end
  end

 
  @objective(model, Min, s + sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) + sum(PROBOD2[i]*t[i] for i in 1:length(DESTC)) )

  @constraint(model, const1[i in 1:length(DESTC)], t[i]  >=  -10^4*d[i]) #be careful with cplex integrality tolerance
  @constraint(model, const2[i in 1:length(DESTC)], t[i] >= u[i] )
  
  @constraint(model, const4[w in 1:length(scenario)], s + sum(scenario[w][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD2[i]*scenario[w][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r,w,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in routesused) >= sum(PRICEOD[i]*scenario[w][i] for i in 1:length(DESTC)))

  #@constraint(model, const4[w in 1:length(scenario)], s + sum(scenario[w][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[w][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r,w,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in routesused) >= 0)

  @constraint(model, const5[i in 1:length(DESTC)], sum(B1[r,i]*z[r] for r in routesused) == 1)

  #@constraint(model, const9[r in routesused], sum(B1[r,i]*z[r] for i in 1:length(DESTC)) <= CAPV)
   #@constraint(model, const11, sum(z[r] for r in routesused) <= REGV)
  #@constraint(model, const9[i in 1:length(DESTC)], d[i] <= 1)

  #@constraint(model,  const10[i in 1:length(DESTC),j in 1:length(DESTC); i != j], x[i,j] == sum(B6[routesused[r],i,j]*z[r] for r in 1:length(routesused)))
  
  @constraint(model, const6[r in routesused,w in 1:length(scenario), i in 1:length(DESTC),j in 1:length(DESTC); i != j], -(B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r,w,i,j] >= 1- scenario[w][i] + 1 - scenario[w][j] - length(DESTC))
  
  @constraint(model, const7[r in routesused,w in 1:length(scenario),j in 1:length(DESTC)], -(0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j))*z[r] + y[r,w,length(DESTC1),j] >=  1 - scenario[w][j] - length(DESTC) + 1)

  @constraint(model, const8[r in routesused,w in 1:length(scenario), i in 1:length(DESTC)], -(0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[r] + y[r,w,i,length(DESTC1)] >= 1- scenario[w][i]  - length(DESTC) + 1)

  status = solve(model)
  global TimeMain += toq()
  println("Solved Main problem ", getobjectivevalue(model), "; Time is ",TimeMain)
  result = getobjectivevalue(model)
  resultz = getvalue(z)
  dsol = getvalue(d)
  tic()
  #read(STDIN,Char)
  #println("d= ",getvalue(d))
  println("number of routes=  ",sum(resultz))
  #for r in routesused
  #  if resultz[r] > 0.01
  #    println("route=  ",route[r])
  #  end
  #end
  #read(STDIN,Char)
  #for r in routesused
  #  for w in 1:length(scenario)
  #    for i in 1:length(DESTC1)
  #      for j in 1:length(DESTC1)
  #        if i != j
  #          if getvalue(y[r,w,i,j])>0.0001
  #            #print("y [",r,",",w,",",i,",",j," ]:",getvalue(y[r,w,i,j]))
  #            #read(STDIN,Char)
  #          end
  #        end
  #      end
  #    end
  #  end 
  #end
  #println()
  #for r in routesused
  #  #println(r," ;",route[r])
  #  if getvalue(z[r]) != 0
  #    println(route[r],": z[ ",r," ]= ",getvalue(z[r])," ;") 
  #    read(STDIN,Char)
  #  end
  #end  #r 
  global NODECOUNT = MathProgBase.getnodecount(model)
  #println(MathProgBase.getnodecount(model)) 
  global NROUTES = sum(resultz)
  global PERHIGH = sum(dsol)/length(dsol) 
return result,resultz 
end #end N1 function


#DDUEXACTV3: Extended formulation. Now decide compensation at master problem, using dynamic programming to calculate routes (in this first moment just list all routes at once, no column generation)
#correlated marginals, worst case approach, decision dependent probabilistic tsp; recourse skips outsourced customers and no detour, BOOLEAN COMPENSATION LEVEL, extended formulation with exact column generation
function solveSTODDUEXACTV3(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
 
  
  #########################################################################
  # MAIN PROGRAM for function V3
  #Generate Prices, Routes and Scenarios to be used
  #Will also generate all costs at this moment
  ###############################################
  #Create vectors for all possible scenarios
  global scenario = Vector{Vector{Int}}() 
  fillscenario(scenario,zeros(length(DESTC)),length(DESTC),1)
  println("Scenario vectors created of size ", length(scenario))
  ################################################
  #compensation = []   #This one in case we want to run with just one compensation
  #push!(compensation,[0,0,0,1,0,1,0,0])
  #push!(compensation,[0,1,1,1,0,0,0,0])
  #push!(compensation,zeros(1:length(DESTC)))
  #println(compensation)
  #println(length(compensation))
  #read(STDIN,Char)
  ###############################################
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  for s in 1:length(scenario)
    if sum(scenario[s]) <= CAPV
      groupe = []
      for c in 1:length(scenario[s])
        if scenario[s][c] !=0
          push!(groupe,c)
        end
      end
      if groupe != []
        fillroute(groupe, 1, length(groupe),route)  #symetry breaking included
        #fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
      end
    end
  end
  println("Route vectors created of size ", length(route))
  #for r in 1:length(route)
  #  if length(route[r]) == 3 
  #    println(route[r])
  #    read(STDIN,Char)
  #  end
  #end
  ################################################
  println(length(route))
  #################################################
  #Bypass REWOD and define new prices to pay ODs
  #XENIA: Note that at this moment I am not using the compensation fee (rewod and rewod2) that comes from the data. I am creating my own compensation fee, and I pay one time that (low level) or twice times that (high level).  See why below...
  DESTC1=union(DESTC,[1])
  global PRICEOD
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance. 
  #XENIA: To start with 
  global PRICEOD = Inf*ones(1:length(DESTC))
  for i in DESTC
    for j in DESTC1
      for r in DESTC1
        if (j !=r && j != i && r != i) || (j==1 && r==1)
          if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
          end
        end
      end
    end
  end
  #PRICEOD[1] =0
  #PRICEOD[2] =10
  #PRICEOD[3] =0
  ################################################
  #Calculate needed parameters for routes and store:
  #B1[r,i]=1, if customer i is in route r; 0 otherwise
  #B2[r,i,j]=1, if cust i is served before cust j is in route r; 0 otherwise
  #B3[r,i,j,r]=1 if cust r is between cust i and j in route r; 0 otherwise
  #B4[r,i,j,r]=1 if cust r is before cust i and j in route r; 0 otherwise
  #B5[r,i,j,r]=1 if cust r is after cust i and j in route r; 0 otherwise
  global B1 = zeros(length(route),length(DESTC))
  global B2 = zeros(length(route),length(DESTC),length(DESTC))
  global B3 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B4 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC))
  global B5 = zeros(length(route),length(DESTC),length(DESTC),length(DESTC)) 


  for r in 1:length(route)
    #println("route= ",route[r])
    #read(STDIN,Char)
    for i in 1:length(DESTC)
      if findfirst(x -> x==i, route[r]) != 0
        B1[r,i] = 1
        #println("B1[ ",r,",",i," ]= ", B1[r,i])
        #read(STDIN,Char)
      end 
      for j in 1:length(DESTC)
        if  findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && (findfirst(x -> x==j, route[r]) > findfirst(x -> x==i, route[r]))
           B2[r,i,j] = 1
           #print("B2[ ",r," ,",i," ,",j," ]= ",B2[r,i,j], "; ")
           #read(STDIN,Char)
        end
        #println()
        for k in 1:length(DESTC)
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B3[r,i,j,k] = 1
            #print("B3[ ",r," ,",i,",",j," ,",k," ]= ",B3[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) < findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) < findfirst(x -> x==j, route[r])
            B4[r,i,j,k] = 1
            #print("B4[ ",r," ,",i,",",j," ,",k," ]= ",B4[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
          if i!=j && j!=k && k!=i &&findfirst(x -> x==i, route[r]) != 0 && findfirst(x -> x==j, route[r]) != 0 && findfirst(x -> x==k, route[r]) != 0 && findfirst(x -> x==k, route[r]) > findfirst(x -> x==i, route[r]) && findfirst(x -> x==k, route[r]) > findfirst(x -> x==j, route[r]) 
            B5[r,i,j,k] = 1
            #print("B5[ ",r," ,",i,",",j," ,",k," ]= ",B5[r,i,j,k]," ;")
            #read(STDIN,Char)
          end
          #println()
        end
      end
    end
  end 
  #####################################################################
 
    #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  routesused = []
  for r in 1:length(route)
    if length(route[r]) <= 1
       push!(routesused,r)
    end
  end
 
  scenarioused = [1:length(scenario)]
 

  result = 0
  resultz= []
  while true
  #Now formulate problem (relaxed problem now)
   
  model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=3600))
  @variable(model,  s >= 0)
  #@variable(model,  z[r in routesused], Bin)
  @variable(model,  z[r in routesused]>=0)
  @variable(model,  y[r in routesused,1:length(scenario),i in 1:length(DESTC1),j in 1:length(DESTC1); i != j] >= 0)
  @variable(model,  u[1:length(DESTC)]<=0)
  @variable(model,  t[1:length(DESTC)]<=0)
  #@variable(model,  d[1:length(DESTC)], Bin)
  @variable(model,  d[1:length(DESTC)]>=0)  
 
  @objective(model, Min, s + sum(PROBOD[i]*u[i] for i in 1:length(DESTC)) + sum(PROBOD2[i]*t[i] for i in 1:length(DESTC)) )

  @constraint(model, const1[i in 1:length(DESTC)], t[i]  >=  -10^4*d[i]) #be careful with cplex integrality tolerance
  @constraint(model, const2[i in 1:length(DESTC)], t[i] >= u[i] )
  
  @constraint(model, const4[w in 1:length(scenario)], s + sum(scenario[w][i]*u[i] for i in 1:length(DESTC)) - sum(PRICEOD[i]*scenario[w][i]*d[i] for i in 1:length(DESTC)) - sum(sum(sum(REWV*A[DESTC1[i],DESTC1[j]]*y[r,w,i,j] for i in 1:length(DESTC1) if i != j) for j in 1:length(DESTC1)) for r in routesused) >= sum(PRICEOD[i]*scenario[w][i] for i in 1:length(DESTC)))

  @constraint(model, const5[i in 1:length(DESTC)], sum(B1[r,i]*z[r] for r in routesused) == 1)

  #@constraint(model, const9[r in routesused], sum(B1[r,i]*z[r] for i in 1:length(DESTC)) <= CAPV)
  @constraint(model, const9[i in 1:length(DESTC)], d[i] <= 1)
  
  @constraint(model, const6[r in routesused,w in 1:length(scenario), i in 1:length(DESTC),j in 1:length(DESTC); i != j], -(B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*z[r] + y[r,w,i,j] >= 1- scenario[w][i] + 1 - scenario[w][j] - length(DESTC))
  
  @constraint(model, const7[r in routesused,w in 1:length(scenario),j in 1:length(DESTC)], -(0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j))*z[r] + y[r,w,length(DESTC1),j] >=  1 - scenario[w][j] - length(DESTC) + 1)

  @constraint(model, const8[r in routesused,w in 1:length(scenario), i in 1:length(DESTC)], -(0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*z[r] + y[r,w,i,length(DESTC1)] >= 1- scenario[w][i]  - length(DESTC) + 1)

  #@constraint(model, const10[r in routesused,w in 1:length(scenario), i in 1:length(DESTC1),j in 1:length(DESTC1); i != j], y[r,w,i,j] <= 1)

  

  status = solve(model)
  global TimeMain += toq()
  println("Solved Main problem ", getobjectivevalue(model), "; Time is",TimeMain)
  result = getobjectivevalue(model)
  resultz = getvalue(z)
  tic()
    #read(STDIN,Char)
    #println("d= ",getvalue(d))
    #read(STDIN,Char) 
    #println("dual= ",getdual(const4))
    #read(STDIN,Char)
    #print("y= ",getvalue(y))
    #read(STDIN,Char) 
    #println()
    #read(STDIN,Char)
    #println("d*t= ",getvalue(t).*getvalue(d))
    #read(STDIN,Char)
    #println("u= ",getvalue(u))
    #read(STDIN,Char) 
  gapsvalue = []
  for r in routesused
      #println(r," ;",route[r])
      #if getvalue(z[r]) >= 1
      #  println(route[r],": z[ ",r," ]= ",getvalue(z[r])," ;") 
      #  read(STDIN,Char)
      #end
    for w in 1:length(scenario)
    # #println("scenario = ", scenario[w])
      for i in 1:length(DESTC1)
        for j in 1:length(DESTC1)
          if i != j
    #       if i != length(DESTC1) && j != length(DESTC1)
    #           #if B2[r,i,j] != 0
    #           #println("B2[ ",r,",",i,",",j," ]=", B2[r,i,j])
    #           #read(STDIN,Char)
    #           #end
    #           end
    #           if getvalue(y[r,w,i,j]) > 0 #&& getvalue(z[r]) > 3 
    #             println("y[ ",r,",",w,",",i,",",j,"]= ",getvalue(y[r,w,i,j]),";",getvalue(z[r])," ;" )
            if DESTC1[i] == 1
                   #println("F[ ",r,",",w,",1,",j,"]= ",getdual(const7[r,w,j]) )
                   #println("gap= ", (0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j))* getvalue(z[r]) +  1 - scenario[w][j] - length(DESTC) + 1, " K1=",   1 - scenario[w][j] - length(DESTC) + 1," K2=", (0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j))," z= ", getvalue(z[r]))
               if !in((0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j)), gapsvalue)
                 push!(gapsvalue,(0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j)))
               end
            elseif DESTC1[j] == 1
                   #println("F[ ",r,",",w,",",i,",1]= ",getdual(const8[r,w,i]) )
                   #println("gap= ", 1- scenario[w][i]  - length(DESTC) + 1 +(0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*getvalue(z[r]), " K2=",(0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))," K1=",1- scenario[w][i]  - length(DESTC) + 1," z= ", getvalue(z[r]))
               if !in((0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i)), gapsvalue)
                 push!(gapsvalue,(0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i)))
               end
            else
                   #println("F[ ",r,",",w,",",i,",",j,"]= ",getdual(const6[r,w,i,j]) )
                   #println("gap= ", 1- scenario[w][i] + 1 - scenario[w][j] - length(DESTC) +(B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*getvalue(z[r]), " K2=",(B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i)), " K1= ", 1- scenario[w][i] + 1 - scenario[w][j] - length(DESTC)," z= ", getvalue(z[r]))
               if !in((B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC)  if k != j && k != i)), gapsvalue)
                 push!(gapsvalue,(B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i)))
               end
            end
          end
        end
      end
    end  #w
  end  #r
  #println(gapsvalue)
  E1dual = getdual(const5)
  println(getdual(const5)) 
  read(STDIN,Char) 
  #Formulate and solve dual to check properties
  model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=3600))
  @variable(model,  D1[1:length(scenario)]>=0)
  @variable(model,  F1[r in routesused,1:length(scenario),i in 1:length(DESTC1),j in 1:length(DESTC1); i != j] >= 0)
  @variable(model,  A1[1:length(DESTC)]>=0)
  @variable(model,  G1[1:length(DESTC)]>=0)
  @variable(model,  C1[1:length(DESTC)]<=0)
  @variable(model,  E1[1:length(DESTC)])    
 
  @objective(model, Max, sum(C1[i] for i in 1:length(DESTC)) + sum( (sum( PRICEOD[i]*scenario[w][i] for i in 1:length(DESTC)))*D1[w]for w in 1:length(scenario)) + sum(E1[i] for i in 1:length(DESTC)) + sum(sum(sum(sum( (1-scenario[w][i]+1-scenario[w][j]-length(DESTC))*F1[r,w,i,j]for i in 1:length(DESTC) if i !=j) for j in 1:length(DESTC) ) for w in 1:length(scenario)) for r in routesused) + sum(sum(sum( (1-scenario[w][i]+1-length(DESTC))*F1[r,w,i,length(DESTC1)] for i in 1:length(DESTC))  for w in 1:length(scenario)) for r in routesused) + sum(sum(sum( (1+1-scenario[w][j]-length(DESTC))*F1[r,w,length(DESTC1),j] for j in 1:length(DESTC) ) for w in 1:length(scenario)) for r in routesused) )

  @constraint(model, const1, sum(D1[w] for w in 1:length(scenario)) <= 1)
  @constraint(model, const2[i in 1:length(DESTC)], -G1[i] + sum(scenario[w][i]*D1[w] for w in 1:length(scenario)) >= PROBOD[i])

  @constraint(model, const3[i in 1:length(DESTC)], A1[i] + G1[i] >= PROBOD2[i] )
  #@constraint(model, const7[i in 1:length(DESTC)], E1[i] == E1dual[i] ) 

  @constraint(model, const4[i in 1:length(DESTC)], 10^4*A1[i] + C1[i] - sum(PRICEOD[i]*scenario[w][i]*D1[w] for w in 1:length(scenario))<=0 ) 

  @constraint(model, const5[r in routesused], sum(B1[r,i]*E1[i] for i in 1:length(DESTC))- sum(sum(sum( (B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*F1[r,w,i,j]for i in 1:length(DESTC) if i !=j) for j in 1:length(DESTC) ) for w in 1:length(scenario))  - sum(sum( (0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*F1[r,w,i,length(DESTC1)] for i in 1:length(DESTC))  for w in 1:length(scenario)) - sum(sum( (0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j))*F1[r,w,length(DESTC1),j] for j in 1:length(DESTC) ) for w in 1:length(scenario))<=0 )  

  @constraint(model, const6[r in routesused,w in 1:length(scenario),i in 1:length(DESTC1),j in 1:length(DESTC1); i != j ],  F1[r,w,i,j]<=A[DESTC1[i],DESTC1[j]]*D1[w])

  #status = solve(model)
  Timetemp = toq()
  global TimeMain += Timetemp
  #println("Solved Dual problem ", getobjectivevalue(model), "; Time is",Timetemp)
  tic() 

  #println("D1 =", getvalue(D1))
  #read(STDIN,Char)
  #println("E1 =", getvalue(E1))
  #read(STDIN,Char)
  ##println("F1 =", getvalue(F1))
  ##read(STDIN,Char)
 
  #println("comecou")
  for r in routesused
    for w in 1:length(scenario)
      for i in 1:length(DESTC1)
        for j in 1:length(DESTC1)
          if i != j
             #if getvalue(F1[r,w,i,j]) != 0
               #println("F1[ ",r,",",w,",",i,",",j, "]=",getvalue(F1[r,w,i,j]))
               #read(STDIN,Char)
             #end
          end
        end
      end
    end
  end
  #println("terminou")

  #checke to see if all routes feasible and optimal does not change
  stop = true
  for r in 1:length(route)
    if !in(r,routesused)
    if sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)) != 0
     println("vou inserir route ", r, " valor constraint= ", sum( B1[r,i]*E1dual[i] for i in 1:length(DESTC)))
     push!(routesused,r)
     stop = false
     break
     read(STDIN,Char)
    end
    end
  end
  println("size routes used = ", length(routesused))
  stop && break
  end  #end while true
    #if sum(B1[r,i]*getvalue(E1[i]) for i in 1:length(DESTC)) !=  sum(sum(sum( (B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*getvalue(F1[r,w,i,j]) for i in 1:length(DESTC) if i !=j) for j in 1:length(DESTC) ) for w in 1:length(scenario))  + sum(sum( (0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*getvalue(F1[r,w,i,length(DESTC1)])for i in 1:length(DESTC))  for w in 1:length(scenario)) + sum(sum( (0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j))*getvalue(F1[r,w,length(DESTC1),j]) for j in 1:length(DESTC) ) for w in 1:length(scenario)) 
    
     # println("ERROR route ",r," :",sum(B1[r,i]*getvalue(E1[i]) for i in 1:length(DESTC)) ," ; ",  sum(sum(sum( (B2[r,i,j] + sum(scenario[w][k]*B3[r,i,j,k]+B4[r,i,j,k]+B5[r,i,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j && k != i))*getvalue(F1[r,w,i,j]) for i in 1:length(DESTC) if i !=j) for j in 1:length(DESTC) ) for w in 1:length(scenario))  + sum(sum( (0 + sum(scenario[w][k]*B2[r,i,k]+B2[r,k,i]+0 + 1-B1[r,k] for k in 1:length(DESTC) if k != i))*getvalue(F1[r,w,i,length(DESTC1)])for i in 1:length(DESTC))  for w in 1:length(scenario)) + sum(sum( (0 + sum(scenario[w][k]*B2[r,k,j]+0+B2[r,j,k] + 1-B1[r,k] for k in 1:length(DESTC) if k != j))*getvalue(F1[r,w,length(DESTC1),j]) for j in 1:length(DESTC) ) for w in 1:length(scenario)) )
             
       
  read(STDIN,Char)    
return result,resultz 




end

#DDUEXACTV2: Generate routes and compensations (options) by column generation
#DDUEXACTV2: independent, decision dependent probabilistic tsp; recourse skips outsourced customers plus detour, BOOLEAN COMPENSATION LEVEL, extended formulation with exact column generation
function solveSTODDUEXACTV2(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
  
  function detour(custpos,orderedcust,PROBCB,PROBC1B) #function used as part of calculation of average cost of a given route and compensation. This function calculates the cost of detour to depot and back to next available customer
    #PROBCB outsourced
    #PROBC1B  not outsourced
    function ff(m,r,orderedcust,PROBCB,PROBC1B) #exactly r customers among the 1,..,m customers are present
      if m == r
        return prod(PROBC1B[orderedcust[1:m]])
      elseif r==0
        return prod(PROBCB[orderedcust[1:m]]) 
      else
        tot = 0
        tot +=PROBC1B[orderedcust[m]]*ff(m-1,r-1,orderedcust,PROBCB,PROBC1B)
        tot +=PROBCB[orderedcust[m]]*ff(m-1,r,orderedcust,PROBCB,PROBC1B)
        return tot
      end
    end

    if custpos<=CAPV-1
      return 0
    else
      summ = 0
      for k in 1:floor(custpos/CAPV)
        summ += ff(custpos-1,k*CAPV-1,orderedcust,PROBCB,PROBC1B) 
      end
      return summ
    end
  end
  
  #########################################################################
  # MAIN PROGRAM for function V2
  #Generate Prices, Routes and Compensations to be used
  #Will also generate all costs at this moment
  ###############################################
  #Create vectors for all possible compensations
  global compensation = Vector{Vector{Int}}() 
  fillscenario(compensation,zeros(length(DESTC)),length(DESTC),1)
  println("Compensations vector created")
  ################################################
  #compensation = []   #This one in case we want to run with just one compensation
  #push!(compensation,[0,0,0,1,0,1,0,0])
  #push!(compensation,[0,1,1,1,0,0,0,0])
  #push!(compensation,zeros(1:length(DESTC)))
  #println(compensation)
  #println(length(compensation))
  #read(STDIN,Char)
  ###############################################
  #Create vectors for all possible routings
  global route = Vector{Vector{Int}}() 
  fillroute(collect(1:length(DESTC)), 1, length(DESTC),route)
  println("Route vectors created of size ", length(route))
  ################################################
  #println(length(route))
  #################################################
  #Bypass REWOD and define new prices to pay ODs
  #XENIA: Note that at this moment I am not using the compensation fee (rewod and rewod2) that comes from the data. I am creating my own compensation fee, and I pay one time that (low level) or twice times that (high level).  See why below...
  DESTC1=union(DESTC,[1])
  global PRICEOD
  #Calculate minimum price to pay for ODs
  #This is a critical point. If the OD is available, the model always pays the defined price (or     compensation fee) to the OD, even if the solution is not optimal. To mitigate sub-optimization, we define a minimal price for each customer below. Note that it coul be zero in the case one customer is located exactly in between the straight line of two other customers. So we add a fixed ammount at the end. Since the probabilities are already given in the isntance, we just assume that the prices defined below are coherent with the probabilities defined in the instance. 
  #XENIA: To start with 
  global PRICEOD = Inf*ones(1:length(DESTC))
  for i in DESTC
    for j in DESTC1
      for r in DESTC1
        if (j !=r && j != i && r != i) || (j==1 && r==1)
          if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
          end
        end
      end
    end
  end
  #PRICEOD[1] =0
  #PRICEOD[2] =10
  #PRICEOD[3] =0
  ################################################
  #Not yet Calculate costs
  global cost = zeros(length(route),length(compensation))
  
  #####################################################################
 
    #Initialize routes and compensations to be added. With the lack of better heuristic for incubent solutions just defined ad-hoc routes and compensations now 
  routesused = [1,2,3]
  compensationused = [1]

  #calculate here upfront the cost coeficient of all feasible solutions and existent scenarios at this time
  #XENIA; Here is how I calculate the average cost of a given route and compenstion vector
  global cost
  for r in routesused #1:length(route) 
    for c in compensationused
      if cost[r,c] == 0
        #Define compensation and probability Vector 
        #CM = REWOD + REWOD2 .* compensation[c]
        CM = PRICEOD + 10*PRICEOD .* compensation[c]
        PM1 = PROBOD + PROBOD2  .* compensation[c]  #outsourced
        #CM = 0 + PRICEOD .* compensation[c]
        #PM1 = 0 + PROBOD  .* compensation[c]  #outsourced
        PM = ones(1:length(DESTC)).-PM1    #not outsourced
        
        ####################################################
        #Define Cost Vector
        #############################################################################
        LowerIndenpright=0 
        for i in 1:length(route[r])
          cost[r,c] += CM[route[r][i]]*PM1[route[r][i]]
          cost[r,c] += REWV*A[1,DESTC[route[r][i]]]*PM[route[r][i]]*prod(PM1[route[r][1:i-1]])
          cost[r,c] += REWV*A[DESTC[route[r][i]],1]*PM[route[r][i]]*prod(PM1[route[r][i+1:end]])    
          for j in i+1:length(route[r])
            cost[r,c] += REWV*A[DESTC[route[r][i]],DESTC[route[r][j]]]*PM[route[r][i]]*PM[route[r][j]]*prod(PM1[route[r][i+1:j-1]])
            LowerIndenpright += REWV*(A[DESTC[route[r][i]],1]+A[1,DESTC[route[r][j]]]-A[DESTC[route[r][i]],DESTC[route[r][j]]]) *PM[route[r][i]]*PM[route[r][j]]*prod(PM1[route[r][i+1:j-1]])*detour(i,route[r],PM1,PM)
          end
        end
        #############################################################################
        LowerIndenpInv=0
        routeinv=reverse(route[r])
        for i in 1:length(routeinv)
          for j in i+1:length(routeinv)
            LowerIndenpInv += REWV*(A[DESTC[routeinv[i]],1]+A[1,DESTC[routeinv[j]]]-A[DESTC[routeinv[i]],DESTC[routeinv[j]]]) *PM[routeinv[i]]*PM[routeinv[j]]*prod(PM1[routeinv[i+1:j-1]])*detour(i,routeinv,PM1,PM)
          end
        end
        cost[r,c] += min(LowerIndenpright,LowerIndenpInv)
        #print("r= ",r, " c= ", c," , ", cost[r,c], " ;;")
        #read(STDIN,Char)
      end #end cost[r,c] == 0
    end
  end


  #Now formulate problem  
  y = Vector{Variable}(length(routesused))  #define vector for y variables
  z = Matrix{Variable}(length(routesused),length(compensationused)) #define matrix for z variables
  const3 = Vector{ConstraintRef}(length(routesused))  #define vector for constraints
   
  model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=3600))
   
  for r in 1:length(routesused)
    y[r] = @variable(model,  lowerbound = 0, upperbound =1)
  end
  for r in 1:length(routesused), c in 1:length(compensationused)
    z[r,c] = @variable(model,  lowerbound = 0, upperbound =1)
  end

  @objective(model, Min, sum(sum( cost[routesused[r],compensationused[c]]*z[r,c] for r=1:length(routesused)) for c=1:length(compensationused) ))

  @constraint(model, const2, sum( y[r] for r=1:length(routesused))== 1)
   
  for r in 1:length(routesused)
    const3[r]=@constraint(model, sum(z[r,c] for c=1:length(compensationused)) - y[r] == 0)
  end
  stop = false
  while true
    status = solve(model)
    global TimeMain += toq()
    tic()
    trueobjective = getobjectivevalue(model)
    println("Solved node ",  trueobjective, "; ",-getdual(const2))
    #read(STDIN,Char)  
    if status != :Optimal  #If Not Optimal stop because there is error
      println("Stopped...Not Optimal! ")
      read(STDIN,Char)  
      stop = true
      break 
    end
   
    #Do pricing now
    const2_dual = -getdual(const2)
    #const3_dual = [-getdual(const3[r]) for r in 1:length(routesused)]

    #Will consider introducing route columns and  compensation columns
    #First route column              
    routereducedcosts = +Inf #0
    best_r = 0 
    best_c = 0       
    for r in 1:length(route)  #test reduced cost for all new routes 
      routereducedcosts = +Inf #0
      best_r = r
      for c in 1:length(compensation)
        global cost
        if cost[r,c] == 0
          #Define compensation and probability Vector 
          #CM = REWOD + REWOD2 .* compensation[c]
          #CM = PRICEOD + 10*PRICEOD .* compensation[c]
          #PM1 = PROBOD + PROBOD2  .* compensation[c]  #outsourced
          CM = 0 + PRICEOD .* compensation[c]
          PM1 = 0 + PROBOD  .* compensation[c]  #outsourced
          PM = ones(1:length(DESTC))-PM1
          #Define Cost Vector
          #############################################################################
          LowerIndenpright=0 
          for i in 1:length(route[r])
            cost[r,c] += CM[route[r][i]]*PM1[route[r][i]]
            cost[r,c] += REWV*A[1,DESTC[route[r][i]]]*PM[route[r][i]]*prod(PM1[route[r][1:i-1]])
            cost[r,c] += REWV*A[DESTC[route[r][i]],1]*PM[route[r][i]]*prod(PM1[route[r][i+1:end]])
            for j in i+1:length(route[r])
              cost[r,c] += REWV*A[DESTC[route[r][i]],DESTC[route[r][j]]]*PM[route[r][i]]*PM[route[r][j]]*prod(PM1[route[r][i+1:j-1]])
              LowerIndenpright += REWV*(A[DESTC[route[r][i]],1]+A[1,DESTC[route[r][j]]]-A[DESTC[route[r][i]],DESTC[route[r][j]]]) *PM[route[r][i]]*PM[route[r][j]]*prod(PM1[route[r][i+1:j-1]])*detour(i,route[r],PM1,PM)
            end
          end
          #############################################################################
          LowerIndenpInv=0
          routeinv=reverse(route[r])
          for i in 1:length(routeinv)
            for j in i+1:length(routeinv)
              LowerIndenpInv += REWV*(A[DESTC[routeinv[i]],1]+A[1,DESTC[routeinv[j]]]-A[DESTC[routeinv[i]],DESTC[routeinv[j]]]) *PM[routeinv[i]]*PM[routeinv[j]]*prod(PM1[routeinv[i+1:j-1]])*detour(i,routeinv,PM1,PM)
            end
          end
          cost[r,c] += min(LowerIndenpright,LowerIndenpInv)
          #print("r= ",r, " c= ", c," , ", cost[r,c], " ;;")
        end
        if  routereducedcosts >  cost[r,c]
          routereducedcosts =  cost[r,c]
          best_c = c
        end
      end
      routereducedcosts = routereducedcosts -getobjectivevalue(model) #+const2_dual #Why error on DUAL!!!!!
      routereducedcosts < -0.001 && break #We look for the first negative reduced cost in order to reduce time spent in pricing problem
    end
    #println("finding routes reducedcost minimal is for ", best_r," comp is ", best_c," result is ", routereducedcosts, " ; ", cost[best_r,best_c])
    #read(STDIN,Char)
    global TimePric += toq()
    tic()
    if routereducedcosts >= 0
      global NOROUTCOUNT += 1 
    end  
    #global TimeSce += toq()
    #tic()
    if routereducedcosts < -0.0001  #Insert new columns
      stop = false
      if !in(best_r,routesused)
        global ROUTCOUNT += 1
        push!(routesused,best_r)  #update solutions used
        println("route number ",ROUTCOUNT, " added ",best_r)

        #calculate new costs
        #Insert new variables Y and Z
        push!(y,@variable(model,  lowerbound = 0, upperbound =1, objective = 0.0, inconstraints = [const2], coefficients = [1.0]))
        #Insert new constraints
        const3new =  @constraint(model,  - y[length(routesused)] == 0)
        const3 = [const3;const3new]
        znew = Array{Variable}(1,length(compensationused))
        for c in 1:length(compensationused)
          if cost[best_r,compensationused[c]] == 0
            #Define compensation and probability Vector 
            #CM = REWOD + REWOD2 .* compensation[compensationused[c]]
            #CM = PRICEOD + 10*PRICEOD .* compensation[compensationused[c]]
            #PM1 = PROBOD + PROBOD2  .* compensation[compensationused[c]]  #outsourced
            CM = 0 + PRICEOD .* compensation[compensationused[c]]
            PM1 = 0 + PROBOD  .* compensation[compensationused[c]]  #outsourced
            PM = ones(1:length(DESTC)).-PM1    #not outsourced
            #Define Cost Vector
            #############################################################################
            LowerIndenpright=0 
            for i in 1:length(route[best_r])
              cost[best_r,compensationused[c]] += CM[route[best_r][i]]*PM1[route[best_r][i]]
              cost[best_r,compensationused[c]] += REWV*A[1,DESTC[route[best_r][i]]]*PM[route[best_r][i]]*prod(PM1[route[best_r][1:i-1]])
              cost[best_r,compensationused[c]] += REWV*A[DESTC[route[best_r][i]],1]*PM[route[best_r][i]]*prod(PM1[route[best_r][i+1:end]])
              for j in i+1:length(route[r])
                cost[best_r,compensationused[c]] += REWV*A[DESTC[route[best_r][i]],DESTC[route[best_r][j]]]*PM[route[best_r][i]]*PM[route[best_r][j]]*prod(PM1[route[best_r][i+1:j-1]])
                LowerIndenpright += REWV*(A[DESTC[route[best_r][i]],1]+A[1,DESTC[route[best_r][j]]]-A[DESTC[route[best_r][i]],DESTC[route[best_r][j]]]) *PM[route[best_r][i]]*PM[route[best_r][j]]*prod(PM1[route[best_r][i+1:j-1]])*detour(i,route[best_r],PM1,PM)
              end
            end
            #############################################################################
            LowerIndenpInv=0
            routeinv=reverse(route[best_r])
            for i in 1:length(routeinv)
              for j in i+1:length(routeinv)
                LowerIndenpInv += REWV*(A[DESTC[routeinv[i]],1]+A[1,DESTC[routeinv[j]]]-A[DESTC[routeinv[i]],DESTC[routeinv[j]]]) *PM[routeinv[i]]*PM[routeinv[j]]*prod(PM1[routeinv[i+1:j-1]])*detour(i,routeinv,PM1,PM)
              end
            end
            cost[best_r,compensationused[c]] += min(LowerIndenpright,LowerIndenpInv)
          end #end cost[r,c] == 0 
          #println("inserir custo r= ",best_r," custo = ",compensationused[c], " : ",cost[best_r,compensationused[c]]) 
          znew[1,c]=@variable(model,  lowerbound = 0, upperbound =1, objective = cost[best_r,compensationused[c]], inconstraints = [const3[length(routesused)]], coefficients = [1.0])
        end 
        z = [z;znew]
      end
      if  !in(best_c,compensationused)
        global COMPCOUNT += 1
        push!(compensationused,best_c)  #update solutions used
        println("comp number ",COMPCOUNT, " added ",best_c)
        #Insert new variables  Z
        znew = Vector{Variable}(length(routesused))
        for r in 1:length(routesused)
          if cost[routesused[r],best_c] == 0
            #Define compensation and probability Vector 
            #CM = REWOD + REWOD2 .* compensation[best_c]
            #CM = PRICEOD + 10*PRICEOD .* compensation[best_c]
            #PM1 = PROBOD + PROBOD2  .* compensation[best_c]  #outsourced
            CM = 0 + PRICEOD .* compensation[best_c]
            PM1 = 0 + PROBOD  .* compensation[best_c]  #outsourced
            PM = ones(1:length(DESTC)).-PM1    #not outsourced
            #Define Cost Vector
            #############################################################################
            LowerIndenpright=0 
            for i in 1:length(route[routesused[r]])
              cost[routesused[r],best_c] += CM[route[routesused[r]][i]]*PM1[route[routesused[r]][i]]
              cost[routesused[r],best_c] += REWV*A[1,DESTC[route[routesused[r]][i]]]*PM[route[routesused[r]][i]]*prod(PM1[route[routesused[r]][1:i-1]])
              cost[routesused[r],best_c] += REWV*A[DESTC[route[routesused[r]][i]],1]*PM[route[routesused[r]][i]]*prod(PM1[route[routesused[r]][i+1:end]])
              for j in i+1:length(route[routesused[r]])
                cost[routesused[r],best_c] += REWV*A[DESTC[route[routesused[r]][i]],DESTC[route[routesused[r]][j]]]*PM[route[routesused[r]][i]]*PM[route[routesused[r]][j]]*prod(PM1[route[routesused[r]][i+1:j-1]])
                LowerIndenpright += REWV*(A[DESTC[route[routesused[r]][i]],1]+A[1,DESTC[route[routesused[r]][j]]]-A[DESTC[route[routesused[r]][i]],DESTC[route[routesused[r]][j]]]) *PM[route[routesused[r]][i]]*PM[route[routesused[r]][j]]*prod(PM1[route[routesused[r]][i+1:j-1]])*detour(i,route[routesused[r]],PM1,PM)
              end
            end
            #############################################################################
            LowerIndenpInv=0
            routeinv=reverse(route[routesused[r]])
            for i in 1:length(routeinv)
              for j in i+1:length(routeinv)
                LowerIndenpInv += REWV*(A[DESTC[routeinv[i]],1]+A[1,DESTC[routeinv[j]]]-A[DESTC[routeinv[i]],DESTC[routeinv[j]]]) *PM[routeinv[i]]*PM[routeinv[j]]*prod(PM1[routeinv[i+1:j-1]])*detour(i,routeinv,PM1,PM)
              end
            end
            cost[routesused[r],best_c] += min(LowerIndenpright,LowerIndenpInv)
          end #end cost[r,c] == 0   
          #println("inserir custo r= ",routesused[r]," custo = ",best_c, " : ",cost[routesused[r],best_c])
          znew[r]=@variable(model,  lowerbound = 0, upperbound =1, objective = cost[routesused[r],best_c], inconstraints = [const3[r]], coefficients = [1.0])
        end 
        z = hcat(z,znew) 
      end
    else
      println("solved = ", getobjectivevalue(model))
      stop = true 
    end
    stop && break 
  end
  zresult = getvalue(z)
  for r in 1:length(routesused)
    for c in 1:length(compensationused)
       if zresult[r,c] != 0
         println("SELECT = ", compensation[compensationused[c]] )
       end
    end
  end
return getobjectivevalue(model),getvalue(y) 

end

#Function Solve Deterministic VRP
         
function solvedetVRP(REGV,V,A,DESTCPAR,REWV,CAPV,DESTOD2,REWOD2,CAPOD2,PR)
  V2 = union(DESTCPAR,DESTOD2)
  V3 = union([1],V2)  
  model = Model(solver = CplexSolver(CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=7200))
  
  @variable(model, x[k=1:REGV+length(DESTOD2),i in V3,j in V3; i != j], Bin)
  @variable(model, y[k=1:REGV+length(DESTOD2),i in DESTCPAR] >= 0 )

  @objective(model, Min, sum(sum(sum( REWV*A[i,j]*x[k,i,j] for k in 1:REGV) for i  in V3 if i!=j) for j  in V3 ) ) 

  @constraint(model, flow1[j in DESTCPAR], sum(sum(x[k,i,j] for k in 1:REGV+length(DESTOD2)) for i in union(DESTCPAR,[1]) if i != j) == 1)
  @constraint(model, flow2[k in 1:REGV+length(DESTOD2)], sum(x[k,1,i] for i in DESTCPAR ) <= 1)
  @constraint(model, flow3[k= 1:REGV,i in union(DESTCPAR,[1]) ], sum(x[k,i,j] for j in union(DESTCPAR,[1]) if i !=j ) - sum(x[k,j,i] for j in union(DESTCPAR,[1]) if i !=j) == 0 )
  #@constraint(model, flow4[k= 1:REGV,i in DESTOD2 ], sum(x[k,i,j] for j in union(DESTC,[1]) if i !=j ) == 0 )
  #@constraint(model, flow5[k= 1+REGV:REGV+length(DESTOD2)], sum(x[k,1,i] for i in DESTC ) ==  sum(x[k,i,DESTOD2[k-REGV]] for i in DESTC if i != DESTOD2[k-REGV] ))
  #@constraint(model, flow6[k= 1+REGV:REGV+length(DESTOD2),i in DESTC ], sum(x[k,i,j] for j in union(DESTC,DESTOD2) if i !=j) - sum(x[k,j,i] for j in union(DESTC,[1]) if i !=j ) == 0 )


  @constraint(model, cap1[k in 1:REGV,i in DESTCPAR, j in DESTCPAR; i != j], y[k,i] - y[k,j] + CAPV*x[k,i,j] <= CAPV -1)
  #@constraint(model, cap2[k in 1+REGV:REGV+length(DESTOD2),i in DESTC, j in DESTC; i != j], y[k,i] - y[k,j] + CAPOD2[k-REGV]*x[k,i,j] <= CAPOD2[k-REGV] -1)
  @constraint(model, cap3[k in 1:REGV,i in DESTCPAR], y[k,i] <= CAPV)
  #@constraint(model, cap4[k in 1+REGV:REGV+length(DESTOD2),i in DESTC], y[k,i] <= CAPOD2[k-REGV])


  solve(model)
  result = getobjectivevalue(model)
  resultx = getvalue(x)
  resulty = getvalue(y)
  #println("Resolvi DET para ", DESTCPAR)
  #for k in  1:REGV
  #  println("Para o veiculo ", k)
  #  sol = []
  #  for i in V3
  #    for j in V3
  #       if i != j
  #         if resultx[k,i,j] != 0
  #            push!(sol,i)
  #            push!(sol,j)
  #         end
  #       end
  #     end
  #   end
  #   println(sol)
  # end
  # read(STDIN,Char)
  global NROUTES = 0
  for k in 1:REGV+length(DESTOD2) 
     NROUTES += sum(resultx[k,1,i] for i in V3 if i != 1)
  end
  #global PERHIGH = sum(globalincumbentsolutiond)/length(globalincumbentsolutiond)
  return result,resultx
end

#Function Solve GODWSKA Algorithm
function solveCATVRP(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,CAPOD)
  V2 = union(DESTC,DESTOD)
  V3 = union(V2,[1])
  DESTC1=union(DESTC,[1])
  #################################
  #Calculate minimum price to pay for ODs
  global PRICEOD = Inf*ones(1:length(DESTC))
  for i in DESTC
    for j in DESTC1
      for r in DESTC1
        if (j !=r && j != i && r != i) || (j==1 && r==1)
          if PRICEOD[findfirst(x -> x==i, DESTC)] > REWV*(A[j,i] + A[i,r] - A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] = REWV*(A[j,i] + A[i,r]- A[j,r])
             PRICEOD[findfirst(x -> x==i, DESTC)] +=0.1
          end
        end
      end
    end
  end
  PRICEOD[1] =0
  PRICEOD[2] =10
  PRICEOD[3] =0
  
  function solveavgCAT(REGV,V,A,DESTC,SELECT,PROBOD,REWOD,REWV,CAPV)
    #generate all possible scenarios for SELECT
    #Create vectors with all possible scenarios...only for toy examples
    #println("entrei de AVGCAT sendo A igual a ", SELECT) 
    scenarioforSELECT = Vector{Vector{Int}}() 
    fillscenario(scenarioforSELECT,zeros(length(SELECT)),length(SELECT),1)
    solution = 0
    multsol = 0
    println("tamanho de scenario for select ", length(scenarioforSELECT)) 
    for w in 1:length(scenarioforSELECT)
      #Define who is U
      U =[]
      mult=1
      for i in 1:length(scenarioforSELECT[w])
        if scenarioforSELECT[w][i] != 0
          push!(U,SELECT[i])
          mult = mult * PROBOD[findfirst(x -> x==SELECT[i], DESTC)]
        else
         mult = mult * (1-PROBOD[findfirst(x -> x==SELECT[i], DESTC)]) 
        end
      end
      println("calculando solpartial para cada scenario ",scenarioforSELECT[w]," U neste = ", U) 
      #calculate VRP without U
      solpartial,solpartialx = solvedetVRP(REGV,V,A,setdiff(DESTC,U),REWV,CAPV,[],[],[],0)
      multsol += mult
      UL =[findfirst(x -> x==U[i], DESTC) for i=1:length(U)]
      println("Select= ",SELECT," cenario= ",scenarioforSELECT[w]," prob= ",mult, "solpartial= ",solpartial,",",  sum(PRICEOD[UL[1:end]]), "probcum= ",multsol)
      #println("UL foi ", UL, " e preco foi ", sum(PRICEOD[UL[1:end]]))
      #read(STDIN,Char)
      solpartial = mult*(solpartial +  sum(PRICEOD[UL[1:end]]) ) #length(U)*REWOD 
      solution += solpartial
    end
    #println("cost of ODs", solution)
    #println("sai de AVGCAT", solution)
    #read(STDIN,Char)
    return solution
  end
  
  #Calculate Deterministic no OD
  zot,zotx = solvedetVRP(REGV,V,A,DESTC,REWV,CAPV,[],[],[],0)
  iot = 0
  println("comecei calculando ZOT supondo A vazio, zot = ", zot)
  SELECT = []
  while true
    #read(STDIN,Char)
    B = copy(SELECT)
    #println("B=",B)
    for i in setdiff(DESTC,B)
      zpar= solveavgCAT(REGV,V,A,DESTC,union(SELECT,[i]),PROBOD,REWOD,REWV,CAPV)
      println("Calculei para A igual a ", union(SELECT,[i]), " deu ",zpar)
      if zpar < zot
        iot=i
        zot = zpar
        SELECT = union(SELECT,[i])
        #println("Resultado melhor SELECT agora,", SELECT)
      end
    end
     #println("B=",B,", ",setdiff(SELECT,B))
    if setdiff(SELECT,B) == []
      #println("Vou break")
      break
    else
      #println("Vou ver se alguem sai, iot =", iot)
      for i in setdiff(SELECT,[iot])
        zpar= solveavgCAT(REGV,V,A,DESTC,setdiff(SELECT,[i]),PROBOD,REWOD,REWV,CAPV)
        ##println("para i = ",i," zpar = ",zpar, "e zot=", zot) 
        if zpar < zot
          zot = zpar
          SELECT = setdiff(SELECT,[i])
          #println("Resultado diminuir SELECT agora,", SELECT)
        end
      end
    end
  end
  
   
return zot,SELECT

end
