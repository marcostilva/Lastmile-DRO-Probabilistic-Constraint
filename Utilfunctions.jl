#########################################################################
# Module containing utility functions  for VRPOD
#
#

#Sequence of Functions
function fillscenario(where, temp,size,position)
  if position == size+1
    push!(where,temp)
  else
    temp[position]=0
    fillscenario(where, temp,size,position+1)
    temp[position]=1
    fillscenario(where,temp,size,position+1)
  end
end


function troca(vetor, nodei, nodej)
  aux = vetor[nodei]
  vetor[nodei] = vetor[nodej]
  vetor[nodej] = aux;
end

function fillroute(vetor, inf, sup,route)
   
  if inf == sup
    data4 = copy(vetor)
    if data4[1] <= data4[end]
      push!(route,data4)
    end
  else
    for i = inf:sup
      troca(vetor, inf, i)
      fillroute(vetor, inf + 1, sup,route)
      troca(vetor, inf, i) #backtracking
    end
  end
end



