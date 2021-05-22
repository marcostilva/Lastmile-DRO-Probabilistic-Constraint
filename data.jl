

function read_data(instance, inst_id)

#Define size of grid to define Graph
grid= 10  #will be a 10 x 10 euclidean grid to include depot, customers and ODs

#Define graph G(V,A)

V = collect(1:grid * grid )  #Depot is vertice/node 1 locate in grid position (1,1)
A = 1.0*zeros(grid * grid, grid * grid)
for i in V                    #i is a vertice/node within the grid
  for j in V    #j is a vertice/node within the grid  
     A[i,j] = sqrt( (div(i-1,grid) - div(j-1,grid))^2 + (mod(i-1,grid) - mod(j-1,grid))^2)  #building complete graph
  end
end


DESTC = []
PROBC= []
REGV = 0
REWV = 0
CAPV = 0
DESTOD= []
PROBOD2= []
REWOD2= []
PROBOD= []
REWOD= []
CAPOD= []
if  inst_id == "6CSA"  #Only one level of Prob and Price
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,44,55,88,95]
  PROBC = [0.4,0.35,0.3,0.25,0.2,0.1]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0,0,0,0,0,0]
  REWOD2 = [0,0,0,0,0,0]
  PROBOD = [0.4,0.35,0.3,0.25,0.2,0.1]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 3
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif  inst_id == "6CSB"  #First level no outsource, second level minimim PRICE
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,44,55,88,95]
  PROBC = [0,0,0,0,0,0]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.35,0.3,0.25,0.2,0.1]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0,0,0,0,0,0]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0,0,0,0,0,0]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 3
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif  inst_id == "6CSC"  #First level and second level equal
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,44,55,88,95]
  PROBC = [0.2,0.2,0.2,0.2,0.2,0.2]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.4,0.4,0.4,0.4,0.4]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0.2,0.2,0.2,0.2,0.2,0.2]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 3
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif  inst_id == "7CSA"  #Only one level of Prob and Price
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,38,44,55,88,95]
  PROBC = [0.4,0.35,0.32,0.3,0.25,0.2,0.1]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0,0,0,0,0,0,0]
  REWOD2 = [0,0,0,0,0,0,0]
  PROBOD = [0.4,0.35,0.32,0.3,0.25,0.2,0.1]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif  inst_id == "7CSB"  #First level no outsource, second level minimim PRICE
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,38,44,55,88,95]
  PROBC = [0,0,0,0,0,0,0]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.35,0.32,0.3,0.25,0.2,0.1]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0,0,0,0,0,0,0]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0,0,0,0,0,0,0]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif  inst_id == "7CSC"  #First level and second level equal
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,38,44,55,88,95]
  PROBC = [0.2,0.2,0.2,0.2,0.2,0.2,0.2]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.4,0.4,0.4,0.4,0.4,0.4]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0.2,0.2,0.2,0.2,0.2,0.2,0.2]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif  inst_id == "8CSA"  #Only one level of Prob and Price
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,38,44,55,63,88,95]
  PROBC = [0.4,0.35,0.32,0.3,0.25,0.22,0.2,0.1]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0,0,0,0,0,0,0,0]
  REWOD2 = [0,0,0,0,0,0,0,0]
  PROBOD = [0.4,0.35,0.32,0.3,0.25,0.22,0.2,0.1]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif  inst_id == "8CSB"  #First level no outsource, second level minimim PRICE
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,38,44,55,63,88,95]
  PROBC = [0,0,0,0,0,0,0]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.35,0.32,0.3,0.25,0.22,0.2,0.1]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0,0,0,0,0,0,0,0]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0,0,0,0,0,0,0,0]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif  inst_id == "8CSC"  #First level and second level equal
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,38,44,55,63,88,95]
  PROBC = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif  inst_id == "10CSA"  #Only one level of Prob and Price
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,35,36,38,44,55,63,88,95]
  PROBC = [0.4,0.35,0.32,0.3,0.3,0.3,0.25,0.22,0.2,0.1]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0,0,0,0,0,0,0,0,0,0]
  REWOD2 = [0,0,0,0,0,0,0,0,0,0]
  PROBOD = [0.4,0.35,0.32,0.3,0.3,0.3,0.25,0.22,0.2,0.1]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 5
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif  inst_id == "10CSB"  #First level no outsource, second level minimim PRICE
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,35,36,38,44,55,63,88,95]
  PROBC = [0,0,0,0,0,0,0,0,0]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.35,0.32,0.3,0.3,0.3,0.25,0.22,0.2,0.1]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0,0,0,0,0,0,0,0,0,0]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0,0,0,0,0,0,0,0,0,0]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 5
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif  inst_id == "10CSC"  #First level and second level equal
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,35,36,38,44,55,63,88,95]
  PROBC = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 5
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif  inst_id == "12CSA"  #Only one level of Prob and Price
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,23,24,35,36,38,44,55,63,88,95]
  PROBC = [0.4,0.35,0.32,0.32,0.32,0.3,0.3,0.3,0.25,0.22,0.2,0.1]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0,0,0,0,0,0,0,0,0,0,0,0]
  REWOD2 = [0,0,0,0,0,0,0,0,0,0,0,0]
  PROBOD = [0.4,0.35,0.32,0.32,0.32,0.3,0.3,0.3,0.25,0.22,0.2,0.1]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 6
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif  inst_id == "12CSB"  #First level no outsource, second level minimim PRICE
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,23,24,35,36,38,44,55,63,88,95]
  PROBC = [0,0,0,0,0,0,0,0,0,0,0]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.35,0.32,0.32,0.32,0.3,0.3,0.3,0.25,0.22,0.2,0.1]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0,0,0,0,0,0,0,0,0,0,0,0]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0,0,0,0,0,0,0,0,0,0,0,0]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 6
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif  inst_id == "12CSC"  #First level and second level equal
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,23,24,35,36,38,44,55,63,88,95]
  PROBC = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4]
  REWOD2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
  PROBOD = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 6
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif  inst_id == "6A"
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,44,55,88,95]
  PROBC = [0.1,0.2,0.3,0.5,0.8,0.95]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.9,0.8,0.5,0.3,0.2,0.05]
  REWOD2 = [0.9,0.8,0.5,0.3,0.2,0.1]
  PROBOD = [0.1,0.2,0.3,0.5,0.8,0.95]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.2,0.3,0.5,0.8,0.95]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 3
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif  inst_id == "IT1"
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [3,5,7,47,58,69]
  PROBC = [0.3,0.3,0.3,0.3,0.3,0.3]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.3,0.3,0,0.3,0.3,0.3]
  REWOD2 = [0.9,0.8,0.5,0.3,0.2,0.1]
  PROBOD = [0.3,0.3,0,0.3,0.3,0.3]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.2,0.3,0.5,0.8,0.95]  #REWARD WILL BE RECALCULATED WITHIN PROGRAM
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 3
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit
elseif inst_id == "IT2"
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,44,55,58,66,88,95,77,79]
  PROBC = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]
  REWOD2 = [0.05,0.2,0.3,0.5,0.5,0.7,0.8,0.95,0.3,0.3]
  PROBOD = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.2,0.3,0.5,0.5,0.7,0.8,0.95,0.3,0.3]
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 5
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif inst_id == "4A"
 #Define customers Destination nodes, marginal probabilities, capacities 
  #DESTC = [11,21,44,55,95,88]
  DESTC = [2,3,7,77]
  PROBC = [0.3,0.3,0.3,0.3]
  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.3,0.3,0.0,0.3]
  REWOD2 = [0.3,0.3,0.3,0.3]
  PROBOD = [0.3,0.3,0.0,0.3]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.3,0.3,0.3,0.3]
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 2
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit


elseif inst_id == "6C"
 #Define customers Destination nodes, marginal probabilities, capacities 
  #DESTC = [11,21,44,55,95,88]
  DESTC = [2,3,7,77,95,99]
  PROBC = [0.3,0.3,0.3,0.3,0.3,0.3]
  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.3,0.3,0.0,0.3,0.3,0.3]
  REWOD2 = [0.3,0.3,0.3,0.3,0.3,0.3]
  PROBOD = [0.3,0.3,0.0,0.3,0.3,0.3]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.3,0.3,0.3,0.3,0.3,0.3]
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 3
  REWV = 2
  CAPV = 6   #demand is given in demands unit, each customer is fixed at one unit



elseif inst_id == "8A"
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,44,55,58,66,88,95]
  PROBC = [0.1,0.2,0.3,0.5,0.5,0.7,0.8,0.95]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.1,0.2,0.3,0.05,0.05,0.05,0.05,0.05]
  REWOD2 = [0.05,0.2,0.3,0.5,0.5,0.7,0.8,0.95]
  PROBOD = [0.1,0.2,0.3,0.5,0.5,0.7,0.8,0.95]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.2,0.3,0.5,0.5,0.7,0.8,0.95]
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 4   #demand is given in demands unit, each customer is fixed at one unit

elseif inst_id == "7B"
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [22,44,55,58,66,88,95]
  PROBC = [0,0.7,0,0.3,0.3,0.3,0.3]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.2,0.3,0.05,0.05,0.05,0.05,0.05]
  REWOD2 = [0.2,0.3,0.5,0.5,0.7,0.8,0.95]
  PROBOD = [0,0.7,0,0.3,0.3,0.3,0.3]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.2,0.3,0.5,0.5,0.7,0.8,0.95]
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 2   #demand is given in demands unit, each customer is fixed at one unit

elseif inst_id == "9A"
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,44,55,58,66,71,95,88]
  PROBC = [0.1,0.2,0.3,0.5,0.5,0.6,0.7,0.95,0.8]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.1,0.2,0.3,0.2,0.2,0.2,0.2,0.05,0.1]
  REWOD2 = [0.1,0.2,0.3,0.5,0.5,0.6,0.7,0.95,0.8]
  PROBOD = [0.1,0.2,0.3,0.5,0.5,0.6,0.7,0.95,0.8]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.2,0.3,0.5,0.5,0.6,0.7,0.95,0.8]
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 3
  REWV = 2
  CAPV = 9
   #demand is given in demands unit, each customer is fixed at one unit

elseif inst_id == "10A"
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,39,44,55,58,66,71,95,88]
  PROBC = [0.1,0.2,0.3,0.3,0.5,0.5,0.6,0.7,0.95,0.8]

  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.3,0.2,0.3,0.2,0.2,0.1,0.1,0.1,0,0]
  REWOD2 = [0.1,0.2,0.3,0.3,0.5,0.5,0.6,0.7,0.95,0.8]
  PROBOD = [0.1,0.2,0.3,0.3,0.5,0.5,0.6,0.7,0.95,0.8]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.1,0.2,0.3,0.3,0.5,0.5,0.6,0.7,0.95,0.8]
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 4
  REWV = 2
  CAPV = 3   #demand is given in demands unit, each customer is fixed at one unit

elseif inst_id == "11C"
 #Define customers Destination nodes, marginal probabilities, capacities 
  DESTC = [11,21,39,44,49,55,58,66,71,95,88]
  PROBC = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]
  #Define OD Destination nodes and marginal probabilities and fixed reward per distance traveled
  DESTOD = [] #ODs nodes must be different then customers!!!
  PROBOD2 = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]
  REWOD2 = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]
  PROBOD = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]    #For decision dependent solutioning. two levels of probability and reward
  REWOD = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]
  CAPOD = []   #capacity is given in demands unit, each customer is fixed at one unit


  #Define number of regular vehicles and fixed reward per distance traveled

  REGV = 6
  REWV = 2
  CAPV = 11 #2   #demand is given in demands unit, each customer is fixed at one unit




else
  println("NO INSTANCES FOUND")

end
return V,A,DESTC,PROBC,REGV,REWV,CAPV,DESTOD,PROBOD2,REWOD2,PROBOD,REWOD,CAPOD
end

