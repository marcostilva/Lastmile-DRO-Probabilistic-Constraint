######################################################
# Module to  run DRO Crowdshipping VRP with probabilistic constraint using column-dependent row pricing algorithm

using JuMP, CPLEX, MathProgBase, StatsBase, DataStructures 

original_sum = sum

include("OptBBfunctions.jl")
include("Utilfunctions.jl")


#Print results to file
appfolder = dirname(@__FILE__) * "/../"   #get current folder
filname=appfolder* "/result/"*"RESULT"*string(now())*".txt"  #define filename to be written
resultado = open(filname, "a")                            #opens file

#Generate instance


#Ðefine size of grid to define Graph
grid= 10 #will be a 10 x 10 euclidean grid to include depot, customers and ODs

#Define graph G(V,A)
V = collect(1:grid * grid )  #Depot is vertice/node 1 locate in grid position (1,1)
A = 1.0*zeros(grid * grid, grid * grid)
for i in V                    #i is a vertice/node within the grid
  for j in V    #j is a vertice/node within the grid  
     A[i,j] = sqrt( (div(i-1,grid) - div(j-1,grid))^2 + (mod(i-1,grid) - mod(j-1,grid))^2)  #building complete graph
  end
end
#Create Instance Variables
DESTC = []
PROBC= []
REWV = 0
CAPV = 0
REWOD= []
PROBOD= []
  
siz= 7 #Number of Customers
inst_id = siz
#Define customers Destination nodes, marginal probabilities, capacities

DESTC = sample(2:grid*grid, siz, replace = false)
DESTC1=union(DESTC,[1])
#print(DESTC)
PROBC = sample(linspace(.1,.3,3), siz, replace = true)
#print(PROBC)
REWV = 2
CAPV = floor(siz/3) #demand is given in demands unit, each customer is fixed at one unit
REWOD= []
PROBOD= sample(linspace(.1,.3,3), siz, replace = true)


##################################################################################################################################

#Define Control variables to be reported at end
global heursol = 0
global initialsol = []
global NODECOUNT = 0
global ROUTCOUNT = 0
global COMPCOUNT = 0
global NOROUTCOUNT = 0
global NOCOMPCOUNT = 0
global TimeMain = 0
global TimePric = 0
global TimeSce = 0
global NROUTES = 0
global PERHIGH = 0
global ROOTSOL = 0
global ROOTTIME = 0


#Run main routine 
tic()
result,initialsol= solveSTODROPCEXACT(V,A,DESTC,PROBC,REWV,CAPV,PROBOD,REWOD)
TimeMain += toq()
println("Final Solution DRO Stochastic EXACT= ",inst_id," ,", result,"  Time spent = ",TimeMain+TimePric+TimeSce,"  Time Main = ",TimeMain,"  Time Pric = ",TimePric,"  Time Comp = ",TimeSce," Nodes =",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#read(STDIN,Char)
println(resultado,"Final Solution DRO Stochastic EXACT=  ",inst_id," ,", result,"  Time spent = ",TimeMain+TimePric,"  Time Main = ",TimeMain,"  Time Pric = ",TimePric,"  Time Comp = ",TimeSce," Nodes =",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#####################################################################################################################################
readline()
global heursol = 0
global initialsol = []
global NODECOUNT = 0
global ROUTCOUNT = 0
global COMPCOUNT = 0
global NOROUTCOUNT = 0
global NOCOMPCOUNT = 0
global TimeMain = 0
global TimePric = 0
global TimeSce = 0
global NROUTES = 0
global PERHIGH = 0
global ROOTSOL = 0
global ROOTTIME = 0

tic()
result,initialsol= solveSTODROPCHEUR(V,A,DESTC,PROBC,REWV,CAPV,PROBOD,REWOD)
TimeMain += toq()
println("Final SolutionDRO Stochastic HEUR=  ",inst_id," ,", result,"  Time spent = ",TimeMain+TimePric+TimeSce,"  Time Main = ",TimeMain,"  Time Pric = ",TimePric,"  Time Comp = ",TimeSce," Nodes =",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#read(STDIN,Char)
println(resultado,"Final Solution DDU Stochastic HEUR= ",inst_id," ,", result,"  Time spent = ",TimeMain+TimePric,"  Time Main = ",TimeMain,"  Time Pric = ",TimePric,"  Time Comp = ",TimeSce," Nodes =",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#####################################################################################################################################
readline()


close(resultado)



