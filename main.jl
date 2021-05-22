#########################################################################
# Module to  run DDU STOCHASTIC VRP OD

using JuMP, CPLEX, MathProgBase, DataStructures #, Distributions
#using LightGraphs,Cairo, Fontconfig,Compose,GraphPlot, Colors,Plots

original_sum = sum

include("../src/data.jl")
include("../src/OptBBfunctions.jl")
include("../src/Utilfunctions.jl")


#Print results to file
appfolder = dirname(@__FILE__) * "/../"   #get current folder
filname=appfolder* "/result/"*"RESULT"*string(now())*".txt"  #define filename to be written
resultado = open(filname, "a")                            #opens file

#Read instance
#XENIA:  INSERT HERE INSTANCES TO RUN. INSTANCES ARE HARD CODED IN data.jl
instances =["12CSA","12CSB","12CSC"] #"6CSA","6CSB","6CSC","7CSA","7CSB","7CSC","8CSA","8CSB","8CSC",]#"6A","6C","7B","IT1","IT2"] #"7B"]#,"6A","6C","7B","8A"]#,"6B","6C","6D","7A","7B","7C","7D"]#,"8A","8B","8C","8D"]#,"9A","9B","9C","9D"]#,"10A","10B","10C","10D","11A","11B","11C","11D","12A","12B","12C","12D"] 
instance = []

for inst_id in instances
V,A,DESTC,PROBC,REGV,REWV,CAPV,DESTOD,PROBOD2,REWOD2,PROBOD,REWOD,CAPOD = read_data(instance, inst_id)

global inst_idf= inst_id

##################################################################################################################################


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

#result,initialsol= solvedetVRP(REGV,V,A,DESTC,REWV,CAPV,DESTOD2,REWOD2,CAPOD2,PR)
#result,initialsol= solvedetVRP(REGV,V,A,DESTC,REWV,CAPV,[],[],0,1)
TimeMain += toq()
#println(resultado,inst_id,";","DETVRP= ",";",result,";",TimeMain+TimePric,";",TimeMain,";",TimePric,";",";",TimeSce,";",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#println("DETVRP IS ", result, " under a time of ", TimeMain, " and number of routes = ", NROUTES)
#read(STDIN,Char)
#####################################################################################################################################

##################################################################################################################################


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
#result,initialsol= solveSTODDUEXACTN1(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
TimeMain += toq()
#println(resultado,inst_id,";","DDUEXACTN1= ",";",result,";",TimeMain+TimePric,";",TimeMain,";",TimePric,";",";",TimeSce,";",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#println("Final Solution DDU Stochastic EXACTN1= ",inst_id," ,", result,"  Time spent = ",TimeMain+TimePric+TimeSce,"  Time Main = ",TimeMain,"  Time Pric = ",TimePric,"  Time Comp = ",TimeSce," Nodes =",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#read(STDIN,Char)
#####################################################################################################################################




##################################################################################################################################


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

#result,initialsol= solveCATVRP(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,CAPOD)
#TimeMain += toq()
#println("GDOWSKA ALG SOLUTION IS ", result," and selected customers to outsource are ", initialsol, " under a time of ", TimeMain)
#read(STDIN,Char)
#####################################################################################################################################

##################################################################################################################################


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
#result,initialsol= solveSTODDUEXACTN10(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
TimeMain += toq()
#println(resultado,inst_id,";","DDUEXACTN10= ",";",result,";",TimeMain+TimePric,";",TimeMain,";",TimePric,";",";",TimeSce,";",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#println("Final Solution DDU Stochastic EXACTN10= ",inst_id," ,", result,"  Time spent = ",TimeMain+TimePric+TimeSce,"  Time Main = ",TimeMain,"  Time Pric = ",TimePric,"  Time Comp = ",TimeSce," Nodes =",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#read(STDIN,Char)
#####################################################################################################################################



##################################################################################################################################


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
result,initialsol= solveSTODDUEXACTN11(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
TimeMain += toq()
println(resultado,inst_id,";","DDUEXACTN11= ",";",result,";",TimeMain+TimePric,";",TimeMain,";",TimePric,";",";",TimeSce,";",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
println("Final Solution DDU Stochastic EXACTN11= ",inst_id," ,", result,"  Time spent = ",TimeMain+TimePric,"  Time Main = ",TimeMain,"  Time Pric = ",TimePric,"  Time Comp = ",TimeSce," Nodes =",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#read(STDIN,Char)
#####################################################################################################################################

##################################################################################################################################


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
#result,initialsol= solveSTODDUEXACTN3(REGV,V,A,DESTC,PROBC,REWV,CAPV,DESTOD,PROBOD,REWOD,PROBOD2,REWOD2,CAPOD)
TimeMain += toq()
#println(resultado,inst_id,";","DDUEXACTN3= ",";",result,";",TimeMain+TimePric,";",TimeMain,";",TimePric,";",";",TimeSce,";",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#println("Final Solution DDU Stochastic EXACTN3= ",inst_id," ,", result,"  Time spent = ",TimeMain+TimePric,"  Time Main = ",TimeMain,"  Time Pric = ",TimePric,"  Time Comp = ",TimeSce," Nodes =",NODECOUNT,"; ",ROUTCOUNT,"; ",COMPCOUNT,"; ",NOROUTCOUNT,"; ",NOCOMPCOUNT,"; ",NROUTES,"; ",PERHIGH,";",ROOTSOL,";",ROOTTIME)
#read(STDIN,Char)
#####################################################################################################################################



end
close(resultado)




