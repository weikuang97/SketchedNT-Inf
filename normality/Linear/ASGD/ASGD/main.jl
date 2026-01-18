using LinearAlgebra
using Random
using Distributions
using DataFrames
using CSV

# cd("/Users/weikuang/Desktop/UChicago/StoOpt/IMA/github/NewtonInf/Regression/Plot/Linear/ASGD/ASGD")



#using Main.Parameter
include("../Parameter/Param.jl")

include("ASGDMain.jl")
#######################################
#########  run main file    ###########
#######################################
function main()
    ASGDMain(ASGDSet)
end

main()
