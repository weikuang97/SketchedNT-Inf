using LinearAlgebra
using Random
using Distributions
using DataFrames
using CSV

# cd("/.../normality/Linear/Toep/ASGD/ASGD")



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
