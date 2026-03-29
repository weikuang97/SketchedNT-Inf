using LinearAlgebra
using Random
using Distributions
using DataFrames
using CSV

# cd("/.../normality/Linear/Toep/SGD/SGD")



#using Main.Parameter
include("../Parameter/Param.jl")

include("SGDMain.jl")
#######################################
#########  run main file    ###########
#######################################
function main()
    SGDMain(SGDSet)
end

main()
