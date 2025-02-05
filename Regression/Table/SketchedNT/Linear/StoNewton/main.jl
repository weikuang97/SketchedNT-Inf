using LinearAlgebra
using Random
using Distributions
using JLD2

cd("/.../Regression/Table/SketchedNT/Linear/StoNewton")



#using Main.Parameter
include("../Parameter/Param.jl")

include("StoNewtonMain.jl")
#######################################
#########  run main file    ###########
#######################################
function main()
    StoNewtonMain(StoNewtonSet)
end

main()
