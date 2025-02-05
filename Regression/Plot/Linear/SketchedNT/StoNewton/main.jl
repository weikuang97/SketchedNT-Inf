using LinearAlgebra
using Random
using Distributions
using DataFrames
using CSV


cd("/.../Regression/Plot/Linear/SketchedNT/StoNewton")



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
