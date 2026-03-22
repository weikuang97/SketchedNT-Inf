Code for Section 5 in paper "Online Covariance Matrix Estimation in Sketched Newton Methods".
Author: Wei Kuang, Mihai Anitescu, and Sen Na

Please update the path and parameters in the scripts before you run the code. 


Sec 5.1 Linear regression
Plot:
generate data: include(".../Regression/Plot/Linear/SketchedNT/StoNewton/main.jl")
summarize data: include(".../Regression/Plot/Linear/SketchedNT/summarize.jl")
generate plot: run ".../Regression/Plot/Linear/plot.R"

Table:
generate data: include(".../Regression/Table/SketchedNT/Linear/StoNewton/main.jl")
generate table: include(".../Regression/Table/SketchedNT/Linear/summarize.jl")


Sec 5.2 Logistic regression
Plot:
generate data: include(".../Regression/Plot/Logistic/SketchedNT/StoNewton/main.jl")
summarize data: include(".../Regression/Plot/Logistic/SketchedNT/summarize.jl")
generate plot: run ".../Regression/Plot/Logistic/plot.R"

Table:
generate data: include(".../Regression/Table/SketchedNT/Logistic/StoNewton/main.jl")
generate table: include(".../Regression/Table/SketchedNT/Logistic/summarize.jl")


Sec 5.3 Inference under different skectching configurations
generate data: run "main.jl" in the corresponding folder for all sketching configurations
generate table: run "summarize.jl" in the corresponding folder for all sketching configurations


Sec 5.4 Hessian preconditioning
generate data: run "main.jl" and "summarize.jl" in ASGD, SGD, and SketchedNT respectively
generate plot: move all summarized data for all algoritms and settings to one folder, set up path, and run "plot.R"


Sec 5.5 impact of sketching iteration number tau
similar to Sec 5.1 and 5.2
generate data: run "main.jl"
summarize data: run "summarize.jl"
generate plot: run "plot.R"


Sec 5.6 CUTEst problems
generate data: include(".../CUTEst/StoSQP/main.jl")
generate table: include(".../CUTEst/summarize.jl")