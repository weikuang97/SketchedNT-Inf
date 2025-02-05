using JLD2
using Statistics

workdir = "/.../Regression/Table/SketchedNT/Logistic"

include(string(workdir, "/Parameter/Param.jl"))
include(string(workdir, "/StoNewton/StoNewtonMain.jl"))

# load parameters
R_Tail=StoNewtonSet.MaxIter; Rep=StoNewtonSet.Rep; tau=StoNewtonSet.tau;
c_1=StoNewtonSet.c_1; c_2=StoNewtonSet.c_2; c_3=StoNewtonSet.c_3;
SigToe=StoNewtonSet.SigToe; SigEqui=StoNewtonSet.SigEqui;
D=StoNewtonSet.D;
LenTau=length(tau); LenDim=length(D);
LenSigToe=length(SigToe);LenSigEqui = length(SigEqui)

# initialize
Err_PI = zeros(Rep) # error in variance estimation of parameter mean (plug-in)
Err_SC = zeros(Rep) # error in variance estimation of parameter mean (weighted sample cov)
RateCov_PI = zeros(Rep) # coverage rate of confidence interval (plug-in)
RateCov_SC = zeros(Rep) # coverage rate of confidence interval (weighted sample cov)



# go over all cases
for IdDim=1:LenDim, IdTau=1:LenTau
    path1 = string(workdir,"/Solution/Dim",IdDim,"Tau",IdTau)

    # Toeplitz Covariance
    for IdSigToe = 1:LenSigToe
        # reset
        Err_PI = zeros(Rep)
        Err_SC = zeros(Rep)
        RateCov_PI = zeros(Rep)
        RateCov_SC = zeros(Rep)


        # load experiment results
        for IdRep = 1:Rep
            println("IDRep:",IdRep)
            path = string(path1,"/Toe",IdSigToe,"/rep",IdRep,".jld2")

            @load path Result
            COV_value_orc = mean(Result.trueCov)
            Err_PI[IdRep] = (Result.COV_value_PI - COV_value_orc)/COV_value_orc
            Err_SC[IdRep] = (Result.COV_value_SC - COV_value_orc)/COV_value_orc
            RateCov_PI[IdRep] = Result.IdCov_PI
            RateCov_SC[IdRep] = Result.IdCov_SC
        end

        # summarize results
        path2 = string(workdir, "/Figures")
        if !isdir(path2)
            mkpath(path2)
        end

        io = open(string(path2, "/Error.txt"), "a");
        write(io, string("Dim-",D[IdDim],"-Tau-",tau[IdTau],"-Toe-",SigToe[IdSigToe],"\n"))
        write(io, string("Err_PI:", mean(Err_PI),'\n'))
        write(io, string("Err_SC:", mean(Err_SC),'\n'))
        close(io)


        io = open(string(path2, "/Rate.txt"), "a");
        write(io, string("Dim-",D[IdDim],"-Tau-",tau[IdTau],"-Toe-",SigToe[IdSigToe],"\n"))
        write(io, string("CovRate_PI:", mean(RateCov_PI),'\n'))
        write(io, string("CovRate_SC:", mean(RateCov_SC),'\n'))
        close(io)
    end


    # Equi-corr Covariance
    for IdSigEqui = 1:LenSigEqui
        # reset
        Err_PI = zeros(Rep)
        Err_SC = zeros(Rep)
        RateCov_PI = zeros(Rep)
        RateCov_SC = zeros(Rep)


        # load experiment results
        for IdRep = 1:Rep
            println("IDRep:",IdRep)
            path = string(path1,"/Equ",IdSigEqui,"/rep",IdRep,".jld2")

            @load path Result
            COV_value_orc = mean(Result.trueCov)
            Err_PI[IdRep] = (Result.COV_value_PI - COV_value_orc)/COV_value_orc
            Err_SC[IdRep] = (Result.COV_value_SC - COV_value_orc)/COV_value_orc
            RateCov_PI[IdRep] = Result.IdCov_PI
            RateCov_SC[IdRep] = Result.IdCov_SC
        end

        # summarize results
        path2 = string(workdir, "/Figures")
        if !isdir(path2)
            mkpath(path2)
        end

        io = open(string(path2, "/Error.txt"), "a");
        write(io, string("Dim-",D[IdDim],"-Tau-",tau[IdTau],"-Equ-",SigEqui[IdSigEqui],"\n"))
        write(io, string("Err_PI:", mean(Err_PI),'\n'))
        write(io, string("Err_SC:", mean(Err_SC),'\n'))
        close(io)


        io = open(string(path2, "/Rate.txt"), "a");
        write(io, string("Dim-",D[IdDim],"-Tau-",tau[IdTau],"-Equ-",SigEqui[IdSigEqui],"\n"))
        write(io, string("CovRate_PI:", mean(RateCov_PI),'\n'))
        write(io, string("CovRate_SC:", mean(RateCov_SC),'\n'))
        close(io)
    end
end
