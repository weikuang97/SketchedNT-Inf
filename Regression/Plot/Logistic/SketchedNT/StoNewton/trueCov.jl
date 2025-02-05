# calculate the true covariance matrix
function trueCov(Omegastar, Bstar, d, tau)
    if (tau==0)
        Xistar = Omegastar/2
    else
        Bstar_norm = mapslices(normalize, Bstar, dims=1)
        Dstar_j = (Bstar_norm*Bstar_norm')./d
        Cstar_j = Matrix(1.0I,d,d)-Dstar_j
        Cstar = (Cstar_j^tau).*(-1)

        eigen_result = eigen(Hermitian(Matrix(1.0I,d,d)+Cstar))
        sigma = eigen_result.values
        U = eigen_result.vectors
        Theta = inv.(sigma*ones(d)'+ones(d)*sigma')

        mid_temp = Omegastar
        for i = 1:tau
            mid_temp = onestep(mid_temp, Bstar_norm, Dstar_j, d)
        end
        exp_term = Omegastar + Cstar*Omegastar + Omegastar*Cstar' + mid_temp

        Xistar = U*(Theta.*(U'*exp_term*U))*U'
    end

    return Xistar
end


function onestep(Omegastar, Bstar_norm, Dstar_j, d)
    wt = mapslices(vec -> weight(vec, Omegastar), Bstar_norm, dims=1)[1,:]
    term = (Bstar_norm*Matrix(Diagonal(wt))*Bstar_norm')./d
    res = Omegastar - Dstar_j*Omegastar - Omegastar*Dstar_j' + term
    return res
end

function weight(vec, Omegastar)
    return vec'*Omegastar*vec
end

function Bstar_Omegastar(Sigma, X_true, d, mu)
    S = zeros(d,d)
    for t=1:Int(3e5)
        a_t = rand(MvNormal(zeros(d),Sigma))
        bar_nab_x2f_t = a_t*a_t'/((1+exp(a_t'*X_true))*(1+exp(-a_t'*X_true)))
        S = (t-1)/t*S +1/t*bar_nab_x2f_t
    end
    pert = mu*Matrix(1.0I,d,d)
    Bstar = Matrix(S+pert)
    Omegastar = inv(Bstar)*S*inv(Bstar)
    return Bstar, Omegastar
end
