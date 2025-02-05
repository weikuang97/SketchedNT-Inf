include("trueCov.jl")

function StoSQP(nlp,c1,c2,c3,sigma2,Max_Iter,ttau,X_true,Lam_true,EPS=1e-8)
    # set Buffer size
    Buff= Int(1e5)

    nx, nlam = nlp.meta.nvar, nlp.meta.ncon
    G_star, B_star = jac(nlp,X_true), hess(nlp,X_true,Lam_true)
    K_star = hcat(vcat(B_star,G_star),vcat(G_star',zeros(nlam, nlam)))
    K_star = Matrix(K_star)
    CovM = sigma2*(Matrix(I,nx,nx)+ones(nx,nx)) # sampling covariance
    Xistar = trueCov(K_star, CovM, nx, nlam, ttau) # true limiting covariance

    G_star = Matrix(G_star)
    if rank(G_star*G_star') == size(G_star,1)
        # Index for inference
        proj = G_star'*inv(G_star*G_star')*G_star
        Index = findall(x -> x >= 0.05, norm.(eachcol(proj-Matrix(1.0I,nx,nx)), 2))
    end

    # Initialize variables
    t, X_t, Lam_t, Step_size = 1, nlp.meta.x0, nlp.meta.y0, 0 # iterate

    NewDir_t, K_t = zeros(nx+nlam), zeros(nx+nlam, nx+nlam)
    cum_bar2x_f_t = Matrix(I,nx,nx)

    W_t,V_t,u_t,barxlam_t = zeros(nx+nlam,nx+nlam),zeros(nx+nlam),0,zeros(nx+nlam)
    Xi_t_SC = zeros(nx+nlam,nx+nlam) # weighted sample cov

    ErrXLam = [] # iterate error
    IDCov_SC = [] # 1 or 0, true parameter covered by confidence interval or not
    COV_value_SC = [] # variance estimation
    push!(ErrXLam, norm([X_t-X_true;Lam_t-Lam_true]))

    # Start the loop
    Time = time()
    while t <= Max_Iter
        # Step 0: compute deterministic quantities
        c_t, G_t = consjac(nlp, X_t)
        cum_bar2x_L_t = cum_bar2x_f_t + hess(nlp,X_t,Lam_t,obj_weight=0.0)

        # Constraints value and Jacobian
        # Compute the reduced Hessian and do Hessian modification
        Q_t, R_t = qr(Matrix(G_t'))
        if abs(R_t[end,end])< EPS || ErrXLam[end]>1e5
            println("not converge")
            return 0,zeros(nx,nx),[],[],[],1
        end
        delta_tt = eigmin(Hermitian(Q_t[:,nlam+1:end]'cum_bar2x_L_t*Q_t[:,nlam+1:end],:L))
        if nlam<nx && delta_tt <= EPS
            delta_t = abs(delta_tt) + 0.1
        else
            delta_t = 0
        end
        K_t = hcat(vcat(cum_bar2x_L_t+delta_t*Matrix(I,nx,nx),G_t),vcat(G_t',zeros(nlam,nlam)))

        # Step 1: generate random gradient and Hessian
        g_t,G_tTLam_t,nab_x2f_t = grad(nlp,X_t),jtprod(nlp,X_t,Lam_t),hess(nlp,X_t,zeros(nlam))
        barg_t = rand(MvNormal(g_t,CovM))
        bar_nab_L_t = [barg_t+G_tTLam_t; c_t]
        Rand_t = rand(Normal(0,sigma2^(1/2)),nx,nx)
        bar_nab_x2f_t = Hermitian(nab_x2f_t,:L) + (Rand_t+Rand_t')/2
        cum_bar2x_f_t = t/(t+1)*cum_bar2x_f_t + 1/(t+1)*bar_nab_x2f_t


        # Step 2: solve Newton system via Randomized Solve
        NewDir_t = zeros(nx+nlam)
        if ttau == 0
            NewDir_t = lu(K_t)\-bar_nab_L_t
        else
            for inner_iter = 1:ttau
                # random pick a index
                j = sample(1:(nx+nlam))
                NewDir_t = NewDir_t - (K_t[j,:]'*NewDir_t+bar_nab_L_t[j])/norm(K_t[j,:])^2*K_t[:,j]
            end
        end


        # Step 3: update the iterate
        beta_t = c1/t^c2
        chi_t = beta_t^c3
        phi_t = beta_t+chi_t/2
        alpha_t = rand(Uniform(beta_t,beta_t+chi_t))
        X_t = X_t +alpha_t*NewDir_t[1:nx]
        Lam_t = Lam_t + alpha_t*NewDir_t[nx+1:end]
        push!(ErrXLam, norm([X_t-X_true;Lam_t-Lam_true]))


        # Step 4: covariance estimation
        if t > Buff
            W_t = (t-Buff-1)/(t-Buff)*W_t + 1/(t-Buff)*[X_t;Lam_t]*[X_t;Lam_t]'/phi_t
            V_t = (t-Buff-1)/(t-Buff)*V_t + 1/(t-Buff)*[X_t;Lam_t]/phi_t
            u_t = (t-Buff-1)/(t-Buff)*u_t + 1/(t-Buff)*(1/phi_t)
            barxlam_t = (t-Buff-1)/(t-Buff)*barxlam_t + 1/(t-Buff)*[X_t;Lam_t]
            Xi_t_SC = W_t-V_t*barxlam_t'-barxlam_t*V_t'+u_t*barxlam_t*barxlam_t'


            # confidence interval
            Cov_value_SC = sum(Xi_t_SC[Index, Index])/(length(Index))^2
            if Cov_value_SC<1e-6
                Cov_value_SC = 0
            end
            radius_SC = 1.96*sqrt(alpha_t)*sqrt(Cov_value_SC)

            push!(COV_value_SC, Cov_value_SC)


            Xmean = mean(X_t[Index])
            IdCov_SC = (mean(X_true[Index])>=(Xmean-radius_SC)) && (mean(X_true[Index])<=(Xmean+radius_SC))

            push!(IDCov_SC, IdCov_SC)
        end

        t = t + 1
    end
    Time = time() - Time
    if ErrXLam[end]>1
        println("not converge")
        return 0,zeros(nx,nx),[],[],[],1
    end

    return Time,Xistar,Index,IDCov_SC,COV_value_SC,0

end
