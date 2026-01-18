function ASGD(c1,c2,Max_Iter,nx,X_true,Sigma,Xistar,sigma)
    ## Initialize variables
    t,X_t,barx_t = 1,zeros(nx),zeros(nx),zeros(nx) # iterates
    beta_t = 0 # stepsize
    diff_vec = [] # normality

    ## Start the loop
    Time = time()
    while t <= Max_Iter
        # Step 1: generate random gradient and Hessian
        a_t, eps_t = rand(MvNormal(zeros(nx),Sigma)), rand(Normal(0,sigma))
        barg_t = a_t*(a_t'*(X_t-X_true)) - eps_t*a_t

        # Step 2: update the iterate
        beta_t = c1/t^c2
        X_t = X_t - beta_t*barg_t
        barx_t = ((t-1)/t)*barx_t + (1/t)*X_t
        if mod(t, 50000) == 0
            diff_vec = push!(diff_vec, sqrt(t) * sum(barx_t-X_true) / sqrt(sum(Xistar)))
        end

        t = t + 1
    end
    Time = time()-Time

    return Time,diff_vec
end
