
using Random, Distributions
using LinearAlgebra   # For the identity matrix, I()
using DataFrames, DataFramesMeta        # The latter works like dplyr in R with different verb names
using MixedModels, GLM
using ProgressMeter


# Set parameters

# N = 1000        # Total sample size
# K = 5           # Number of groups
# n = Int(N / K)       # Number of observations per group
p = 2           # Number of (mixed-effects) covariates

all_Ns = [500, 1000, 2000, 5000, 10000]
all_Ks = [5, 10, 25, 50, 100]

beta = [1.0, 2.0, 3.0]           # Fixed effects
Gamma = 0.2 * I(1 + p)      # Covariance matrix of random-effects

sigma = 1.0     # Standard deviation of the residuals


# Vector of all GLMM parameters
theta = [beta; sqrt.(diag(Gamma)); fill(0, 1 + p)]



Random.seed!(123)    


M = 100        # Number of Monte Carlo simulations

all_theta_hats = []

all_parameter_pairs = vec(collect(Iterators.product(all_Ns, all_Ks)))

@showprogress for (N, K) in all_parameter_pairs

    n = Int(N / K)       # Number of observations per group

    some_theta_hats = []

    for _ in 1:M

        # Generate data
        all_Xs = [randn(n, p) for _ in 1:K]

        ## Fixed-effects
        all_eta_vecs_fixed = [beta[1] .+ X * beta[2:end] for X in all_Xs]

        ## Add random effects
        all_bs = [rand(MvNormal(zeros(1 + p), Gamma)) for _ in 1:K]
        all_eta_vecs = [eta_vec_fixed .+ fill(b[1], n) .+ X * b[2:end] for (eta_vec_fixed, X, b) in zip(all_eta_vecs_fixed, all_Xs, all_bs)]

        ## Construct Y
        all_p_vecs = [1 ./ (1 .+ exp.(-eta_vec)) for eta_vec in all_eta_vecs]
        all_Ys = [rand.(Bernoulli.(p_vec)) for p_vec in all_p_vecs]

        ## Combine data from all groups
        X = vcat(all_Xs...)
        Y = vcat(all_Ys...)
        group = repeat(1:K, inner=n)

        data = DataFrame(Y = Y, X1 = X[:,1], X2 = X[:,2], group = group)



        # Fit the model
        model_formula = @formula(Y ~ 1 + X1 + X2 + (1 + X1 + X2 | group))

        model = fit(MixedModel, model_formula, data, Bernoulli(), verbose=false, progress=false)



        # Extract fitted parameters
        ## Fixed-effects
        beta_hat = coef(model)

        ## Random-effects' parameters
        RE_info = VarCorr(model).σρ[1]
        RE_SDs = setdiff(RE_info[1])
        RE_covs = setdiff(RE_info[2])


        theta_hat = [beta_hat; RE_SDs; RE_covs]
        push!(some_theta_hats, theta_hat)


    end

    push!(all_theta_hats, some_theta_hats)
end


# Get error in estimating theta using the 2-norm
all_err_norms = []

for i in 1:length(all_parameter_pairs)
    some_err_norms = []

    for j in 1:M
        theta_hat = all_theta_hats[i][j]
        if length(theta_hat) != length(theta) continue end
        # println("$i,$j - " * "$(length(theta_hat))")
        err = theta_hat - theta
        err_norm = norm(err, 2)

        push!(some_err_norms, err_norm)
    end
    push!(all_err_norms, some_err_norms)
end


# Summarize estimation error with mean and SD for each parameter configuration
all_err_means = [mean(errs) for errs in all_err_norms]
all_err_SDs = [std(errs) for errs in all_err_norms]



err_summary = DataFrame(N = get.(all_parameter_pairs, 1, -1), K = get.(all_parameter_pairs, 2, -1), mean_err = all_err_means, SD_err = all_err_SDs)

@transform!(err_summary, :n_per_group = Int.(:N ./ :K))



@chain err_summary begin
    @groupby(:N)
    @combine :mean_err = mean(:mean_err)
end

@chain err_summary begin
    @groupby(:K)
    @combine :mean_err = mean(:mean_err)
end


err_simp_reg_N = lm(@formula(log(mean_err) ~ log(N)), err_summary)
err_simp_reg_K = lm(@formula(log(mean_err) ~ log(K)), err_summary)
err_simp_reg_n = lm(@formula(log(mean_err) ~ log(n_per_group)), err_summary)

err_reg_N = lm(@formula(log(mean_err) ~ log(N) + log(K)), err_summary)
err_reg_n = lm(@formula(log(mean_err) ~ log(n_per_group) + log(K)), err_summary)