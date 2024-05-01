#* For now, I'm going to just focus on parameter estimation. This is an easier problem than random effects prediction
#* Later, I can introduce prediction of random effects

using Distributions, Random, Statistics
using LogExpFunctions   # For the logit and expit functions
using DataFrames
using MixedModels
using ProgressMeter, DrWatson
using JLD2
using GLM   # For linear and logistic regression with fixed effects

# Assume regression models for M and Y:
##  M = a_0 + a_1 * X + A_2 * W + d
## Y = b_0 + b_1 * M + b_2 * X + B_3 * W + e
## Note: In general, W is a matrix of confounders, so A_2 and B_3 are vectors of coefficients


Random.seed!(1)

num_reps = 1000    # Number of datasets to generate

n = 1000

a_0 = 1
a_1 = 1
a_2 = 1
sigma_d = 0.2   # Residual SD for M

b_0 = 1
b_1 = 1
b_2 = 1
b_3 = 1
sigma_e = 0.2   # Residual SD for Y




# X_dist = Bernoulli(0.5)
# W_dist = Bernoulli(0.5)

X_dist = Normal(1, 1)
W_dist = Normal(1, 1)



all_a_hats = []
all_b_hats = []

all_a_SEs = []
all_b_SEs = []

all_med_hats = []
all_med_SEs = []

# @showprogress Threads.@threads for _ in 1:M
@showprogress for _ in 1:num_reps

    
# Make X and W
X = rand(X_dist, n)
W = rand(W_dist, n)


    # Make data
    ## M
    delta = rand(Normal(0, sigma_d), n)
    M = a_0 .+ a_1 .* X .+ a_2 .* W .+ delta

    ## Y
    epsilon = rand(Normal(0, sigma_e), n)
    Y = b_0 .+ b_1 .* M .+ b_2 .* X .+ b_3 .* W .+ epsilon




    ## Fit M
    M_model = lm(@formula(M ~ X + W), DataFrame(M = M, X = X, W = W))
    a_hat = coef(M_model)
    a_SE = stderror(M_model)
    a_cov = vcov(M_model)

    push!(all_a_hats, a_hat)
    push!(all_a_SEs, a_SE)

    push!(all_a_hats, coef(M_model))
    push!(all_a_SEs, stderror(M_model))

    ## Fit Y
    Y_model = lm(@formula(Y ~ M + X + W), DataFrame(Y = Y, M = M, X = X, W = W))
    b_hat = coef(Y_model)
    b_SE = stderror(Y_model)
    b_cov = vcov(Y_model)

    push!(all_b_hats, b_hat)
    push!(all_b_SEs, b_SE)


    ## Mediation effect
    b2_hat = b_hat[3]
    b1_hat = b_hat[2]
    a1_hat = a_hat[2]
    
    ### Estimate
    med_hat = b2_hat + b1_hat * a1_hat

    push!(all_med_hats, med_hat)

    ### Standard error
    var_a1 = n * a_SE[2]^2
    var_b1 = n * b_SE[2]^2
    var_b2 = n * b_SE[3]^2
    cov_b1_b2 = n * b_cov[2, 3]

    med_var = (b1_hat^2 * var_a1 + var_b2 + a1_hat^2 * var_b1 + 2 * a1_hat * cov_b1_b2) / n
    med_SE = sqrt(med_var)

    push!(all_med_SEs, med_SE)
end



par_list_name = @savename num_reps n 
@save datadir("Cont-Resp, Cont-Med, Fixed-$par_list_name.jld2") all_a_hats all_a_SEs all_b_hats all_b_SEs all_med_hats all_med_SEs
println("Work complete!")






#* Multiple confounders


Random.seed!(1)

num_reps = 1000    # Number of datasets to generate

n = 1000
p_conf = 3      # Number of confounders

a_0 = 1
a_1 = 1
A_2 = repeat([1], p_conf)
sigma_d = 0.2   # Residual SD for M

b_0 = 1
b_1 = 1
b_2 = 1
B_3 = repeat([1], p_conf)
sigma_e = 0.2   # Residual SD for Y





# X_dist = Bernoulli(0.5)
# W_dist = Bernoulli(0.5)

X_dist = Normal(1, 1)
W_dist = Normal(1, 1)



all_a_hats = []
all_b_hats = []

all_a_SEs = []
all_b_SEs = []

all_med_hats = []
all_med_SEs = []

# @showprogress Threads.@threads for _ in 1:M
@showprogress for _ in 1:num_reps

    
    # Make X and W
    X = rand(X_dist, n)
    W = rand(W_dist, n, p_conf)
    
    # Make data
    ## M
    delta = rand(Normal(0, sigma_d), n)
    M = a_0 .+ a_1 .* X .+ W * A_2 .+ delta

    ## Y
    epsilon = rand(Normal(0, sigma_e), n)
    Y = b_0 .+ b_1 .* M .+ b_2 .* X .+ W * B_3 .+ epsilon




    ## Fit M
    M_data = DataFrame(M=M, X=X)
    for p in 1:p_conf
        M_data[!,Symbol("W$p")] = W[:, p]
    end

    M_formula = term(:M) ~ sum(term.(names(M_data, Not(:M))))

    M_model = lm(M_formula, M_data)
    a_hat = coef(M_model)
    a_SE = stderror(M_model)
    a_cov = vcov(M_model)

    push!(all_a_hats, a_hat)
    push!(all_a_SEs, a_SE)

    push!(all_a_hats, coef(M_model))
    push!(all_a_SEs, stderror(M_model))

    ## Fit Y
    Y_data = DataFrame(Y=Y, M=M, X=X)
    for p in 1:p_conf
        Y_data[!,Symbol("W$p")] = W[:, p]
    end

    Y_formula = term(:Y) ~ sum(term.(names(Y_data, Not(:Y))))

    Y_model = lm(Y_formula, Y_data)
    b_hat = coef(Y_model)
    b_SE = stderror(Y_model)
    b_cov = vcov(Y_model)

    push!(all_b_hats, b_hat)
    push!(all_b_SEs, b_SE)


    ## Mediation effect
    b2_hat = b_hat[3]
    b1_hat = b_hat[2]
    a1_hat = a_hat[2]
    
    ### Estimate
    med_hat = b2_hat + b1_hat * a1_hat

    push!(all_med_hats, med_hat)

    ### Standard error
    var_a1 = n * a_SE[2]^2
    var_b1 = n * b_SE[2]^2
    var_b2 = n * b_SE[3]^2
    cov_b1_b2 = n * b_cov[2, 3]

    med_var = (b1_hat^2 * var_a1 + var_b2 + a1_hat^2 * var_b1 + 2 * a1_hat * cov_b1_b2) / n
    med_SE = sqrt(med_var)

    push!(all_med_SEs, med_SE)
end



par_list_name = @savename num_reps n p_conf
@save datadir("Cont-Resp, Cont-Med, Fixed-$par_list_name.jld2") all_a_hats all_a_SEs all_b_hats all_b_SEs all_med_hats all_med_SEs
println("Work complete!")




