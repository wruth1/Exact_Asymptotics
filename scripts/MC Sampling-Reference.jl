#* For now, I'm going to just focus on parameter estimation. This is an easier problem than random effects prediction
#* Later, I can introduce prediction of random effects

using Distributions, Random, Statistics
using LogExpFunctions   # For the logit and expit functions
using DataFrames
using MixedModels
using ProgressMeter, DrWatson
using JLD2


Random.seed!(10000)

M = 1000    # Number of datasets to generate

all_beta_hats = []
all_fixef_SEs = []

n = 100
K = 5


beta = [-1, 1, -1]
sigma = 1


X_dist = Bernoulli(0.5)
Z_dist = Bernoulli(0.5)

# Make X and Z for each group
all_X_vecs = []
all_Z_vecs = []

for k = 1:K
    push!(all_X_vecs, rand(X_dist, n))
    push!(all_Z_vecs, rand(Z_dist, n))
end

@showprogress Threads.@threads for _ in 1:M

    
    # Generate random effects
    U = rand(Normal(0, sigma), K)


    # Compute linear predictors in each group

    all_eta_vecs = []

    for k = 1:K
        ## Fixed effects
        eta = beta[1] .+ beta[2] .* all_X_vecs[k] .+ beta[3] .* all_Z_vecs[k]

        ## Random effects
        eta += all_Z_vecs[k] .* U[k]

        push!(all_eta_vecs, eta)
    end


    # Generate responses

    ## Get probabilities
    all_p_vecs = []
    for k = 1:K
        push!(all_p_vecs, logistic.(all_eta_vecs[k]))
    end

    ## Generate responses
    all_Y_vecs = []
    for k = 1:K
        push!(all_Y_vecs, rand.(Bernoulli.(all_p_vecs[k])))
    end


    # Compile dataset

    ## Concatenate data across groups
    Y_data = vcat(all_Y_vecs...)
    X_data = vcat(all_X_vecs...)
    Z_data = vcat(all_Z_vecs...)
    group_data = repeat(1:K, inner = n)

    data = DataFrame(Y = Y_data, X = X_data, Z = Z_data, group = group_data)




    # Fit GLMM


    model_form = @formula(Y ~ X + Z + (0 + Z | group))      # The 0 + Z syntax is used to omit a random effect for the intercept
    model = fit(MixedModel, model_form, data, Bernoulli(), verbose=false, progress=false)


    # Extract and store fixed effects
    beta_hat = fixef(model)
    push!(all_beta_hats, beta_hat)

    # Extract and store fixed effects' standard errors
    fixef_SE = stderror(model)
    push!(all_fixef_SEs, fixef_SE)

end


par_list_name = @savename M n K beta sigma
@save datadir("Fixed Covariate Resampling-$par_list_name.jld2") all_beta_hats all_fixef_SEs
println("Work complete!")