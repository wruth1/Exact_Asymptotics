

using Distributions, Random, Statistics
using LogExpFunctions   # For the logit and expit functions
using DataFrames
using MixedModels
using ProgressMeter, DrWatson
using JLD2
using GLM   # For linear and logistic regression with fixed effects




num_reps = 1000
n = 100
p_conf = 3

par_list_name = @savename num_reps n p_conf
@load datadir("Cont-Resp, Bin-Med, Fixed-$par_list_name.jld2") all_a_hats all_a_SEs all_b_hats all_b_SEs all_med_hats all_med_SEs




# Empirical SE
emp_SE = std(all_med_hats)


# Analytical SEs
mean_SE = mean(all_med_SEs)
median_SE = median(all_med_SEs)

std(all_med_SEs)
