
# Compute mediation effect

## Increment in conditional expectation of M for a unit increase in X
function get_delta(eta, a_x)
    Q1 = exp(-eta)
    Q2 = (1 - exp(-a_x))
    Q3 = (1 + exp(-eta - a_x))
    Q4 = (1 + exp(-eta))

    return Q1 * Q2 / (Q3 * Q4)
end

## Increment in conditional expectation of Y for a unit increase in X
## I.e. Mean over M of increment in conditional expectation for Y given M and X
function get_gamma(delta, b_m, b_x)
    return b_m * delta + b_x
end


# Compute SE of estimated mediation effect
# We use the delta method, so many derivatives are required


## Helpful intermedeate quantities
function d_delta_d_eta(eta, a_x)
    delta = get_delta(eta, a_x)

    Q1 = 1 - exp(-2*eta - a_x)
    Q2 = 1 + exp(-eta)
    Q3 = 1 + exp(-eta - a_x)

    return (-1) * delta * Q1 / (Q2 * Q3)
end

# Old, incorrect version
# function d_delta_d_a_x(eta, a_x)
#     Q1 = exp(-2*eta - a_x)
#     Q2 = 1 + exp(-eta)
#     Q3 = 1 + exp(-eta - a_x)

#     return Q1 / (Q2 * Q3^2)
# end

# Alternative expression
# function d_delta_d_a_x(eta, a_x)
#     delta = get_delta(eta, a_x)

#     Q1 = exp(-a_x)
#     Q2 = 1 + exp(-eta)
#     Q3 = 1 + exp(-eta - a_x)
#     Q4 = 1 - exp(-a_x)

#     return delta * Q1 * Q2 / (Q3 * Q4)
# end

function d_delta_d_a_x(eta, a_x)
    Q1 = exp(-2*eta - a_x)
    Q2 = exp(-eta - a_x)
    Q3 = 1 + exp(-eta)
    Q4 = 1 + exp(-eta - a_x)

    return (Q1 + Q2)/(Q3 * Q4^2)
end


## Derivative of gamma wrt each parameter
function d_gamma_d_a_0(dd_de, b_m)
    return dd_de * b_m
end

function d_gamma_d_a_1(dd_de, dd_da_x, x, b_m)
    return b_m * x * dd_de + b_m * dd_da_x
end

function d_gamma_d_A_2(dd_de, W, b_m)
    return dd_de * b_m .* W
end


function d_gamma_d_b_0()
    return 0
end

function d_gamma_d_b_1(eta, a_x)
    return get_delta(eta, a_x)
end

function d_gamma_d_b_2()
    return 1
end

function d_gamma_d_B_3(W)
    return 0 .* W
end



## Gradient of gamma wrt all parameters
function d_gamma_d_theta(eta, x, W, a, b)
    a_0, a_x = a[1:2]
    A_2 = a[3:end]

    b_0, b_m, b_x = b[1:3]
    B_3 = b[4:end]

    dd_de = d_delta_d_eta(eta, a_x)
    dd_da_x = d_delta_d_a_x(eta, a_x)


    dg_da_0 = d_gamma_d_a_0(dd_de, b_m)
    dg_da_1 = d_gamma_d_a_1(dd_de, dd_da_x, x, b_m)
    dg_dA_2 = d_gamma_d_A_2(dd_de, W, b_m)

    dg_db_0 = d_gamma_d_b_0()
    dg_db_1 = d_gamma_d_b_1(eta, a_x)
    dg_db_2 = d_gamma_d_b_2()
    dg_dB_3 = d_gamma_d_B_3(W)

    return [dg_da_0; dg_da_1; dg_dA_2; dg_db_0; dg_db_1; dg_db_2; dg_dB_3]
end

