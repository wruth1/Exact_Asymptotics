
function expit(x)
    return 1 / (1 + exp(-x))
end

function logit(p)
    return log(p / (1 - p))
end