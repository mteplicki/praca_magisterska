struct SingleMachineDueDates
    n::Int
    p::Vector{Int}
    phat::Vector{Int}
    r::Vector{Int}
    d::Vector{Int}
    Γ::Int
end


function SingleMachineDueDates(n::Int, R::Float64, T::Float64, G::Float64)
    p_min = 1
    p_max = 100

    p = rand(p_min:p_max, n)

    phat = [rand(ceil(p_i * 2 /G):floor(p_i * 7 /G)) for p_i in p]

    P = sum(p)

    d = rand(ceil(P*(1-T-R/2)):floor(P*(1-T+R/2)), n)

    Γ = rand(ceil(2*10^-3*G*n):floor(18*10^-3*G*n))

    r = zeros(Int, n)

    return SingleMachineDueDates(n, p, phat, r, d, Γ)
    
end