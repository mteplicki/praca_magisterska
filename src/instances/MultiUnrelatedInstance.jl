struct MultiUnrelatedInstance <: AbstractInstance
    n::Int
    m::Int
    p::Matrix{Int}
    phat::Matrix{Int}
    Γ::Int
end

function MultiUnrelatedInstance(n::Int, m::Int, G::Float64)
    p_min = 1
    p_max = 100

    p = rand(p_min:p_max, m, n)

    phat = [rand(floor(p_i * 2 /G):ceil(p_i * 7 /G)) for p_i in p]

    Γ = rand(floor(5*10^-3*G*n):ceil(9*10^-3*G*n))

    return MultiUnrelatedInstance(n, m, p, phat, Γ)
end

function MultiUnrelatedInstance(n::Int, m::Int, G::Float64, p::Matrix{Int})
    phat = [rand(floor(p_i * 2 /G):ceil(p_i * 7 /G)) for p_i in p]

    Γ = rand(floor(5*10^-3*G*n):ceil(9*10^-3*G*n))

    return MultiUnrelatedInstance(n, m, p, phat, Γ)
end 