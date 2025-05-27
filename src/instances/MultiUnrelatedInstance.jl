struct MultiUnrelatedInstance <: AbstractInstance
    n::Int
    m::Int
    p::Matrix{Int}
    phat::Matrix{Int}
    Γ::Int
    p_min::Int
    p_max::Int
    G::Float64
end

function MultiUnrelatedInstance(n::Int, m::Int, G::Float64, p_min=1, p_max=100)

    p = rand(p_min:p_max, m, n)

    phat = [rand(floor(p_i * 2 /G):ceil(p_i * 7 /G)) for p_i in p]

    Γ = rand(floor(5*10^-3*G*n):ceil(9*10^-3*G*n))

    return MultiUnrelatedInstance(n, m, p, phat, Γ, p_min, p_max, G)
end

function MultiUnrelatedInstance(n::Int, m::Int, p::Matrix{Int}, G::Float64)
    p_min = minimum(p)
    p_max = maximum(p)

    phat = [rand(floor(p_i * 2 /G):ceil(p_i * 7 /G)) for p_i in p]

    Γ = rand(floor(5*10^-3*G*n):ceil(9*10^-3*G*n))

    return MultiUnrelatedInstance(n, m, p, phat, Γ, p_min, p_max)
end 

function summary(instance::MultiUnrelatedInstance)
    return "MultiUnrelatedInstance(n=$(instance.n), m=$(instance.m), Γ=$(instance.Γ), p_min=$(instance.p_min), p_max=$(instance.p_max), G=$(instance.G))"
end

function Base.show(io::IO, instance::MultiUnrelatedInstance)
    print(io, summary(instance))
end
