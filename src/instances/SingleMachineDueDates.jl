struct SingleMachineDueDates <: AbstractInstance
    n::Int
    p::Vector{Int}
    phat::Vector{Int}
    r::Vector{Int}
    d::Vector{Int}
    Γ::Int
    R::Float64
    T::Float64
    G::Float64
end


function SingleMachineDueDates(n::Int, R::Float64, T::Float64, G::Float64)
    p_min = 1
    p_max = 100

    p = rand(p_min:p_max, n)

    phat = [rand(ceil(p_i * 2 /G):ceil(p_i * 7 /G)) for p_i in p]

    P = sum(p)

    d = rand(max(ceil(P*(1-T-R/2)),0):ceil(P*(1-T+R/2)), n)

    Γ = rand(ceil(5*10^-3*G*n):ceil(9*10^-3*G*n))

    r = zeros(Int, n)

    return SingleMachineDueDates(n, p, phat, r, d, Γ, R, T, G)
end

function summary(instance::SingleMachineDueDates)
    return "SingleMachineDueDates(n=$(instance.n), r=$(instance.r), d=$(instance.d), Γ=$(instance.Γ), R=$(instance.R), T=$(instance.T), G=$(instance.G))"
end

function Base.show(io::IO, instance::SingleMachineDueDates)
    print(io, summary(instance))
end

