mutable struct SingleTardinessDominanceRulesVariables <: AbstractVariableRef
    x::Vector
    t::Matrix{VariableRef}
    z::VariableRef
    feasible_substraction_set::Vector{Vector{Int}}
    A::Vector{Vector{Int}}
    B::Vector{Vector{Int}}
end

struct SingleTardinessDominanceRules <: AbstractColumnGenerationModel
    n::Int
    p::Vector{Int}
    phat::Vector{Int}
    d::Vector{Int}
    Γ::Int
    variables::SingleTardinessDominanceRulesVariables
    Λ::Vector{BitVector}
    model::Model
end

function domination(i, j, E_j, L_i, p, phat, d, Γ)
    if p[i] + phat[i] < p[j] && d[i] <= max(E_j, d[j])
        return true
    elseif d[j] >= max(L_i - p[j], d[i])
        return true
    elseif d[j] >= L_i
        return true
    else
        return false
    end
end

E_func(i, B_i, p) = sum(p[B_i]) + p[i]

function L_func(i, n, A_i, p, phat, Γ)
    N_minus_A_i = setdiff(1:n, A_i)
    sorted_phat = sort([phat[k] for k in N_minus_A_i]; rev=true)
    Γ_effective = min(Γ, length(sorted_phat))
    return sum(p[N_minus_A_i]) + sum(sort(phat[N_minus_A_i]; rev=true)[1:Γ_effective])
end

function compute_dominated_sets(n, p, phat, d, Γ)
    A = [[] for _ in 1:n]
    B = [[] for _ in 1:n]
    
    changed = true
    while changed
        changed = false
        for i in 1:n
            Ei = E_func(i, B[i], p)
            Li = L_func(i, n, A[i], p, phat, Γ)
            
            for j in 1:n
                if i == j 
                    continue 
                end
                @show i,j
                
                Ej = E_func(j, B[j], p)
                Lj = L_func(j, n, A[j], p, phat, Γ)
                
                if domination(i, j, Ej, Li, p, phat, d, Γ) && (!(i in A[j]) || !(j in B[i]))
                    if !(j in A[i])
                        push!(A[i], j)
                        changed = true
                    end
                    if !(i in B[j])
                        push!(B[j], i)
                        changed = true
                    end
                end
            end
        end
    end
    return A, B
end
"""
Creates a SingleTardinessDominanceRules model from an instance.
- `optimizer`: The optimizer to use (e.g., `GLPK.Optimizer`).
- `instance`: An instance of `SingleMachineDueDates`.
"""
SingleTardinessDominanceRules(optimizer, instance::SingleMachineDueDates) = SingleTardinessDominanceRules(optimizer, instance.n, instance.p, instance.phat, instance.d, instance.Γ)

function SingleTardinessDominanceRules(optimizer, n::Int, p::Vector{Int}, phat::Vector{Int}, d::Vector{Int}, Γ::Int)
    model = Model(optimizer)

    Λ = [falses(n)]

    A, B = compute_dominated_sets(n, p, phat, d, Γ)

    x = [@variable(model, [(length(B[i]) + 1):(n-length(A[i]))], base_name = "x") for i in 1:n]

    @variable(model, t[1:n, 1:length(Λ)] >= 0)

    @variable(model, z >= 0)

    @constraint(model, [(λ, δ) in enumerate(Λ)], sum(t[j,λ] for j in 1:n) <= z)

    feasible_substraction_set = [[i for i in 1:n if length(B[i]) < k && k <= n-length(A[i])] for k in 1:n]

    @constraint(model, [k in 1:n, (λ, δ) in enumerate(Λ)], 
        sum((p[i] + phat[i] * δ[i]) * sum(x[i][u] for u in (length(B[i])+1):min(k, n-length(A[i]))) for i in 1:n) - sum(d[i]* x[i][k] for i in feasible_substraction_set[k]) <= t[k,λ])

    @constraint(model, [i in 1:n], sum(x[i][k] for k in (length(B[i])+1):(n-length(A[i]))) == 1)

    @constraint(model, [k in 1:n], sum(x[i][k] for i in feasible_substraction_set[k]) == 1)

    for i in 1:n
        @constraint(model, [j in A[i]], sum(k*x[j][k] for k in (length(B[j])+1):(n-length(A[j]))) >= sum(k*x[i][k] for k in (length(B[i])+1):(n-length(A[i]))) + 1)
    end
    @objective(model, Min, z)

    return SingleTardinessDominanceRules(n, p, phat, d, Γ, SingleTardinessDominanceRulesVariables(x,t,z, feasible_substraction_set, A, B), Λ, model)
end

function update_model!(model::SingleTardinessDominanceRules, new_Λ::Vector{BitVector}, LB::Float64)
    t_new = @variable(model.model, [1:model.n, (length(model.Λ) + 1):(length(model.Λ) + length(new_Λ))], lower_bound = 0, base_name = "t")
    # concatenate the new t variables to the old ones
    model.variables.t = hcat(model.variables.t, t_new)
    t = model.variables.t
    x = model.variables.x
    z = model.variables.z
    d = model.d
    p = model.p
    phat = model.phat
    feasible_substraction_set = model.variables.feasible_substraction_set
    n = model.n
    A = model.variables.A
    B = model.variables.B
    λ_old = length(model.Λ)
    @constraint(model.model, [(λ, δ) in enumerate(new_Λ)], sum(t[j,λ_old+λ] for j in 1:model.n) <= z)
    @constraint(model.model, [k in 1:model.n, (λ, δ) in enumerate(new_Λ)], 
    sum((p[i] + phat[i] * δ[i]) * sum(x[i][u] for u in length(B[i]):min(k, n-length(A[i]))) for i in 1:n) - sum(d[i]* x[i][k] for i in feasible_substraction_set[k]) <= t[k,λ_old+λ])
    append!(model.Λ, new_Λ)
    # @constraint(model.model, z >= LB)
    return model
end

function oracle_subproblem(model::SingleTardinessDominanceRules, permutation, kwargs)

    r = zeros(Int, model.n)
    A = model.variables.A
    B = model.variables.B
    n = model.n

    # create a permutation of the jobs based on the x variables
    σ = permutation

    #check if in kwargs is a key "subproblem" and if it is equal to 1
    worst_value, δ = if haskey(kwargs, :subproblem_method) && kwargs[:subproblem_method] == 2
        oracle_subproblem_single_tardiness2(model.n, model.Γ, model.p, model.phat, model.d ,σ)
    else
        oracle_subproblem_single_tardiness1(model.n, model.Γ, model.p, model.phat, model.d ,σ)
    end

    if haskey(kwargs, :subproblem_method) && kwargs[:subproblem_method] == 2
        worst_value_second, δ_second = oracle_subproblem_single_tardiness1(model.n, model.Γ, model.p, model.phat, model.d ,σ)
        if worst_value_second != worst_value
            @warn "Expected value $worst_value, got $worst_value_second"
            @info "δ: $δ"
            @info "δ_second: $δ_second"
        end
    end

    check_worst_case(model, σ, δ, worst_value)

    return worst_value, δ
    
end

function find_permutation(model::SingleTardinessDominanceRules)
    # find the permutation of the jobs based on the x variables
    σ = zeros(Int, model.n)

    x = value.(model.variables.x)

    x = x .> 0.5

    for i in 1:model.n
        for k in (length(model.variables.B[i])+1):(model.n-length(model.variables.A[i]))
            if x[i][k] == 1
                σ[k] = i
            end
        end
    end
    return σ
end