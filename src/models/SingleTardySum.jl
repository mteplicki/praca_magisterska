mutable struct SingleTardinessVariables <: AbstractVariableRef
    x::Matrix{VariableRef}
    t::Matrix{VariableRef}
    z::VariableRef
end

struct SingleTardiness <: AbstractColumnGenerationModel
    instance::SingleMachineDueDates
    n::Int
    p::Vector{Int}
    phat::Vector{Int}
    d::Vector{Int}
    Γ::Int
    variables::SingleTardinessVariables
    Λ::Vector{BitVector}
    model::Model
end

SingleTardiness(optimizer, instance::SingleMachineDueDates) = SingleTardiness(optimizer, instance.n, instance.p, instance.phat, instance.d, instance.Γ)

function SingleTardiness(optimizer, n::Int, p::Vector{Int}, phat::Vector{Int}, d::Vector{Int}, Γ::Int)
    model = Model(optimizer)

    Λ = [falses(n)]

    @variable(model, x[1:n, 1:n], Bin)

    @variable(model, t[1:n, 1:length(Λ)] >= 0)

    @variable(model, z >= 0)

    @constraint(model, [(λ, δ) in enumerate(Λ)], sum(t[j,λ] for j in 1:n) <= z)

    @constraint(model, [k in 1:n, (λ, δ) in enumerate(Λ)], 
        sum((p[i] + phat[i] * δ[i]) * sum(x[i,u] for u in 1:k) for i in 1:n) - sum(d[i]* x[i,k] for i in 1:n) <= t[k,λ])

    @constraint(model, [k in 1:n], sum(x[i,k] for i in 1:n) == 1)

    @constraint(model, [i in 1:n], sum(x[i,k] for k in 1:n) == 1)

    @objective(model, Min, z)

    return SingleTardiness(SingleMachineDueDates(n,p,phat,r,d,Γ),n, p, phat, d, Γ, SingleTardinessVariables(x,t,z), Λ, model)
end

function update_model!(model::SingleTardiness, new_Λ::Vector{BitVector}, LB::Float64)
    t_new = @variable(model.model, [1:model.n, (length(model.Λ) + 1):(length(model.Λ) + length(new_Λ))], lower_bound = 0, base_name = "t")
    # concatenate the new t variables to the old ones
    model.variables.t = hcat(model.variables.t, t_new)
    t = model.variables.t
    x = model.variables.x
    z = model.variables.z
    λ_old = length(model.Λ)
    @constraint(model.model, [(λ, δ) in enumerate(new_Λ)], sum(t[j,λ_old+λ] for j in 1:model.n) <= z)
    @constraint(model.model, [k in 1:model.n, (λ, δ) in enumerate(new_Λ)], 
        sum((model.p[i] + model.phat[i] * δ[i]) * sum(x[i,u] for u in 1:k) for i in 1:model.n) - sum(model.d[i]* x[i,k] for i in 1:model.n) <= t[k,λ_old+λ])
    append!(model.Λ, new_Λ)
    # @constraint(model.model, z >= LB)
    return model
end

max_zero(x...) = max(0, x...)

struct FMemStruct
    F::OffsetArrays.OffsetArray{Int64, 3, Array{Int64, 3}}
    set_F
    r::Vector{Int}
    p::Vector{Int}
    phat::Vector{Int}
    d::Vector{Int}
    last_F::OffsetArrays.OffsetArray{Tuple{Int64, Int64, Int64}, 3, Array{Tuple{Int64, Int64, Int64}, 3}}
end

function F_mem(F_mem_struct::FMemStruct, i::Int, γ::Int, β::Int)
    F = F_mem_struct.F
    set_F = F_mem_struct.set_F
    last_F = F_mem_struct.last_F
    r = F_mem_struct.r
    p = F_mem_struct.p
    phat = F_mem_struct.phat
    d = F_mem_struct.d

    value_δ0_1 = if F[i-1, γ, β + 1] != typemin(Int) F[i-1, γ, β + 1] + (β + 1) * p[i] - d[i] else typemin(Int) end
    value_δ0_2 = if F[i-1, γ, 0] != typemin(Int) F[i-1, γ, 0] + (β + 1) * r[i] + (β + 1) * p[i] - d[i] else typemin(Int) end
    value_δ0_3 = if F[i-1, γ, β] != typemin(Int) F[i-1, γ, β] + β * p[i] else typemin(Int) end
    value_δ0_4 = if F[i-1, γ, 0] != typemin(Int) F[i-1, γ, 0] + β * r[i] + β * p[i] else typemin(Int) end
    value_δ1_1 = if F[i-1, γ-1, β + 1] != typemin(Int) F[i-1, γ-1, β + 1] + (β + 1) * (p[i] + phat[i]) - d[i] else typemin(Int) end
    value_δ1_2 = if F[i-1, γ-1, 0] != typemin(Int) F[i-1, γ-1, 0] + (β + 1) * r[i] + (β + 1) * (p[i] + phat[i]) - d[i] else typemin(Int) end
    value_δ1_3 = if F[i-1, γ-1, β] != typemin(Int) F[i-1, γ-1, β] + β * (p[i] + phat[i]) else typemin(Int) end
    value_δ1_4 = if F[i-1, γ-1, 0] != typemin(Int) F[i-1, γ-1, 0] + β * r[i] + β * (p[i] + phat[i]) else typemin(Int) end
    value_table = [value_δ0_1, value_δ0_2, value_δ0_3, value_δ0_4, value_δ1_1, value_δ1_2, value_δ1_3, value_δ1_4]
    last_f_table = [(i-1, γ, β + 1), (i-1, γ, 0), (i-1, γ, β), (i-1, γ, 0),
                    (i-1, γ-1, β + 1), (i-1, γ-1, 0), (i-1, γ-1, β), (i-1, γ-1, 0)]
    max_index = argmax(value_table)
    last_F[i,γ,β] = last_f_table[max_index]
    F[i,γ,β] = value_table[max_index]
    return F[i,γ,β]
end

function check_worst_case(model::SingleTardiness, σ::Vector{Int}, δ::BitVector, expected_value::Int)
    T = []
    C = [0]
    p = model.p
    phat = model.phat
    d = model.d
    for i in 1:model.n
        if δ[σ[i]]
            push!(C, C[end] + p[σ[i]] + phat[σ[i]])
        else
            push!(C, C[end] + p[σ[i]])
        end
        push!(T, max_zero(C[end] - d[σ[i]]))
    end

    sum_T = sum(T)

    C = C[2:end]
    
    if sum_T != expected_value
        @warn "Expected value $expected_value, got $sum_T"
    end
end

function oracle_subproblem(model::SingleTardiness, permutation, kwargs)
    r = zeros(Int, model.n)

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

function oracle_subproblem_single_tardiness1(n, Γ, p, phat, d, σ::Vector{Int})
    p = p[σ]
    phat = phat[σ]
    d = d[σ]
    r = zeros(Int, n)


    F = zeros(Int, n, Γ + 2, n + 2)
    F = OffsetArray(F, 1:n, -1:Γ, 0:(n+1))

    fill!(F, typemin(Int))

    set_F = falses(n, Γ + 2, n + 2)

    set_F = OffsetArray(set_F, 1:n, -1:Γ, 0:(n+1))

    last_F = [(-1,0,-1) for _ in F]

    for γ in 0:1
        for β in 0:n
            F[1, γ, β] = max_zero(r[1] + p[1] + γ * phat[1] - d[1]) + β * (r[1] + p[1] + γ * phat[1])
            # set_F[1, γ, β] = true
        end
    end
    # @show typeof(F)
    # @show typeof(last_F)

    F_mem_struct = FMemStruct(F, set_F, r, p, phat, d, last_F)

    for i in 2:n
        for γ in 0:min(i, Γ)
            for β in 0:n
                F_mem(F_mem_struct, i, γ, β)
            end
        end
    end

    δ = falses(n)

    #pretty-print table F using show function, do not skip any rows, use type MIME
    # show(IOContext(stdout, :limit => false), "text/plain", F)

    worst_value = F_mem_struct.F[n, Γ, 0]

    i = n
    last_case = (n, Γ, 0)
    while i > 0 && last_case != (-1,0,-1)
        i, γ, β = last_case
        last_case = last_F[i,γ,β]
        if γ - last_case[2] == 1
            δ[σ[i]] = true
        end
    end

    return worst_value, δ
end

function oracle_subproblem_single_tardiness2(n, Γ, p, phat, d, σ::Vector{Int})
    
    p = p[σ]
    phat = phat[σ]
    d = d[σ]

    Φ_hat = sort(phat, rev = true)[1:Γ] |> sum

    α = zeros(Int, n, Γ + 1, Φ_hat + 1)

    α = OffsetArray(α, 1:n, 0:Γ, 0:Φ_hat)

    fill!(α, typemin(Int))

    set_α = falses(n, Γ + 1, Φ_hat + 1)

    set_α = OffsetArray(set_α, 1:n, 0:Γ, 0:Φ_hat)

    last_α = [(-1,0,-1) for _ in α]

    α[1,0,0] = max_zero(p[1] - d[1])
    set_α[1,0,0] = true

    for γ in 1:Γ
        α[1,γ,phat[1]] = max_zero(p[1] + phat[1] - d[1])
        
    end
    set_α[1,0:Γ,0:Φ_hat] .= true

    for k in 2:n
        α[k,0,0] = max_zero(sum(p[i] for i in 1:k) - d[k]) + α[k-1,0,0]
        set_α[k,0,0:Φ_hat] .= true
    end

    max_value = -Inf
    max_value_args = (-1, 0,-1)
    for k in 2:n
        k_last = k-1
        for γ in 1:model.Γ
            for Θ in 0:Φ_hat
                if Θ - phat[k] < 0
                    if !set_α[k-1,γ,Θ]
                        @warn "coś nie zadziałało, odwoływanie się do nieistniejącego przypadku 1"
                    end
                    if α[k-1,γ,Θ] != typemin(Int)
                        α[k,γ,Θ] = max_zero(Θ + sum(p[i] for i in 1:k) - d[k]) + α[k-1,γ,Θ]
                        last_α[k,γ,Θ] = (k_last, γ, Θ)
                    else
                        α[k,γ,Θ] = typemin(Int)
                    end
                else
                    if !set_α[k-1,γ,Θ]
                        @warn "coś nie zadziałało, odwoływanie się do nieistniejącego przypadku 2: $(k-1), $γ, $Θ"
                    end
                    opt_1 = α[k-1,γ,Θ]
                    if !set_α[k-1,γ-1,Θ-phat[k]]
                        @warn "coś nie zadziałało, odwoływanie się do nieistniejącego przypadku 3: $(k-1), $(γ-1), $(Θ-phat[k])"
                    end
                    opt_2 = α[k-1,γ-1,Θ-phat[k]]
                    if opt_1 != typemin(Int) || opt_2 != typemin(Int)
                        α[k,γ,Θ] = max_zero(Θ + sum(p[i] for i in 1:k) - d[k]) + max(opt_1, opt_2)
                        if opt_1 > opt_2
                            last_α[k,γ,Θ] = (k_last, γ, Θ)
                        else
                            last_α[k,γ,Θ] = (k_last, γ-1, Θ-phat[k])
                        end
                    else
                        α[k,γ,Θ] = typemin(Int)
                    end
                end
                set_α[k,γ,Θ] = true
            end
        end
    end

    δ = falses(n)

    # find worst last_case

    @show all(set_α)

    max_value_args = (-1,0,-1)
    max_value = typemin(Int)

    for Φ in 0:Φ_hat
        if α[n,Γ,Φ] > max_value
            max_value = α[n,Γ,Φ]
            max_value_args = (n,Γ,Φ)
        end
    end

    i = n
    last_case = max_value_args
    while i > 0 && last_case != (-1,0,-1)
        i, γ, α = last_case
        last_case = last_α[i,γ,α]
        if γ - last_case[2] == 1
            δ[σ[i]] = true
        end
    end

    if i != 1
        @warn "Something went wrong, i = $i"
    end

    return max_value, δ
end

function find_permutation(model::SingleTardiness)
    x = value.(model.variables.x)
    x = x .> 0.5
    # create a permutation of the jobs based on the x variables
    σ = zeros(Int, model.n)
    for i in 1:model.n
        for k in 1:model.n
            if x[i,k] == 1
                σ[k] = i
            end
        end
    end
    return σ
end

save_permutation_variable(model::SingleTardiness) = value.(model.variables.x)
set_start_value_for_model(model::SingleTardiness, permutation_variable) = begin
    if !isnothing(permutation_variable)
        set_start_value.(model.variables.x, permutation_variable)
    end
end