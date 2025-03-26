mutable struct SingleTardinessVariables <: AbstractVariableRef
    x::Matrix{VariableRef}
    t::Matrix{VariableRef}
    z::VariableRef
end

struct SingleTardiness <: AbstractColumnGenerationModel
    n::Int
    p::Vector{Int}
    phat::Vector{Int}
    d::Vector{Int}
    Γ::Int
    variables::SingleTardinessVariables
    Λ::Vector{BitVector}
    model::Model
end

SingleTardiness(instance::SingleMachineDueDates) = SingleTardiness(instance.n, instance.p, instance.phat, instance.d, instance.Γ)

function SingleTardiness(n::Int, p::Vector{Int}, phat::Vector{Int}, d::Vector{Int}, Γ::Int)
    model = Model(GLPK.Optimizer)

    Λ = [falses(n)]

    @variable(model, x[1:n, 1:n], Bin)

    @variable(model, t[1:n, 1:length(Λ)] >= 0)

    @variable(model, z >= 0)

    @constraint(model, [j in 1:n, (λ, δ) in enumerate(Λ)], t[j,λ] <= z)

    @constraint(model, [k in 1:n, (λ, δ) in enumerate(Λ)], 
        sum((p[i] + phat[i] * δ[i]) * sum(x[i,u] for u in 1:k) for i in 1:n) - sum(d[i]* x[i,k] for i in 1:n) <= t[k,λ])

    @constraint(model, [k in 1:n], sum(x[i,k] for i in 1:n) == 1)

    @constraint(model, [i in 1:n], sum(x[i,k] for k in 1:n) == 1)

    @objective(model, Min, z)

    return SingleTardiness(n, p, phat, d, Γ, SingleTardinessVariables(x,t,z), Λ, model)
end

function update_model!(model::SingleTardiness, new_Λ::Vector{BitVector})
    t_new = @variable(model.model, [1:model.n, (length(model.Λ) + 1):(length(model.Λ) + length(new_Λ))], lower_bound = 0, base_name = "t")
    # concatenate the new t variables to the old ones
    model.variables.t = hcat(model.variables.t, t_new)
    t = model.variables.t
    x = model.variables.x
    z = model.variables.z
    λ_old = length(model.Λ)
    @constraint(model.model, [j in 1:model.n, (λ, δ) in enumerate(new_Λ)], t[j,λ_old+λ] <= z)
    @constraint(model.model, [k in 1:model.n, (λ, δ) in enumerate(new_Λ)], 
        sum((model.p[i] + model.phat[i] * δ[i]) * sum(x[i,u] for u in 1:k) for i in 1:model.n) - sum(model.d[i]* x[i,k] for i in 1:model.n) <= t[k,λ_old+λ])
    append!(model.Λ, new_Λ)
    return model
    
end

max_zero(x...) = max(0, x...)

struct FMemStruct
    F::OffsetArrays.OffsetArray{Int64, 3, Array{Int64, 3}}
    r::Vector{Int}
    p::Vector{Int}
    phat::Vector{Int}
    d::Vector{Int}
    last_F::OffsetArrays.OffsetArray{Tuple{Int64, Int64, Int64}, 3, Array{Tuple{Int64, Int64, Int64}, 3}}
end

function F_mem(F_mem_struct::FMemStruct, i::Int, γ::Int, β::Int)
    F = F_mem_struct.F
    last_F = F_mem_struct.last_F
    r = F_mem_struct.r
    p = F_mem_struct.p
    phat = F_mem_struct.phat
    d = F_mem_struct.d

    value_δ0_1 = if F[i-1, γ, β + 1] != typemin(Int) F[i-1, γ, β + 1] + (β + 1) * p[i] - d[i] else typemin(Int) end
    value_δ0_2 = if F[i-1, γ, 0] != typemin(Int) F[i-1, γ, 0] + (β + 1) * r[i] + (β + 1) * p[i] - d[i] else typemin(Int) end
    value_δ0_3 = if F[i-1, γ, β] != typemin(Int) F[i-1, γ, β] + (β + 1) * r[i] + β * p[i] else typemin(Int) end
    value_δ0_4 = if F[i-1, γ, 0] != typemin(Int) F[i-1, γ, 0] + β * r[i] + β * p[i] else typemin(Int) end
    value_δ1_1 = if F[i-1, γ-1, β + 1] != typemin(Int) F[i-1, γ-1, β + 1] + (β + 1) * (p[i] + phat[i]) - d[i] else typemin(Int) end
    value_δ1_2 = if F[i-1, γ-1, 0] != typemin(Int) F[i-1, γ-1, 0] + (β + 1) * r[i] + (β + 1) * (p[i] + phat[i]) - d[i] else typemin(Int) end
    value_δ1_3 = if F[i-1, γ-1, β] != typemin(Int) F[i-1, γ-1, β] + (β + 1) * r[i] + β * (p[i] + phat[i]) else typemin(Int) end
    value_δ1_4 = if F[i-1, γ-1, 0] != typemin(Int) F[i-1, γ-1, 0] + β * r[i] + β * (p[i] + phat[i]) else typemin(Int) end
    value_table = [value_δ0_1, value_δ0_2, value_δ0_3, value_δ0_4, value_δ1_1, value_δ1_2, value_δ1_3, value_δ1_4]
    last_f_table = [(i-1, γ, β + 1), (i-1, γ, 0), (i-1, γ, β), (i-1, γ, 0),
                    (i-1, γ-1, β + 1), (i-1, γ-1, 0), (i-1, γ-1, β), (i-1, γ-1, 0)]
    max_index = argmax(value_table)
    last_F[i,γ,β] = last_f_table[max_index]
    F[i,γ,β] = value_table[max_index]
    return F[i,γ,β]
end

function oracle_subproblem(model::SingleTardiness, kwargs)
    x = value.(model.variables.x)
    r = zeros(Int, model.n)
    # x[i,k] = 1 if job i is assigned to position każdej

    # create a permutation of the jobs based on the x variables
    σ = zeros(Int, model.n)
    for i in 1:model.n
        for k in 1:model.n
            if x[i,k] == 1
                σ[k] = i
            end
        end
    end

    #check if in kwargs is a key "subproblem" and if it is equal to 1
    if haskey(kwargs, :subproblem_method) && kwargs[:subproblem] == 2
        return oracle_subproblem2(model, σ)
    else
        return oracle_subproblem1(model, σ)
    end
end

function oracle_subproblem1(model::SingleTardiness, σ::Vector{Int})
    p = model.p[σ]
    phat = model.phat[σ]
    d = model.d[σ]
    r = zeros(Int, model.n)


    F = zeros(Int, model.n, model.Γ + 2, model.n + 2)
    F = OffsetArray(F, 1:model.n, -1:model.Γ, 0:(model.n+1))

    fill!(F, typemin(Int))

    # set_F = falses(model.n, model.Γ + 1, model.n + 1)

    last_F = [(-1,0,-1) for _ in F]

    for γ in 0:1
        for β in 0:model.n
            F[1, γ, β] = max_zero(r[1] + p[1] + γ * phat[1] - d[1]) + β * (r[1] + p[1] + γ * phat[1])
            # set_F[1, γ, β] = true
        end
    end
    @show typeof(F)
    @show typeof(last_F)

    F_mem_struct = FMemStruct(F,#= set_F,=# r, p, phat, d, last_F)

    for i in 2:model.n
        for γ in 0:min(i, model.Γ)
            for β in 0:model.n
                F_mem(F_mem_struct, i, γ, β)
            end
        end
    end

    δ = falses(model.n)

    #pretty-print table F using show function, do not skip any rows, use type MIME
    show(IOContext(stdout, :limit => false), "text/plain", F)

    worst_value = F_mem_struct.F[model.n, model.Γ, 0]

    i = model.n
    last_case = (model.n, model.Γ, 0)
    while i > 0 && last_case != (-1,0,-1)
        i, γ, β = last_case
        last_case = last_F[i,γ,β]
        if γ - last_case[2] == 1
            δ[σ[i]] = true
        end
    end

    return worst_value, δ
end

function oracle_subproblem_2(model::SingleTardiness, σ::Vector{Int})
    n = model.n
    Γ = model.Γ
    
    p = model.p[σ]
    phat = model.phat[σ]
    d = model.d[σ]

    Φ_hat = sort(phat, rev = true)[1:Γ] |> sum

    α = zeros(Int, n, Γ + 1, Φ_hat + 1)

    α = OffsetArray(α, 1:n, 0:Γ, 0:Φ_hat)

    fill!(α, -Inf)

    set_α = falses(n, Γ + 1, Φ_hat + 1)

    set_α = OffsetArray(set_α, 1:n, 0:Γ, 0:Φ_hat)

    last_α = [(-1,0,-1) for _ in α]

    α[1,0,0] = max_zero(p[1] - d[1])
    set_α[1,0,0] = true

    for γ in 1:Γ
        α[1,γ,p[1]] = max_zero(p[1] + phat[1] - d[1])
        set_α[1,γ,0:Φ_hat] .= true
    end

    for k in 2:n
        α[k,0,0] = max_zero(sum(p[i] for i in 1:k)) - d[k] + α[k-1,0,0]
        set_α[k,0,0:Φ_hat] .= true
    end

    max_value = -Inf
    max_value_args = (-1, 0,-1)
    for k in 2:n
        k_last = k-1
        for γ in 1:model.Γ
            for Θ in 0:Φ_hat
                if Θ - p[k] < 0
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
                        @warn "coś nie zadziałało, odwoływanie się do nieistniejącego przypadku 2"
                    end
                    opt_1 = α[k-1,γ,Θ]
                    if !set_α[k-1,γ-1,Θ-p[k]]
                        @warn "coś nie zadziałało, odwoływanie się do nieistniejącego przypadku 3"
                    end
                    opt_2 = α[k-1,γ-1,Θ-p[k]]
                    if opt_1 != typemin(Int) || opt_2 != typemin(Int)
                        α[k,γ,Θ] = max_zero(Θ + sum(p[i] for i in 1:k) - d[k]) + max(opt_1, opt_2)
                        if opt_1 > opt_2
                            last_α[k,γ,Θ] = (k_last, γ, Θ)
                        else
                            last_α[k,γ,Θ] = (k_last, γ-1, Θ-p[k])
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

    max_value_args = (-1,0,-1)
    max_value = typemin(Int)

    for Φ in 0:Φ_hat
        if α[n,Γ,Φ] > max_value
            max_value = α[n,γ,Φ]
            max_value_args = (n,γ,Φ)
        end
    end

    i = n
    last_case = max_value_args
    while i > 0
        i, γ, α = last_case
        if last_V[i,γ,α] == (-1, 0,-1)
            @warn "breaking xd coś nie zadziałało"
            break
        end
        if γ - last_case[2] == 1
            δ[σ[i]] = true
        end
        last_case = last_V[i,γ,α]
    end

    return max_value, δ
end

