mutable struct SingleTardyJobsVariables <: AbstractVariableRef
    Z::Matrix{VariableRef}
    S::Matrix{VariableRef}
    U::Matrix{VariableRef}
    y::VariableRef
end

struct SingleTardyJobsModel <: AbstractColumnGenerationModel
    instance::SingleMachineDueDates
    n::Int
    p::Vector{Int}
    phat::Vector{Int}
    d::Vector{Int}
    Γ::Int
    w::Vector{Int}
    Λ::Vector{BitVector}
    M::Int
    variables::SingleTardyJobsVariables
    model::Model
end

SingleTardyJobsModel(optimizer, instance::SingleMachineDueDates) = SingleTardyJobsModel(optimizer,instance.n, instance.p, instance.phat, instance.d, instance.Γ)

SingleTardyJobsModel(optimizer,n::Int, p::Vector{Int}, phat::Vector{Int}, d::Vector{Int}, Γ::Int) = SingleTardyJobsModel(optimizer,n, p, phat, d, Γ, ones(Int, n))

function SingleTardyJobsModel(optimizer, n::Int, p::Vector{Int}, phat::Vector{Int}, d::Vector{Int}, Γ::Int, w::Vector{Int})
    model = Model(optimizer)

    Λ = [falses(n)]

    # big M

    M = sum(p) + sum(phat)

    @variable(model, Z[1:n, 1:n], Bin)

    @variable(model, S[1:n, 1:length(Λ)] >= 0)

    @variable(model, U[1:n, 1:length(Λ)], Bin)

    @variable(model, y >= 0)

    @objective(model, Min, y)

    for i in 1:n
        for j in 1:n
            if i == j
                continue
            end
            @constraint(model, [(λ, δ) in enumerate(Λ)], S[j,λ] >= S[i,λ] + p[i] + phat[i] * δ[i] - M * (1 - Z[i,j]))
        end
    end

    for i in 1:n
        @constraint(model, [j in (i+1):n], Z[i,j] + Z[j,i] == 1)
    end

    @constraint(model, [j in 1:n, (λ, δ) in enumerate(Λ)], S[j,λ] + p[j] + phat[j] * δ[j] <= d[j] + M * U[j,λ])

    @constraint(model, [(λ, δ) in enumerate(Λ)], sum(w[j]*U[j,λ] for j in 1:n) <= y)

    return SingleTardyJobsModel(SingleMachineDueDates(n,p,phat,r,d,Γ),n, p, phat, d, Γ, w, Λ, M, SingleTardyJobsVariables(Z,S,U,y), model)
end


function update_model!(model::SingleTardyJobsModel, new_Λ::Vector{BitVector}, LB::Float64)
    S_new = @variable(model.model, [1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))], lower_bound = 0.0, base_name = "S")   
    U_new = @variable(model.model, [1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))], Bin, base_name = "U")

    model.variables.S = hcat(model.variables.S, S_new)
    model.variables.U = hcat(model.variables.U, U_new)
    S = model.variables.S
    U = model.variables.U
    Z = model.variables.Z
    y = model.variables.y


    λ = length(model.Λ) + 1

    for (δ) in new_Λ
        for i in 1:model.n
            for j in 1:model.n
                if i == j
                    continue
                end
                @constraint(model.model, S[j,λ] >= S[i,λ] + model.p[i] + model.phat[i] * δ[i] - model.M * (1 - Z[i,j]))
            end
        end

        @constraint(model.model, [j in 1:model.n], S[j,λ] + model.p[j] + model.phat[j] * δ[j] <= model.d[j] + model.M * U[j,λ])

        @constraint(model.model, sum(model.w[j]*U[j,λ] for j in 1:model.n) <= y)

        λ += 1
    end

    # @constraint(model.model, y >= LB)

    append!(model.Λ, new_Λ)
    
    return model
end

function find_job_permutation(x_matrix)
    n = size(x_matrix, 1)  # Number of jobs
    in_degree = Dict(j => 0 for j in 1:n)
    adj_list = Dict(i => Int[] for i in 1:n)

    # Build adjacency list and compute in-degrees
    for i in 1:n
        for j in 1:n
            if i != j && x_matrix[i, j]
                push!(adj_list[i], j)
                in_degree[j] += 1
            end
        end
    end

    # Topological sorting using Kahn's Algorithm
    queue = Queue{Int}()
    for j in 1:n
        if in_degree[j] == 0
            enqueue!(queue, j)
        end
    end

    job_sequence = Int[]

    while !isempty(queue)
        job = dequeue!(queue)
        push!(job_sequence, job)
        for neighbor in adj_list[job]
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0
                enqueue!(queue, neighbor)
            end
        end
    end

    return job_sequence  # The job permutation
end

function oracle_subproblem(model::SingleTardyJobsModel, permutation, kwargs)
    σ = permutation

    value_worst_case, δ = if haskey(kwargs, :subproblem_method) && kwargs[:subproblem_method] == 2
        oracle_subproblem2(model, σ)
    else
        oracle_subproblem1(model, σ)
    end

    return value_worst_case, δ
end

function oracle_subproblem1(model::SingleTardyJobsModel, σ::Vector{Int})


    @show σ

    @show "wssa"

    n = model.n
    Γ = model.Γ
    r = zeros(Int, n)
    d = model.d[σ]
    p = model.p[σ]
    phat = model.phat[σ]
    w = model.w[σ]

    V = zeros(Int, model.n, model.Γ + 2, sum(model.w) + 1)
    fill!(V, typemin(Int))
    V = OffsetArray(V, 1:model.n, 0:(model.Γ+1), 0:sum(model.w))
    last_V = [(-1,0,-1) for _ in V]

    # @show last_V

    # @show V

    max_alpha = (-1,0,-1)

    if r[1] + p[1] <= d[1]
        V[1,0,0] = r[1]
        max_alpha = (1,0,0)
    else
        V[1,0,0] = typemin(Int)
    end

    if r[1] + p[1] <= d[1] 
        V[1,0,1] = typemin(Int)
    else
        V[1,0,1] = r[1] + p[1]
        max_alpha = (1,0,1)
    end

    if r[1] + p[1] + phat[1] <= d[1]
        V[1,1,0] = r[1] + p[1] + phat[1]
    else
        V[1,1,0] = typemin(Int)
    end

    if r[1] + p[1] + phat[1] <= d[1]
        V[1,1,1] = typemin(Int)
    else
        V[1,1,1] = r[1] + p[1] + phat[1]
    end

    for i in 2:n
        V[i,0,0] = if V[i-1,0,0] != typemin(Int) && V[i-1,0,0] + p[i] <= d[i]
             max(r[i], V[i-1,0,0] + p[i])
        else
            typemin(Int)
        end 
    end

    for i in 2:n
        cum_arr = [0;cumsum(w[1:i])]
        for (α_before, α) in zip(cum_arr[1:end-1], cum_arr[2:end])
            V[i,0,α] = if V[i-1,0,α_before] + p[i] > d[i]
                V[i-1,0,α_before] + p[i]
            elseif V[i-1,0,α] != typemin(Int) && V[i-1,0,α] + p[i] <= d[i]
                max(r[i], V[i-1,0,α] + p[i])
            else
                typemin(Int)
            end
            if V[i,0,α] != (typemin(Int)) && α > max_alpha[3]
                max_alpha = (i,0,α)
            end
        end
    end

    for i in 2:n
        for γ in 1:min(i, model.Γ+1)
            if V[i-1,γ,0] != typemin(Int) && V[i-1,γ,0] + p[i] <= d[i] && V[i-1,γ-1,0] != typemin(Int) && V[i-1,γ-1,0] + p[i] + phat[i] <= d[i]
                opt_1 = r[i]
                opt_2 = V[i-1,γ,0] + p[i]
                opt_3 = V[i-1,γ-1,0] + p[i] + phat[i]
                V[i,γ,0] = max(opt_1, opt_2, opt_3)
                if opt_2 > opt_3
                    last_V[i,γ,0] = (i-1,γ,0)
                else
                    last_V[i,γ,0] = (i-1,γ-1,0)
                end
            elseif V[i-1,γ-1,0] != typemin(Int) && V[i-1,γ-1,0] + p[i] + phat[i] <= d[i]
                V[i,γ,0] = max(r[i], V[i-1,γ-1,0] + p[i] + phat[i])
                last_V[i,γ,0] = (i-1,γ-1,0)
            elseif V[i-1,γ,0] != typemin(Int) && V[i-1,γ,0] + p[i] <= d[i]
                V[i,γ,0] = max(r[i], V[i-1,γ,0] + p[i])
                last_V[i,γ,0] = (i-1,γ,0)
            else
                V[i,γ,0] = typemin(Int)
            end
        end
    end

    for i in 2:n
        for γ in 1:min(i, model.Γ+1)
            cum_arr = [0;cumsum(w[1:i])]
            for (α_before, α) in zip(cum_arr[1:end-1], cum_arr[2:end])
                if V[i-1,γ-1, α_before] + p[i] + phat[i] > d[i] || V[i-1,γ,α_before] + p[i] > d[i]
                    opt_1 = V[i-1,γ,α_before] + p[i]
                    opt_2 = V[i-1,γ-1,α_before] + p[i] + phat[i]
                    V[i,γ,α] = max(opt_1, opt_2)
                    if V[i,γ,α] == opt_1
                        last_V[i,γ,α] = (i-1,γ,α_before)
                    elseif V[i,γ,α] == opt_2
                        last_V[i,γ,α] = (i-1,γ-1,α_before)
                    end
                    if α > max_alpha[3]
                        max_alpha = (i,γ,α)
                    end
                elseif V[i-1,γ-1,α] != typemin(Int) && V[i-1,γ-1,α] + p[i] + phat[i] <= d[i] && V[i-1,γ,α] != typemin(Int) && V[i-1,γ,α] + p[i] <= d[i]
                    opt_1 = r[i]
                    opt_2 = V[i-1,γ,α] + p[i]
                    opt_3 = V[i-1,γ-1,α] + p[i] + phat[i]
                    V[i,γ,α] = max(opt_1, opt_2, opt_3)
                    if opt_2 > opt_3
                        last_V[i,γ,α] = (i-1,γ,α)
                    else
                        last_V[i,γ,α] = (i-1,γ-1,α)
                    end
                elseif V[i-1,γ-1,α] != typemin(Int) && V[i-1,γ-1,α] + p[i] + phat[i] <= d[i]
                    V[i,γ,α] = max(r[i], V[i-1,γ-1,α] + p[i] + phat[i])
                    last_V[i,γ,α] = (i-1,γ-1,α)
                elseif V[i-1,γ,α] != typemin(Int) && V[i-1,γ,α] + p[i] <= d[i]
                    V[i,γ,α] = max(r[i], V[i-1,γ,α] + p[i])
                    last_V[i,γ,α] = (i-1,γ,α)
                else
                    V[i,γ,α] = typemin(Int)
                end
            end
        end
    end

    worst_case = max_alpha
    @show worst_case
    value_worst_case = worst_case[3]
    δ = falses(n)
    i = n
    last_case = worst_case
    while i > 0 && last_case != (-1,0,-1)
        i, γ, α = last_case
        @show i
        @show last_case
        last_case = last_V[i,γ,α]
        if γ - last_case[2] == 1
            δ[σ[i]] = true
        end
    end

    return value_worst_case, δ
end

indicator_function(x) = if x > 0
    1
else
    0
end

function oracle_subproblem2(model::SingleTardyJobsModel, σ::Vector{Int})
    n = model.n
    Γ = model.Γ
    
    p = model.p[σ]
    phat = model.phat[σ]
    d = model.d[σ]

    Φ_hat = sort(phat, rev = true)[1:Γ] |> sum

    α = zeros(Int, n, Γ + 1, Φ_hat + 1)

    α = OffsetArray(α, 1:n, 0:Γ, 0:Φ_hat)

    fill!(α, typemin(Int))

    set_α = falses(n, Γ + 1, Φ_hat + 1)

    set_α = OffsetArray(set_α, 1:n, 0:Γ, 0:Φ_hat)

    last_α = [(-1,0,-1) for _ in α]

    α[1,0,0] = indicator_function(p[1] - d[1])
    set_α[1,0,0] = true

    for γ in 1:Γ
        α[1,γ,phat[1]] = indicator_function(p[1] + phat[1] - d[1])
        
    end
    set_α[1,0:Γ,0:Φ_hat] .= true

    for k in 2:n
        α[k,0,0] = indicator_function(sum(p[i] for i in 1:k) - d[k]) + α[k-1,0,0]
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
                        α[k,γ,Θ] = indicator_function(Θ + sum(p[i] for i in 1:k) - d[k]) + α[k-1,γ,Θ]
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
                        α[k,γ,Θ] = indicator_function(Θ + sum(p[i] for i in 1:k) - d[k]) + max(opt_1, opt_2)
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

function find_permutation(model::SingleTardyJobsModel)
    Z = value.(model.variables.Z)
    Z = Z .> 0.5
    n = model.n
    σ = zeros(Int, n)
    for i in 1:n
        for j in 1:n
            if Z[i,j] == 1
                σ[j] = i
            end
        end
    end
    return σ
end