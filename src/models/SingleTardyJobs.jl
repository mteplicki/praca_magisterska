struct SingleTardyJobsModel <: AbstractColumnGenerationModel
    n::Int
    p::Vector{Int}
    phat::Vector{Int}
    r::Vector{Int}
    d::Vector{Int}
    Γ::Int
    w::Vector{Int}
    Z::Vector{VariableRef}
    Λ::Vector{BitVector}
    model::Model
end

SingleTardyJobsModel(n::Int, p::Vector{Int}, phat::Vector{Int}, r::Vector{Int}, d::Vector{Int}, Γ::Int) = SingleTardyJobsModel(n, p, phat, r, d, Γ, ones(Int, n))

function SingleTardyJobsModel(n::Int, p::Vector{Int}, phat::Vector{Int}, r::Vector{Int}, d::Vector{Int}, Γ::Int, w::Vector{Int})
    model = Model(GLPK.Optimizer)

    Λ = [falses(n)]

    # big M

    M = sum(p) + sum(phat) + sum(r)

    @variable(model, Z[1:n, 1:n], Bin) 

    @variable(model, S[1:n, 1:length(Λ)] >= 0)

    @variable(model, U[1:n, 1:length(Λ)], Bin)

    @variable(model, y >= 0)

    @objective(model, Min, y)

    @constraint(model, [j in 1:n, (λ, δ) in enumerate(Λ)], S[j,λ] >= r[j])

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

    return SingleTardyJobsModel(n, p, phat, r, d, Γ, w, Z, Λ, model)

end


function update_model!(model::SingleTardyJobsModel, new_Λ::Vector{BitVector})
    @variable(model.model, S[1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))] >= 0)
    @variable(model.model, U[1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))], Bin)

    λ = length(model.Λ) + 1

    for (δ) in new_Λ
        @constraint(model.model, [j in 1:model.n], S[j,λ] >= model.r[j])
        for i in 1:model.n
            for j in 1:model.n
                if i == j
                    continue
                end
                @constraint(model.model, S[j,λ] >= S[i,λ] + model.p[i] + model.phat[i] * δ[i] - M * (1 - model.Z[i,j]))
            end
        end

        @constraint(model.model, [j in 1:model.n], S[j,λ] + model.p[j] + model.phat[j] * δ[j] <= model.d[j] + M * U[j,λ])

        @constraint(model.model, sum(model.w[j]*U[j,λ] for j in 1:model.n) <= model.y)

        λ += 1
    end
    
    return model
end

function oracle_subproblem(model::SingleTardyJobsModel, Z::Matrix{Int})
    σ = [findfirst(Z[i,:]) for i in 1:model.n]

    n = model.n
    Γ = model.Γ
    r = model.r[σ]
    d = model.d[σ]
    p = model.p[σ]
    phat = model.phat[σ]
    w = model.w[σ]

    V = zeros(Int, model.n, model.Γ + 1, sum(model.w) + 1)
    fill!(V, -Inf)
    V = OffsetArray(V, 1:model.n, 0:model.Γ, 0:sum(model.w))

    V[1,0,0] = if r[1] + p[1] <= d[1]
        r[1]
    else
        -Inf
    end

    V[1,0,1] = if r[1] + p[1] <= d[1] 
        -Inf
    else
        r[1] + p[1]
    end

    V[1,1,0] = if r[1] + p[1] + phat[1] <= d[1]
        r[1] + p[1] + phat[1]
    else
        -Inf
    end

    V[1,1,1] = if r[1] + p[1] + phat[1] <= d[1]
        -Inf
    else
        r[1] + p[1] + phat[1]
    end

    for i in 2:n
        V[i,0,0] = if abs(V[i-1,0,0]) + p[i] <= d[i]
             max(r[i], V[i-1,0,0] + p[i])
        else
            -Inf
        end 
    end

    for i in 2:n
        cum_arr = [0;cumsum(w[1:i])]
        for (α_before, α) in zip(cum_arr[1:end-1], cum_arr[2:end])
            V[i,0,α] = if V[i-1,0,α_before] + p[i] > d[i]
                V[i-1,0,α_before] + p[i]
            elseif abs(V[i-1,0,α]) + p[i] <= d[i]
                max(r[i], V[i-1,0,α] + p[i])
            else
                -Inf
            end
        end
    end

    for i in 2:n
        for γ in 1:i
            V[i,γ,0] = if abs(V[i-1,γ,0]) + p[i] <= d[i] && abs(V[i-1,γ-1,0]) + p[i] + phat[i] <= d[i]
                max(r[i], V[i-1,γ,0] + p[i], V[i-1,γ-1,0] + p[i] + phat[i])
            elseif abs(V[i-1,γ-1,0]) + p[i] + phat[i] <= d[i]
                max(r[i], V[i-1,γ-1,0] + p[i] + phat[i])
            elseif abs(V[i-1,γ,0]) + p[i] <= d[i]
                max(r[i], V[i-1,γ,0] + p[i])
            else
                -Inf
            end
        end
    end

    for i in 2:n
        for γ in 1:i
            cum_arr = [0;cumsum(w[1:i])]
            for (α_before, α) in zip(cum_arr[1:end-1], cum_arr[2:end])
                V[i,γ,α] = if V[i-1,γ-1, α_before] + p[i] + phat[i] > d[i] || V[i-1,γ,α_before] + p[i] > d[i]
                    max(V[i-1,γ,α_before] + p[i], V[i-1,γ-1,α_before] + p[i] + phat[i])
                elseif abs(V[i-1,γ-1,α]) + p[i] + phat[i] <= d[i] && abs(V[i-1,γ,α]) + p[i] <= d[i]
                    max(r[i], V[i-1,γ,α] + p[i], V[i-1,γ-1,α] + p[i] + phat[i])
                elseif abs(V[i-1,γ-1,α]) + p[i] + phat[i] <= d[i]
                    max(r[i], V[i-1,γ-1,α] + p[i] + phat[i])
                elseif abs(V[i-1,γ,α]) + p[i] <= d[i]
                    max(r[i], V[i-1,γ,α] + p[i])
                else
                    -Inf
                end
            end
        end
    end


    
end