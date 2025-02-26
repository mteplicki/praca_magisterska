
function generateBinaryArrays(n::Int, Γ::Int)
    return ((i in combo ? 1 : 0 for i in 1:n) |> BitVector for combo in combinations(1:n, Γ) )
end

struct WagnerModel
    n::Int
    p::Matrix{Int}
    phat::Matrix{Int}
    Λ::Vector{Tuple{BitVector, BitVector}}
    model::Model
end

function generate_wagner_model(n::int, p::Matrix{Int}, phat::Matrix{Int}, Λ::Vector{Tuple{BitVector, BitVector}})
    model = Model(GLPK.Optimizer)

    @variable(model, Z[1:n, 1:n], Bin) 

    @variable(model, X[1:n, 1:length(Γ_subsets)] >= 0, Bin)

    @variable(model, Y[1:n, 1:length(Γ_subsets)] >= 0, Bin) 

    @variable(model, y >= 0)

    # (4)
    @objective(model, Min, y)

    # (5)
    λ = 1
    for (λ,(δ1, δ2)) in enumerate(Λ)
        @constraint(model, sum(p[2,i]+ phat[2,i] * δ2[i] for i in 1:n) + sum(X[j, λ] for j in 1:n) <= y)
        λ += 1
    end

    # (6)
    for (λ,(δ1, δ2)) in enumerate(Λ)
        @constraint(model, [j in 1:(n-1)],
            sum((p[1,i]+ phat[1,i] * δ1[i])*Z[i,j+1] for i in 1:n) + Y[j+1, λ] == 
            sum((p[2,i]+ phat[2,i] * δ2[i])*Z[i,j] for i in 1:n) + X[j+1, λ] + Y[j, λ])
    end

    # (7)
    for (λ,(δ1, δ2)) in enumerate(Λ)
        @constraint(model, sum((p[1,i]+ phat[1,i] * δ1[i])*Z[i,1] for i in 1:n) = X[1, λ])
    end

    # (8)
    @constraint(model, [j in 1:n], sum(Z[i,j] for i in 1:n) == 1)

    # (9)
    @constraint(model, [i in 1:n], sum(Z[i,j] for j in 1:n) == 1)


    return WagnerModel(n, p, phat, Λ, model)
end


function update_wagner_model(model::WagnerModel, new_Λ::Vector{Tuple{BitVector, BitVector}})
    @variable(model.model, X[1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))] >= 0, Bin)
    @variable(model.model, Y[1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))] >= 0, Bin)

    # update constaints (5)-(7)
    
    λ = length(model.Λ) + 1
    for (δ1, δ2) in new_Λ
        @constraint(model.model, sum(model.p[2,i]+ model.phat[2,i] * δ2[i] for i in 1:model.n) + sum(X[j, λ] for j in 1:model.n) <= model.y)
        λ += 1
    end

    λ = length(model.Λ) + 1
    for (δ1, δ2) in new_Λ
        @constraint(model.model, [j in 1:(model.n-1)],
            sum((model.p[1,i]+ model.phat[1,i] * δ1[i])*model.Z[i,j+1] for i in 1:model.n) + model.Y[j+1, λ] == 
            sum((model.p[2,i]+ model.phat[2,i] * δ2[i])*model.Z[i,j] for i in 1:model.n) + X[j+1, λ] + Y[j, λ])
        λ += 1
    end

    λ = length(model.Λ) + 1
    for (δ1, δ2) in new_Λ
        @constraint(model.model, sum((model.p[1,i]+ model.phat[1,i] * δ1[i])*model.Z[i,1] for i in 1:model.n) = X[1, λ])
        λ += 1
    end

    model.Λ = vcat(model.Λ, new_Λ)

    return model
end


function oracle_subproblem(model::WagnerModel, Z::Matrix{Int})
    # Z is a permutation matrix. Transform it into a permutation vector
    σ = [findfirst(Z[i,:]) for i in 1:model.n]

    # solve subproblem using dynamic programming - top-down approach
    # initialize the dynamic programming table
    α = fill(typemin(Int), (2, model.n, model.n+1))

    α[1:2, 0, :] .= 0

    for k in 1:model.n
        α[2,k,1] = model.p[2,σ[k]] + max()





    
end