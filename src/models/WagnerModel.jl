struct WagnerModelVariables <: AbstractVariableRef
    Z::Matrix{VariableRef}
    X::Matrix{VariableRef}
    Y::Matrix{VariableRef}
    y::VariableRef
end

struct WagnerModel <: AbstractColumnGenerationModel
    n::Int
    p::Matrix{Int}
    phat::Matrix{Int}
    Λ::Vector{Tuple{BitVector, BitVector}}
    variables::WagnerModelVariables
    Γ1::Int
    Γ2::Int
    model::Model
end

function WagnerModel(n::Int, p::Matrix{Int}, phat::Matrix{Int}, Γ::Tuple{Int, Int})
    model = Model(GLPK.Optimizer)

    Λ = [(falses(n), falses(n))]

    Γ1::Int, Γ2::Int = Γ

    @variable(model, Z[1:n, 1:n], Bin) 

    Z = Z

    @variable(model, X[1:n, 1:length(Λ)] >= 0)

    @variable(model, Y[1:n, 1:length(Λ)] >= 0) 

    @variable(model, y >= 0)

    # (4)
    @objective(model, Min, y)

    variables = WagnerModelVariables(Z, X, Y, y)

    # (5)
    for (λ,(δ1, δ2)) in enumerate(Λ)
        @constraint(model, sum(p[2,i]+ phat[2,i] * δ2[i] for i in 1:n) + sum(X[j, λ] for j in 1:n) <= y)
    end

    # (6)
    for (λ,(δ1, δ2)) in enumerate(Λ)
        @constraint(model, [j in 1:(n-1)],
            sum((p[1,i]+ phat[1,i] * δ1[i])*Z[i,j+1] for i in 1:n) + Y[j+1, λ] == 
            sum((p[2,i]+ phat[2,i] * δ2[i])*Z[i,j] for i in 1:n) + X[j+1, λ] + Y[j, λ])
    end

    # (7)
    for (λ,(δ1, δ2)) in enumerate(Λ)
        @constraint(model, sum((p[1,i]+ phat[1,i] * δ1[i])*Z[i,1] for i in 1:n) == X[1, λ])
    end

    # (8)
    @constraint(model, [j in 1:n], sum(Z[i,j] for i in 1:n) == 1)

    # (9)
    @constraint(model, [i in 1:n], sum(Z[i,j] for j in 1:n) == 1)


    return WagnerModel(n, p, phat, Λ, variables, Γ1, Γ2, model)
end


function update_model!(model::WagnerModel, new_Λ::Vector{Tuple{BitVector, BitVector}})
    X = @variable(model.model, [1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))], lower_bound = 0, base_name = "X")
    Y = @variable(model.model, [1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))], lower_bound = 0, base_name = "Y")

    model.variables.X = hcat(model.variables.X, X)
    model.variables.Y = hcat(model.variables.Y, Y)

    # update constaints (5)-(7)
    
    λ = length(model.Λ) + 1
    for (δ1, δ2) in new_Λ
        @constraint(model.model, sum(model.p[2,i]+ model.phat[2,i] * δ2[i] for i in 1:model.n) + sum(model.variables.X[j, λ] for j in 1:model.n) <= model.y)
        λ += 1
    end

    λ = length(model.Λ) + 1
    for (δ1, δ2) in new_Λ
        @constraint(model.model, [j in 1:(model.n-1)],
            sum((model.p[1,i]+ model.phat[1,i] * δ1[i])*model.variables.Z[i,j+1] for i in 1:model.n) + model.variables.Y[j+1, λ] == 
            sum((model.p[2,i]+ model.phat[2,i] * δ2[i])*model.variables.Z[i,j] for i in 1:model.n) + model.variables.X[j+1, λ] + model.variables.Y[j, λ])
        λ += 1
    end

    λ = length(model.Λ) + 1
    for (δ1, δ2) in new_Λ
        @constraint(model.model, sum((model.p[1,i]+ model.phat[1,i] * δ1[i])*model.variables.Z[i,1] for i in 1:model.n) == model.variables.X[1, λ])
        λ += 1
    end

    model.Λ = vcat(model.Λ, new_Λ)

    return model
end


function oracle_subproblem(model::WagnerModel, kwargs)
    # Z is a permutation matrix. Transform it into a permutation vector

    Z = value.(model.variables.Z)
    
    σ = [findfirst(Z[i,:]) for i in 1:model.n]

    # solve subproblem using dynamic programming - top-down approach
    # initialize the dynamic programming table
    α_offset = Array{Union{Missing, Int}}(missing,2, model.n+1, model.n+2)


    α = OffsetArray(α_offset, 1:2, 0:model.n, -1:model.n)

    α[1:2, 0:model.n, -1] .= typemin(Int)

    α[1:2, 0, 0:model.n] .= 0

    for k in 1:model.n
        α[2,k,0] = model.p[2,σ[k]] + max(α[2,k-1,0], α[1,k,Γ1])
    end


    delay_array = [(0,0) for _ in α]



    # r = 1
    for k in 1:model.n
        for γ in 0:model.Γ1
            no_delay = model.p[1,σ[k]] + α[1,k-1,γ]
            delay = model.p[1,σ[k]] + model.phat[1,σ[k]] + α[1,k-1,γ-1]
            α[1,k,γ] = max(no_delay, delay)
            if delay > no_delay
                delay_array[1,k,γ] = (1,0)
            end
        end
    end

    for k in 1:model.n
        for γ in 0:model.Γ1
            no_delay = model.p[2,σ[k]] + max(α[2,k-1,γ], α[1,k,model.Γ1])
            delay = model.p[2,σ[k]] + model.phat[2,σ[k]] + max(α[2,k-1,γ-1], α[1,k,model.Γ1])
            α[2,k,γ] = max(no_delay, delay)
            if delay > no_delay
                if α[2,k-1,γ-1] > α[1,k,model.Γ1]
                    delay_array[2,k,γ] = (1,0)
                else
                    delay_array[2,k,γ] = (1,1)
                end
            else
                if α[2,k-1,γ] > α[1,k,model.Γ1]
                    delay_array[2,k,γ] = (0,0)
                else
                    delay_array[2,k,γ] = (0,1)
                end
            end
        end
    end

    value = α[2,model.n,model.Γ2]

    λ1 = BitVector(1:model.n)
    λ2 = BitVector(1:model.n)

    γ = model.Γ2

    starting_point = missing

    for k in model.n:-1:1
        if delay_array[2,k,γ] == (1,0)
            λ2[σ[k]] = 1
            γ -= 1
        end
        if delay_array[2,k,model.Γ2] == (1,1)
            λ2[σ[k]] = 1
            if starting_point == missing
                starting_point = k
            end
            γ -= 1
        end
        if delay_array[2,k,γ] == (0,1) && starting_point == missing
            starting_point = k
        end
    end

    γ = model.Γ1
    for k in model.n:-1:1
        if delay_array[1,k,γ] == (1,0)
            λ1[σ[k]] = 1
            γ -= 1
        end
    end

    return value, (λ1, λ2)
end