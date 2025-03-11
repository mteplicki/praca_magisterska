struct WagnerModelBenders <: AbstractBendersDecompositionModel
    model::Model
    subproblem::Model
end

function generateBinaryArrays(n::Int, Γ::Int)
    # Dla każdej kombinacji indeksów o długości γ tworzymy ciąg binarny
    return ([i in combo ? 1 : 0 for i in 1:n] for combo in combinations(1:n, Γ))
end

function generate_all_scenarios(n::int, Γ::Tuple{Int, Int})
    return [(collect(δ1), collect(δ2)) for δ1 in generateBinaryArrays(n, Γ[1]) for δ2 in generateBinaryArrays(n, Γ[2])]
end

function WagnerModelBenders(n::int, p::Matrix{Int}, phat::Matrix{Int}, Γ::Tuple{Int, Int})
        model = Model(GLPK.Optimizer)
        subproblem = Model(GLPK.Optimizer)

        set_attribute(subproblem, "presolve", "off")
    
        model.Λ = generate_all_scenarios(n, Γ)

        M = Inf # M is a large number, ustal konkretne M
    
        @variable(model, Z[1:n, 1:n], Bin)

        @variable(model, Θ >= M)

        @variable(subproblem, Z_sub[1:n, 1:n], Bin)
    
        @variable(subproblem, X[1:n, 1:length(model.Λ)] >= 0)
    
        @variable(subproblem, Y[1:n, 1:length(model.Λ)] >= 0)
    
        @variable(subproblem, y >= 0)

        @objective(subproblem, Min, Θ)
    
        # (4)
        @objective(subproblem, Min, y)
    
        # (5)
        for (λ,(δ1, δ2)) in enumerate(Λ)
            @constraint(subproblem, sum(p[2,i]+ phat[2,i] * δ2[i] for i in 1:n) + sum(X[j, λ] for j in 1:n) <= y)
        end
    
        # (6)
        for (λ,(δ1, δ2)) in enumerate(Λ)
            @constraint(subproblem, [j in 1:(n-1)],
                sum((p[1,i]+ phat[1,i] * δ1[i])*Z[i,j+1] for i in 1:n) + Y[j+1, λ] == 
                sum((p[2,i]+ phat[2,i] * δ2[i])*Z[i,j] for i in 1:n) + X[j+1, λ] + Y[j, λ])
        end
    
        # (7)
        for (λ,(δ1, δ2)) in enumerate(Λ)
            @constraint(subproblem, sum((p[1,i]+ phat[1,i] * δ1[i])*Z[i,1] for i in 1:n) = X[1, λ])
        end
    
        # (8)
        @constraint(model, [j in 1:n], sum(Z[i,j] for i in 1:n) == 1)
    
        # (9)
        @constraint(model, [i in 1:n], sum(Z[i,j] for j in 1:n) == 1)
    
    
        return WagnerModel(n, p, phat, Λ, Γ1, Γ2, model)
    
end