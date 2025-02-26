
function generateBinaryArrays(n::Int, Γ::Int)
    return ((i in combo ? 1 : 0 for i in 1:n) |> BitVector for combo in combinations(1:n, Γ) )
end

function generate_wagner_model(n::int, p::Matrix{Int}, phat::Matrix{Int}, Γ1::Int, Γ2::Int)

    model = Model(GLPK.Optimizer)

    # select every Γ1 and Γ2 element subset of the n elements
    # and create a binary variable for each subset 

    Γ1_subsets = generateBinaryArrays(n, Γ1) |> collect
    Γ2_subsets = generateBinaryArrays(n, Γ2) |> collect

    #carthesian product of Γ1 and Γ2 subsets
    Λ = [(i, j) for i in Γ1_subsets, j in Γ2_subsets]

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

    



    





    

    

end