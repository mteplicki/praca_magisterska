mutable struct WagnerModelVariables <: AbstractVariableRef
    Z::Matrix{VariableRef}
    X::Matrix{VariableRef}
    Y::Matrix{VariableRef}
    y::VariableRef
end

struct WagnerModel <: AbstractColumnGenerationModel
    instance::TwoRPFPInstance
    n::Int
    p::Matrix{Float64}
    phat::Matrix{Float64}
    Λ::Vector{Tuple{BitVector, BitVector}}
    variables::WagnerModelVariables
    Γ1::Int
    Γ2::Int
    model::Model
end

"""
Constructor that creates a Wagner MIP program from an instance of `TwoRPFPInstance`.
"""
WagnerModel(optimizer, instance::TwoRPFPInstance) = WagnerModel(optimizer, instance.n, instance.p, instance.phat, (instance.Γ1, instance.Γ2))


function WagnerModel(optimizer, n::Int, p::Matrix{Float64}, phat::Matrix{Float64}, Γ::Tuple{Int, Int})
    model = direct_model(optimizer)

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

    for (λ,(δ1, δ2)) in enumerate(Λ)
        @constraint(model, Y[1, λ] == 0)
    end

    return WagnerModel(TwoRPFPInstance(n,p,phat,Γ[1],Γ[2]),n, p, phat, Λ, variables, Γ1, Γ2, model)
end


function update_model!(model::WagnerModel, new_Λ::Vector{Tuple{BitVector, BitVector}}, LB::Float64)
    X = @variable(model.model, [1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))], lower_bound = 0, base_name = "X")
    Y = @variable(model.model, [1:model.n, (length(model.Λ)+1):(length(model.Λ)+length(new_Λ))], lower_bound = 0, base_name = "Y")

    model.variables.X = hcat(model.variables.X, X)
    model.variables.Y = hcat(model.variables.Y, Y)

    # @constraint(model.model, model.variables.y >= LB)

    # update constaints (5)-(7)
    
    λ = length(model.Λ) + 1
    for (δ1, δ2) in new_Λ
        @constraint(model.model, sum(model.p[2,i]+ model.phat[2,i] * δ2[i] for i in 1:model.n) + sum(model.variables.X[j, λ] for j in 1:model.n) <= model.variables.y)
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

    λ = length(model.Λ) + 1
    for (δ1, δ2) in new_Λ
        @constraint(model.model, model.variables.Y[1, λ] == 0)
        λ += 1
    end

    # model.Λ = vcat(model.Λ, new_Λ)
    append!(model.Λ, new_Λ)

    return model
end

"""
Solves the oracle subproblem for a given job permutation `permutation`.
"""
function oracle_subproblem(model::WagnerModel, permutation, kwargs)
    # Z is a permutation matrix. Transform it into a permutation vector
    
    σ = permutation

    # solve subproblem using dynamic programming - top-down approach
    # initialize the dynamic programming table
    α_offset = Array{Union{Missing, Float64}}(missing,2, model.n+1, model.n+2)


    α = OffsetArray(α_offset, 1:2, 0:model.n, -1:model.n)

    α[1:2, 0:model.n, -1] .= typemin(Float64)

    α[1:2, 0, 0:model.n] .= 0

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

    # @show "essa"
    for k in 1:model.n
        α[2,k,0] = model.p[2,σ[k]] + max(α[2,k-1,0], α[1,k,model.Γ1])
        if α[2,k-1,0] > α[1,k,model.Γ1]
            delay_array[2,k,0] = (0,0)
        else
            delay_array[2,k,0] = (0,1)
        end
        # @show k
        # @show α[2,k,0]
    end

    for k in 1:model.n
        for γ in 1:model.Γ2
            # @show k
            # @show γ
            # @show α[2,k-1,γ]
            # @show α[1,k,model.Γ1]
            # @show α[2,k-1,γ-1]
            # @show α[1,k,model.Γ1]
            no_delay = model.p[2,σ[k]] + max(α[2,k-1,γ], α[1,k,model.Γ1])
            delay = model.p[2,σ[k]] + model.phat[2,σ[k]] + max(α[2,k-1,γ-1], α[1,k,model.Γ1])
            # @show delay
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

    # @show delay_array

    value_α = α[2,model.n,model.Γ2]

    # @show value_α

    λ1 = falses(model.n)
    λ2 = falses(model.n)

    γ = model.Γ2

    starting_point = -1

    for k in model.n:-1:1
        # @show k
        # @show σ[k]
        # @show delay_array[2,k,γ]
        # @show γ
        # @show starting_point
        if delay_array[2,k,γ] == (1,0)
            λ2[σ[k]] = 1
            γ -= 1
        elseif delay_array[2,k,γ] == (1,1)
            λ2[σ[k]] = 1
            if starting_point == -1
                starting_point = k
            end
            γ -= 1
        elseif delay_array[2,k,γ] == (0,1) && starting_point == -1
            starting_point = k
        end
    end

    γ = model.Γ1
    for k in (starting_point):-1:1
        # @show k
        # @show σ[k]
        # @show delay_array[1,k,γ]
        # @show γ
        if delay_array[1,k,γ] == (1,0)
            λ1[σ[k]] = 1
            γ -= 1
        end
    end

    # @show λ1
    # @show λ2

    # check if the solution λ1, λ2 is feasible
    machine_1_time = 0.
    machine_2_time = 0.
    for k in 1:model.n
        # @show k
        # @show λ1[σ[k]]
        if λ1[σ[k]] == 1
            machine_1_time += model.p[1,σ[k]] + model.phat[1,σ[k]]
        else
            machine_1_time += model.p[1,σ[k]]
        end
        # @show machine_1_time
        # @show λ2[σ[k]]
        if λ2[σ[k]] == 1
            machine_2_time = max(machine_2_time, machine_1_time) + model.p[2,σ[k]] + model.phat[2,σ[k]]
        else
            machine_2_time = max(machine_2_time, machine_1_time) + model.p[2,σ[k]]
        end
        # @show machine_2_time
    end

    if !(machine_2_time ≈ value_α)
        @error("Solution is not feasible")
        @show machine_2_time
        @show value_α
    end

    # Y = value.(model.variables.Y)
    # X = value.(model.variables.X)

    # if sum(λ1) < model.Γ1
    #     # select all the jobs that are not delayed, and sort them by Y value

    #     number_of_jobs_to_delay = model.Γ1 - sum(λ1)

    #     jobs_not_delayed = [i for i in 1:model.n if λ1[i] == 0]

    #     Y_sum = [sum(Y[i, j] for j in 1:length(model.Λ)) for i in 1:model.n]

    #     # select the number_of_jobs_to_delay jobs with the highest Y value
    #     jobs_not_delayed = sort(jobs_not_delayed, by = i -> Y_sum[i], rev = true)[1:number_of_jobs_to_delay]

    #     for i in jobs_not_delayed
    #         λ1[i] = 1
    #     end
    # end

    # if sum(λ2) < model.Γ2
    #     # select all the jobs that are not delayed, and sort them by X value 

    #     number_of_jobs_to_delay = model.Γ2 - sum(λ2)

    #     jobs_not_delayed = [i for i in 1:model.n if λ2[i] == 0]

    #     #sum the X by column
    #     sum_X = [sum(X[i, j] for j in 1:length(model.Λ)) for i in 1:model.n]

    #     # select the number_of_jobs_to_delay jobs with the highest Y value

    #     jobs_not_delayed = sort(jobs_not_delayed, by = i -> sum_X[i])[1:number_of_jobs_to_delay]

    #     for i in jobs_not_delayed
    #         λ2[i] = 1
    #     end
    # end

    return value_α, (λ1, λ2)
end

"""
Extracts the job permutation from the solution of the Wagner model.
"""
function find_permutation(model::WagnerModel)
    Z_float = value.(model.variables.Z)

    #convert Z to a matrix of 0s and 1s
    Z = Z_float .> 0.5
    
    σ = [findfirst(Z[:,i]) for i in 1:model.n]

    return σ
end

"""
Saves the current permutation variable values from the model.
"""
save_permutation_variable(model::WagnerModel) = value.(model.variables.Z)

"""
Sets the start values for the permutation variable Z based on a given permutation.
"""
function set_start_value_for_model(model::WagnerModel, permutation)
    if !isnothing(permutation)
        set_start_value.(model.variables.Z, permutation)
    end
end