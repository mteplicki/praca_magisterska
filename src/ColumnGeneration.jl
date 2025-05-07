struct ColumnGenerationStats
    optimality::Ref{Bool}
    LB::Vector{Float64}
    UB::Vector{Float64}
    iterations::Ref{Int}
    optimizer::String
    kwargs
end

function column_generation(model::T; ϵ=10^-6, max_iterations=-1, kwargs...) where T <: AbstractColumnGenerationModel

    # Pobieramy "prawdziwy" optymalizator spod warstw mostów
    optimizer = backend(model.model).optimizer
    while optimizer isa MathOptInterface.AbstractOptimizer && hasproperty(optimizer, :model)
        optimizer = optimizer.model
    end

    solver_name = string(typeof(optimizer))

    stats = ColumnGenerationStats(Ref(false),Float64[], Float64[], Ref(0), solver_name, kwargs)

    LB = -Inf
    UB = Inf
    v = 1
    iteration=0

    while (UB - LB) > ϵ
        @info "Iteration $iteration"

        @show model.model

        optimize!(model.model)

        if termination_status(model.model) != MOI.OPTIMAL
            throw(ArgumentError("Model is not optimal"))
        end

        y_opt = objective_value(model.model)

        LB = max(LB, y_opt)

        value, λ = oracle_subproblem(model, kwargs)

        UB = min(UB, value)

        push!(stats.LB, LB)
        push!(stats.UB, UB)
        

        @show LB
        @show UB
        @show λ

        if (UB - LB) > ϵ
            update_model!(model, [λ], LB)
        end
        iteration += 1
        if iteration > max_iterations && max_iterations != -1
            @warn("Maximum number of iterations reached")
            break
        end
    end
    stats.iterations[] = iteration
    stats.optimality[] = (UB - LB) <= ϵ

    permutation = find_permutation(model)

    return UB, permutation, stats

end