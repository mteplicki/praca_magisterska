function column_generation(model::T; ϵ=10^-6, max_iterations=100, kwargs...) where T <: AbstractColumnGenerationModel

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

        @show LB
        @show UB
        @show λ

        if (UB - LB) > ϵ
            update_model!(model, [λ])
        end
        iteration += 1
        if iteration > max_iterations
            @warn("Maximum number of iterations reached")
            break
        end
    end

    return UB

end