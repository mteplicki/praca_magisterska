function column_generation(model::T; ϵ=10^-6, max_iterations=100) where T <: AbstractColumnGenerationModel

    LB = -Inf
    UB = Inf
    v = 1
    Z = value.(model.Z)
    iteration=0

    while (UB - LB)/LB > ϵ
        optimize!(model.model)

        if termination_status(model.model) != MOI.OPTIMAL
            throw(ArgumentError("Model is not optimal"))
        end

        Z = value.(model.Z)

        y_opt = objective_value(model.model)

        LB = max(LB, y_opt)

        value, λ = oracle_subproblem(model, Z)

        UB = min(UB, value)

        if (UB - LB)/LB > ϵ
            update_model!(model, [λ])
        end
        iteration += 1
        if iteration > max_iterations
            @warn("Maximum number of iterations reached")
            break
        end

    end

    return UB, Z

end