struct ColumnGenerationStats
    res::Float64
    permutation
    optimality::Bool
    LB::Vector{Float64}
    UB::Vector{Float64}
    iterations::Int
    optimizer::String
    time::Float64
    kwargs
end

function column_generation(model::T; ϵ=10^-6, max_iterations=-1, timeout=7200, kwargs...) where T <: AbstractColumnGenerationModel
    time_start = time()
    # Pobieramy "prawdziwy" optymalizator spod warstw mostów
    optimizer = backend(model.model).optimizer
    while optimizer isa MathOptInterface.AbstractOptimizer && hasproperty(optimizer, :model)
        optimizer = optimizer.model
    end

    solver_name = string(typeof(optimizer))

    LB_vec = Float64[]
    UB_vec = Float64[]

    # stats = ColumnGenerationStats(Ref(false),Float64[], Float64[], Ref(0), solver_name, kwargs)

    LB = -Inf
    UB = Inf
    v = 1
    iteration=0
    permutation = []

    permutation_variable = nothing

    while (UB - LB) > ϵ
        @info "Iteration $iteration"

        println(summary(model.instance))

        @show model.model

        @show length(model.Λ)

        # set timeout for CPLEX 
        if timeout != -1
            current_time = time()
            time_left = timeout - (current_time - time_start)
            if time_left < 0
                @warn("Timeout reached")
                break
            end
            set_optimizer_attribute(model.model, "CPXPARAM_TimeLimit",time_left)
        end

        set_start_value_for_model(model, permutation_variable)

        optimize!(model.model)

        permutation = find_permutation(model)

        permutation_variable = save_permutation_variable(model)

        @show permutation

        if termination_status(model.model) != MOI.OPTIMAL
            @warn("Model is not optimal")
            time_elapsed = time() - time_start
            return ColumnGenerationStats(UB, [], false, LB_vec, UB_vec, iteration, solver_name, time_elapsed, kwargs)
        end

        y_opt = objective_value(model.model)

        LB = max(LB, y_opt)

        value_subproblem, λ = oracle_subproblem(model, permutation, kwargs)

        UB = min(UB, value_subproblem)

        push!(LB_vec, LB)
        push!(UB_vec, UB)
        

        @show LB
        @show UB
        @show λ

        if (UB - LB) > ϵ
            # if iteration == 2
            #     λ = (trues(model.n), trues(model.n))
            # end
            update_model!(model, [λ], LB)
            iteration += 1
            if iteration > max_iterations && max_iterations != -1
                @warn("Maximum number of iterations reached")
                break
            end
            if timeout != -1
                time_elapsed = time() - time_start
                if time_elapsed > timeout
                    @warn("Timeout reached")
                    break
                end
            end
        end
    end

    time_elapsed = time() - time_start

    @show UB
    @show permutation
    @show (UB - LB) <= ϵ
    @show LB_vec
    @show UB_vec
    @show iteration
    @show solver_name
    @show kwargs
    @show time_elapsed

    return ColumnGenerationStats(UB, permutation, (UB - LB) <= ϵ, LB_vec, UB_vec, iteration, solver_name, time_elapsed, kwargs)

end
