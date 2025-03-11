function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end


function solve_subproblem_with_feasibility(subproblem::Model, x)
    fix.(subproblem[:x_copy], x)
    optimize!(model)
    if is_solved_and_feasible(subproblem; dual = true)
        return (
            is_feasible = true,
            obj = objective_value(subproblem),
            y = value.(subproblem[:y]),
            π = reduced_cost.(subproblem[:x_copy]),
        )
    end
    return (
        is_feasible = false,
        v = dual_objective_value(subproblem),
        u = reduced_cost.(subproblem[:x_copy]),
    )
end

function benders_decomposition(model::T; ϵ=10^-6, max_iterations=100, verbose=false ) where T <: AbstractBendersDecompositionModel
    verbose && println("Iteration  Lower Bound  Upper Bound          Gap")
    for k in 1:max_iterations
        optimize!(model.model)
        assert_is_solved_and_feasible(model.model)
        lower_bound = objective_value(model.model)
        x_k = value.(x)
        ret = solve_subproblem_with_feasibility(model, x_k)
        if ret.is_feasible
            # Benders Optimality Cuts
            upper_bound = (objective_value(model.model) - value(θ)) + ret.obj
            gap = abs(upper_bound - lower_bound) / abs(upper_bound)
            verbose && print_iteration(k, lower_bound, upper_bound, gap)
            if gap < ϵ
                verbose && println("Terminating with the optimal solution")
                break
            end
            @constraint(model.model, θ >= ret.obj + sum(ret.π .* (x .- x_k)))
        else
            # Benders Feasibility Cuts
            cut = @constraint(model.model, ret.v + sum(ret.u .* (x .- x_k)) <= 0)
            verbose && @info "Adding the feasibility cut $(cut)"
        end
    end

end