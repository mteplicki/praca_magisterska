struct MultiUnrelatedMakespanVariables <: AbstractVariableRef
    Z::Matrix{VariableRef}
    T::VariableRef
end

struct MultiUnrelatedMakespan <: AbstractColumnGenerationModel
    instance::MultiUnrelatedInstance
    n::Int
    m::Int
    p::Matrix{Int}
    phat::Matrix{Int}
    Λ::Vector{BitVector}
    Γ::Int
    variables::MultiUnrelatedMakespanVariables
    model::Model
end

MultiUnrelatedMakespan(optimizer, instance::MultiUnrelatedInstance) = MultiUnrelatedMakespan(optimizer, instance.n, instance.m, instance.p, instance.phat, instance.Γ, instance.p_min, instance.p_max, instance.G)

function MultiUnrelatedMakespan(optimizer, n::Int, m::Int, p::Matrix{Int}, phat::Matrix{Int}, Γ::Int, p_min, p_max, G::Float64)
    model = Model(optimizer)

    #Lenstra–Shmoys–Tardos linear program
    Λ = [falses(n)]

    @variable(model, Z[1:m, 1:n], Bin)

    @variable(model, T >= 0)

    @constraint(model, [j in 1:n], sum(Z[i,j] for i in 1:m) == 1)

    @constraint(model, [i in 1:m, δ in Λ], sum((p[i,j] + phat[i,j] * δ[j]) * Z[i,j] for j in 1:n) <= T)

    @objective(model, Min, T)

    return MultiUnrelatedMakespan(MultiUnrelatedInstance(n,m,p,phat,Γ, p_min, p_max,G),n, m, p, phat, Λ, Γ, MultiUnrelatedMakespanVariables(Z,T) ,model)
end

function update_model!(model::MultiUnrelatedMakespan, new_Λ::Vector{BitVector}, LB::Float64)
    Z = model.variables.Z
    T = model.variables.T

    @constraint(model.model, [i in 1:model.m, δ in new_Λ], sum((model.p[i,j] + model.phat[i,j] * δ[j]) * Z[i,j] for j in 1:model.n) <= T)
    
    append!(model.Λ, new_Λ)

    # @constraint(model.model, T >= LB)

    return model
end

function check_worst_case(model::MultiUnrelatedMakespan, jobs_to_machines, λ::BitVector, expected_value)
    times = [sum(model.p[i,j] + model.phat[i,j] * λ[j] for j in jobs_to_machines[i]; init=0) for i in 1:model.m]
    result = maximum(times)
    if result != expected_value
        @warn "essa: res $expected_value != $result"
    end
end

function oracle_subproblem(model::MultiUnrelatedMakespan, permutation, kwargs)
    # oracle_subproblem is a upper bound for the problem of minimizing the makespan of unrelated machines
    # given variable Z[i,j] (Z[i,j]=1 if job j is assigned to machine i), it calculates a assigment to machines and returns the value of the objective for all uncertainty scenarios

    # Z is a matrix of size m x n
    # Z[i,j] = 1 if job j is assigned to machine i
    # Z[i,j] = 0 otherwise

    jobs_to_machines = permutation

    nominal_times = [sum(model.p[i,j] for j in jobs_to_machines[i]; init=0) for i in 1:model.m]

    overtimes = [[(model.phat[i,j], j) for j in jobs_to_machines[i]] for i in 1:model.m]

    for i in 1:model.m
        sort!(overtimes[i], by = x -> x[1], rev = true)
        overtimes[i] = overtimes[i][1:min(model.Γ,length(overtimes[i]))]
    end

    sum_of_overtimes = [sum(x[1] for x in overtimes[i]; init=0) for i in 1:model.m]

    times = nominal_times .+ sum_of_overtimes

    # select the machine with the largest time
    machine_worst = argmax(times)

    value_worst = times[machine_worst]

    jobs_deviate = [x[2] for x in overtimes[machine_worst]]

    if length(jobs_deviate) < model.Γ
        times_machines = [(times[machine], machine) for machine in 1:model.m]
        sort!(times_machines, by = x -> x[1], rev = true)
        machines_ordered = [x[2] for x in times_machines[2:end]]

        for machine in machines_ordered
            jobs_deviate_other_machine = [x[2] for x in overtimes[machine]]
            if length(jobs_deviate) + length(jobs_deviate_other_machine) >= model.Γ
                append!(jobs_deviate, jobs_deviate_other_machine[1:(model.Γ - length(jobs_deviate))])
                break
            else
                append!(jobs_deviate, jobs_deviate_other_machine)
            end
        end
    end

    λ = BitVector([j in jobs_deviate for j in 1:model.n])

    check_worst_case(model, jobs_to_machines, λ, value_worst)
    
    return value_worst, λ
end

function oracle_subproblem2(model::MultiUnrelatedMakespan, jobs_to_machines, kwargs)
    # oracle_subproblem is a upper bound for the problem of minimizing the makespan of unrelated machines
    # given variable Z[i,j] (Z[i,j]=1 if job j is assigned to machine i), it calculates a assigment to machines and returns the value of the objective for all uncertainty scenarios

    # Z is a matrix of size m x n
    # Z[i,j] = 1 if job j is assigned to machine i
    # Z[i,j] = 0 otherwise


    nominal_times = [sum(model.p[i,j] for j in jobs_to_machines[i]; init=0) for i in 1:model.m]

    overtimes = [[(model.phat[i,j], j) for j in jobs_to_machines[i]] for i in 1:model.m]

    for i in 1:model.m
        sort!(overtimes[i], by = x -> x[1], rev = true)
        overtimes[i] = overtimes[i][1:min(model.Γ,length(overtimes[i]))]
    end

    sum_of_overtimes = [sum(x[1] for x in overtimes[i]; init=0) for i in 1:model.m]

    times = nominal_times .+ sum_of_overtimes

    # select the machine with the largest time
    machine_worst = argmax(times)

    @show machine_worst

    value_worst = times[machine_worst]

    jobs_deviate = [x[2] for x in overtimes[machine_worst]]

    if length(jobs_deviate) < model.Γ
        times_machines = [(times[machine], machine) for machine in 1:model.m]
        sort!(times_machines, by = x -> x[1], rev = true)
        machines_ordered = [x[2] for x in times_machines[2:end]]

        for machine in machines_ordered
            jobs_deviate_other_machine = [x[2] for x in overtimes[machine]]
            if length(jobs_deviate) + length(jobs_deviate_other_machine) >= model.Γ
                append!(jobs_deviate, jobs_deviate_other_machine[1:(model.Γ - length(jobs_deviate))])
                break
            else
                append!(jobs_deviate, jobs_deviate_other_machine)
            end
        end
    end

    λ = BitVector([j in jobs_deviate for j in 1:model.n])

    expected_value = value_worst    

    times1 = [sum(model.p[i,j] + model.phat[i,j] * λ[j] for j in jobs_to_machines[i]; init=0) for i in 1:model.m]
    result1 = maximum(times1)
    if result1 != expected_value
        @warn "essa: res $expected_value != $result"
    end
    
    return value_worst, λ
end

function find_permutation(model::MultiUnrelatedMakespan)
    # Z is a matrix of size m x n
    # Z[i,j] = 1 if job j is assigned to machine i
    # Z[i,j] = 0 otherwise

    Z = value.(model.variables.Z)

    Z = Z .> 0.5

    jobs_to_machines = [[] for _ in 1:model.m]

    for j in 1:model.n
        for i in 1:model.m
            if Z[i,j] == 1
                push!(jobs_to_machines[i], j)
            end
        end
    end

    return jobs_to_machines
end

save_permutation_variable(model::MultiUnrelatedMakespan) = value.(model.variables.Z)

set_start_value_for_model(model::MultiUnrelatedMakespan, permutation_variable) = begin
    if !isnothing(permutation_variable) 
        set_start_value.(model.variables.Z, permutation_variable)
    end
end
