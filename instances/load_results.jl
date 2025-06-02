using MasterThesis
using Serialization
using DataFrames

# load results from dir and save them in a load_instance
function load_results(dir::AbstractString)
    results = []
    for file in readdir(dir)
        if endswith(file, ".bin")
            path = joinpath(dir, file)
            try
                result = load_instance(path)
                push!(results, result)
            catch e
                @warn "Failed to load instance from $path: $e"
            end
        end
    end
    return results
end

function SingleTardinessSumInstanceResultsToDataFrame(results::Vector)
    df = DataFrame([
        "name" => String[],
        "n" => Int[],
        "R" => Float64[],
        "T" => Float64[],
        "G" => Float64[],
        "time" => Float64[],
        "res" => Float64[],
        "optimality" => Bool[],
        "iterations" => Int[],
        "LB" => String[],
        "UB" => String[],
        "first_LB" => Float64[],
        "first_UB" => Float64[],
        "last_LB" => Float64[],
        "last_UB" => Float64[]
    ])
    for result in results
        row = Dict(
            "name" => result.name,
            "n" => result.instance.n,
            "R" => result.instance.R,
            "T" => result.instance.T,
            "G" => result.instance.G,
            "time" => result.time,
            "res" => result.result_stats.res,
            "optimality" => result.result_stats.optimality,
            "iterations" => result.result_stats.iterations,
            "LB" => string(result.result_stats.LB),
            "UB" => string(result.result_stats.UB),
            "first_LB" => if isempty(result.result_stats.LB)
                0
            else
                first(result.result_stats.LB)
            end,
            "first_UB" => if isempty(result.result_stats.UB)
                0
            else
                first(result.result_stats.UB)
            end,
            "last_LB" => if isempty(result.result_stats.LB)
                Inf
            else
                last(result.result_stats.LB)
            end,
            "last_UB" => if isempty(result.result_stats.UB)
                Inf
            else
                last(result.result_stats.UB)
            end,
        )
        @show row
        push!(df, row)
    end
    return df
end

function MultiUnrelatedInstanceResultsToDataFrame(results::Vector)
    df = DataFrame([
        "name" => String[],
        "n" => Int[],
        "m" => Int[],
        "G" => Float64[],
        "time" => Float64[],
        "res" => Float64[],
        "optimality" => Bool[],
        "iterations" => Int[],
        "LB" => String[],
        "UB" => String[]
    ])
    for result in results
        row = Dict(
            "name" => result.name,
            "n" => result.instance.n,
            "m" => result.instance.m,
            "G" => result.instance.G,
            "time" => result.time,
            "res" => result.result_stats.res,
            "optimality" => result.result_stats.optimality,
            "iterations" => result.result_stats.iterations,
            "LB" => string(result.result_stats.LB),
            "UB" => string(result.result_stats.UB),
        )

        @show row
        push!(df, row)
    end
    return df
end

function load_MultiUnrelatedMakespanResults(dir::AbstractString)
    a1 = load_results("./instances/unrela/test_results_dir")
    a2 = load_results("./instances/unrela_2/test_results_dir")
    a3 = load_results("./instances/unrela_3/test_results_dir")

    append!(a1, a2)

    append!(a1, a3)
    df = MultiUnrelatedInstanceResultsToDataFrame(a1)

    # remove rows with m=n*0.4
    df = filter(row -> !isapprox(row.m,row.n * 0.4), df)
    df
end
