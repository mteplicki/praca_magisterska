struct TestModel{T} where T <: AbstractColumnGenerationModel
    model::T
    time::Float64
    result_value::Float64
    result_permutation
    result_stats::ColumnGenerationStats
    kwargs
end

function test_model(model::T, kwargs...) where T <: AbstractColumnGenerationModel
    # test the model with column_generation

    duration = @elapsed begin
        # solve the model
        value, permutation, stats = column_generation(model, kwargs...)
    end

    return TestModel(model, duration, value, permutation, stats, kwargs)
end

function test_instances(dir_instances::String)
    
end