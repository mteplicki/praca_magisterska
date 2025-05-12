struct TestModel{T} where T <: AbstractInstance
    name::String
    instance::T
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

    return TestModel(string(typeof(model)),model.instance, duration, value, permutation, stats, kwargs)
end

function test_instances(dir_instances::String, optimizer)
    # Get all files in the directory
    files = readdir(dir_instances)
    # for every file in the directory

    #find file with name "model_type"
    model_type = findfirst(x -> occursin("model_type", x), files)
    if model_type == nothing
        error("No file with name 'model_type' found in directory $dir_instances")
    end
    # Read the file
    model_type = read(joinpath(dir_instances, model_type), String)

    test_model_list = []

    for file in files
        # Check if the file is a .jl file
        # Load the instance
        
        #skip if file is a "model_type" 
        if file == "model_type"
            continue
        end
        path = joinpath(dir_instances, file)
        instance = load_instance(path)
        # Create the model
        model = @match model_type begin
            "MultiUnrelatedMakespan" => MultiUnrelatedMakespan(optimizer, instance)
            "SingleTardyJobs" => SingleTardyJobsModel(optimizer, instance)
            "SingleTardinessDominanceRules" => SingleTardinessDominanceRules(optimizer, instance)
            "SingleMachineDueDates" => SingleMachineDueDates(optimizer, instance)
            "2RPFP" => WagnerModel(optimizer, instance)
        end
        # Test the model
        res_test_model = test_model(model)
        # Save the result
        push!(test_model_list, res_test_model)
    end
    # serialize the result
    path = joinpath(dir_instances, "test_results")
    # Check if the directory exists, if not create it
    dir_path = dirname(path)
    if !isdir(dir_path)
        mkpath(dir_path)
    end
    # Check if the file exists, if so, delete it
    if isfile(path)
        rm(path)
    end
    # Open the file for writing
    open(path, "w") do io
        serialize(io, test_model_list)
    end
    # return the result
    return test_model_list
end