struct TestModel{T<:AbstractInstance}
    name::String
    instance::T
    time::Float64
    result_stats::ColumnGenerationStats
    kwargs
end

function test_model(model::T, kwargs...) where T <: AbstractColumnGenerationModel
    # test the model with column_generation

    duration = @elapsed begin
        # solve the model
        stats = column_generation(model; timeout=7200, kwargs...)
    end

    return TestModel(string(typeof(model)),model.instance, duration, stats, kwargs)
end

function test_instances(dir_instances::String, optimizer)
    # Get all files in the directory
    files = readdir(dir_instances)
    # for every file in the directory

    #find file with name "model_type"
    model_type_index = findfirst(x -> occursin("model_type", x), files)
    if model_type_index == nothing
        error("No file with name 'model_type' found in directory $dir_instances")
    end
    model_type = files[model_type_index]
    # Read the file
    model_type = read(joinpath(dir_instances, model_type), String)

    # create folder for test results inside the directory
    # Check if the directory exists, if not create it

    test_dir_path = joinpath(dir_instances, "test_results_dir")
    if !isdir(test_dir_path)
        mkpath(test_dir_path)
    end
    println("Test results will be saved in $test_dir_path")

    test_model_list = []

    error_files = []
    
    for file in files
        try
            # Check if the file is a .jl file
            # Load the instance
        
            #skip if file is a "model_type" 
            if file == "model_type"
                continue
            end

            #check if file is a directory
            if isdir(joinpath(dir_instances, file))
                continue
            end

            path = joinpath(dir_instances, file)
            instance = load_instance(path)
            @show "Loaded instance from file: $path"
            # @show instance
            # Create the model
            model = @match model_type begin
                "MultiUnrelatedMakespan" => MultiUnrelatedMakespan(optimizer, instance)
                "SingleTardyJobs" => SingleTardyJobsModel(optimizer, instance)
                "SingleTardinessDominanceRules" => SingleTardinessDominanceRules(optimizer, instance)
                "SingleMachineDueDates" => SingleMachineDueDates(optimizer, instance)
                "2RPFP" => WagnerModel(optimizer, instance)
            end
            # @show model
            # Test the model
            res_test_model = test_model(model)
            # Save the result
            push!(test_model_list, res_test_model)

            # Save the result to a file inside the test_results_dir
            # Check if the directory exists, if not create it
            path = joinpath(test_dir_path, file)
            # Check if the file exists, if so, delete it
            if isfile(path)
                rm(path)
            end
            # Open the file for writing
            open(path, "w") do io
                serialize(io, res_test_model)
            end
    
        catch e
            @error "Error in file $file: $(e)"
            # Add the file to the error list
            push!(error_files, file)
            # continue to the next file
            continue
        end     
    end

    @show test_model_list

    @show error_files

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
    return test_model_list, error_files
end