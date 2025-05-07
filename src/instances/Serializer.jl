function save_instance(path::AbstractString, inst::T) where T <: AbstractInstance
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
        serialize(io, inst)
    end
end

function load_instance(path::AbstractString, ::Typeof{T})::T where T <: AbstractInstance
    open(path) do io
        return deserialize(io)
    end
end
