"""
Module for saving and loading instances of AbstractInstance to and from disk.
"""
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

"""
Loads an instance of type T from a file.
- `path`: Path to the file containing the instance data.
- `T`: Type of the instance to be loaded (must be a subtype of `AbstractInstance`).
"""
function load_instance(path::AbstractString, ::Type{T})::T where T <: AbstractInstance
    open(path) do io
        return deserialize(io)
    end
end

"""
Loads an instance from a file, inferring the type automatically.
- `path`: Path to the file containing the instance data.
"""
function load_instance(path::AbstractString)
    open(path) do io
        return deserialize(io)
    end
end