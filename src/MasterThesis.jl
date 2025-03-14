module MasterThesis

using DataStructures
using JuMP
using GLPK
using Combinatorics
using OffsetArrays

abstract type AbstractColumnGenerationModel end
abstract type AbstractBendersDecompositionModel end

include("BendersDecomposition.jl")

include("ColumnGeneration.jl")

include("models/WagnerModel.jl")

include("models/WagnerModelBenders.jl")

include("models/SingleTardyJobs.jl")


end