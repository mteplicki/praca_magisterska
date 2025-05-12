module MasterThesis

using DataStructures
using JuMP
using GLPK
using Combinatorics
using OffsetArrays
using Random
using Printf
using Serialization

export AbstractColumnGenerationModel, AbstractBendersDecompositionModel, AbstractVariableRef

export WagnerModel, WagnerModelBenders, SingleTardiness, SingleTardyJobsModel, MultiUnrelatedMakespan

export SingleMachineDueDates, MultiUnrelatedInstance

export SingleTardinessDominanceRules

export column_generation, benders_decomposition

abstract type AbstractColumnGenerationModel end
abstract type AbstractBendersDecompositionModel end

abstract type AbstractInstance end

abstract type AbstractVariableRef end


include("BendersDecomposition.jl")

include("ColumnGeneration.jl")

include("instances/Serializer.jl")

include("instances/SingleMachineDueDates.jl")

include("instances/MultiUnrelatedInstance.jl")

include("models/WagnerModel.jl")

# include("models/WagnerModelBenders.jl")

include("models/SingleTardyJobs.jl")

include("models/MultiUnrelatedMakespan.jl")

include("models/SingleTardySum.jl")

include("models/SingleTardinessDominanceRules.jl")

include("test_scripts/test_script.jl")


end