module MQDT

using LinearAlgebra
using Roots
using SparseArrays
using CGcoefficient
using DataFrames
using PythonCall
using LRUCache

export lsQuantumNumbers,
    jjQuantumNumbers,
    fjQuantumNumbers,
    lsChannels,
    jjChannels,
    fjChannels,
    fModel,
    kModel,
    Parameters,
    EigenStates,
    BasisArray,
    eigenstates,
    basisarray,
    matrix_elements,
    single_channel_models

include("general.jl")
include("boundstates.jl")
include("matrixelements.jl")
include("output.jl")

# parameters as submodules
include("parameters/Sr87.jl")
include("parameters/Sr88.jl")
include("parameters/Yb171.jl")
include("parameters/Yb173.jl")
include("parameters/Yb174.jl")

end
