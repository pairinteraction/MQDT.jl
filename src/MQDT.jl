module MQDT

using LinearAlgebra
using Roots
using SparseArrays
using CGcoefficient
using Memoize
using DataFrames
using PythonCall

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
    DataBaseArray,
    eigenstates,
    basisarray,
    databasearray,
    matrix_element,
    state_data,
    matrix_data,
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
