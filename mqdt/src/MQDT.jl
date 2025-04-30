module MQDT

using LinearAlgebra
using Roots
using SparseArrays
using GSL
using QuadGK
using Memoize
using DataFrames

export 
    Model, Parameters, EigenStates, BasisArray, DataBaseArray, # structs
    eigenstates, basisarray, databasearray, matrix_element, state_data, matrix_data, tri_to_full # functions

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