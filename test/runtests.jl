using MQDT
using Test
using LinearAlgebra
using CGcoefficient

wigner_init_float(13, "Jmax", 9) # initialize Wigner symbol calculation

@testset "MQDT.jl" begin
    include("small_test.jl")
    include("test_model.jl")
end
