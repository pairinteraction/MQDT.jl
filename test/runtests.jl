using MQDT
using Test
using GSL
using LinearAlgebra
using CGcoefficient

wigner_init_float(13, "Jmax", 9) # initialize Wigner symbol caluclation

@testset "MQDT.jl" begin
    # TODO replace the following code with proper tests
    s = MQDT.eigenstates(24, 30, MQDT.Yb171.FMODEL_HIGHN_S15, MQDT.Yb171.PARA)
    b = MQDT.basisarray([s], [MQDT.Yb171.FMODEL_HIGHN_S15])
    m = diag(MQDT.matrix_element(MQDT.Yb171.PARA, b))
end
