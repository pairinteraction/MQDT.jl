using MQDT

@testset "small test" begin
    s = eigenstates(24, 30, MQDT.Yb171.FMODEL_HIGHN_S15, MQDT.Yb171.PARAMETERS)
    b = basisarray([s], [MQDT.Yb171.FMODEL_HIGHN_S15])
    m = [MQDT.multipole_moments(i, i, MQDT.Yb171.PARAMETERS)[3] for i in b.states]
    @test isapprox(m, 1e6*[2.491, 2.887, 3.329, 3.82], atol=1e3)
end
