using MQDT

@testset "small test" begin
    s = MQDT.eigenstates(24, 30, MQDT.Yb171.FMODEL_HIGHN_S15, MQDT.Yb171.PARA)
    b = MQDT.basisarray([s], [MQDT.Yb171.FMODEL_HIGHN_S15])
    m = [MQDT.multipole_moments(i, i, MQDT.Yb171.PARA)[4] for i in b.states]
end
