using LinearAlgebra
using SparseArrays
using CGcoefficient
using DataFrames
using MQDT
using Plots

default(
    dpi = 400,
    fontfamily = "sans-serif",
    titlefontsize = 12,
    framestyle = :box,
    labels = false,
)

# LS g-factors
g_L(J, L, S) = (J*(J+1) + L*(L+1) - S*(S+1)) / (2*J*(J+1))
g_J(J, L, S) = 2 * (J*(J+1) - L*(L+1) + S*(S+1)) / (2*J*(J+1)) + g_L(J, L, S)
function g_F(F, I, J, L, S)
    g = 0
    if J != 0
        g = g_J(J, L, S) * (F*(F+1) - I*(I+1) + J*(J+1)) / (2*F*(F+1))
    end
    return g
end

# wigner3j symbols
wigner_init_float(10, "Jmax", 3)
function wigner3j(k::Int, q::Int, f1::Number, m1::Number, f2::Number, m2::Number)
    a = 0.0
    if abs(q) <= k && abs(m1) <= f1 && abs(m2) <= f2 && abs(f1-f2) <= k && m1-m2 == q
        qn = Vector{Int64}(2*[f1, k, f2, -m1, q, m2])
        a = (-1.0)^(f1-m1) * f3j(qn...)
    end
    return a
end

# S series of Yb171
f = 0.5
s_S05 = eigenstates(24, 129, MQDT.Yb171.FMODEL_HIGHN_S05, MQDT.Yb171.PARA)
b_S05 = basisarray([s_S05], [MQDT.Yb171.FMODEL_HIGHN_S05])
m_S05 = diag(matrix_element(MQDT.Yb171.PARA, b_S05))
g_S05 = -2m_S05 / f * wigner3j(1, 0, f, f, f, f)

f = 1.5
s_S15 = eigenstates(24, 129, MQDT.Yb171.FMODEL_HIGHN_S15, MQDT.Yb171.PARA)
b_S15 = basisarray([s_S15], [MQDT.Yb171.FMODEL_HIGHN_S15])
m_S15 = diag(matrix_element(MQDT.Yb171.PARA, b_S15))
g_S15 = -2m_S15 / f * wigner3j(1, 0, f, f, f, f)

scatter(layout = (2, 1), size = (400, 500), legend = :topleft)
scatter!(
    s_S05.n,
    mod.(-s_S05.nu[7, :], 1),
    title = "¹⁷¹Yb S series",
    subplot = 1,
    label = "F=1/2",
)
scatter!(
    s_S15.n,
    mod.(-s_S15.nu[1, :], 1),
    xlabel = "ν",
    ylabel = "μ",
    subplot = 1,
    label = "F=3/2",
)
scatter!(s_S05.n, g_S05, subplot = 2)
scatter!(s_S15.n, g_S15, xlabel = "ν", ylabel = "g factor", subplot = 2)
hline!([g_F(1/2, 1/2, 0, 0, 0)], l = :dash, label = "¹S₁", subplot = 2)
hline!([g_F(1/2, 1/2, 1, 0, 1)], l = :dash, label = "³S₁", subplot = 2)
hline!([g_F(3/2, 1/2, 1, 0, 1)], l = :dash, subplot = 2)

savefig("Yb171_S_series.pdf")

# P series of Yb171
f = 0.5
s_P05 = eigenstates(10, 70, MQDT.Yb171.FMODEL_HIGHN_P05, MQDT.Yb171.PARA)
b_P05 = basisarray([s_P05], [MQDT.Yb171.FMODEL_HIGHN_P05])
m_P05 = diag(matrix_element(MQDT.Yb171.PARA, b_P05))
g_P05 = -2m_P05 / f * wigner3j(1, 0, f, f, f, f)

f = 1.5
s_P15 = eigenstates(10, 70, MQDT.Yb171.FMODEL_HIGHN_P15, MQDT.Yb171.PARA)
b_P15 = basisarray([s_P15], [MQDT.Yb171.FMODEL_HIGHN_P15])
m_P15 = diag(matrix_element(MQDT.Yb171.PARA, b_P15))
g_P15 = -2m_P15 / f * wigner3j(1, 0, f, f, f, f)

scatter(layout = (2, 2), size = (700, 500))
scatter!(
    s_P05.n,
    mod.(-s_P05.nu[1, :], 1),
    ylabel = "μ",
    title = "¹⁷¹Yb P F=1/2",
    subplot = 1,
)
scatter!(s_P05.n, g_P05, xlabel = "ν", ylabel = "g factor", subplot = 3)
scatter!(
    s_P15.n,
    mod.(-s_P15.nu[1, :], 1),
    ylabel = "μ",
    title = "¹⁷¹Yb P F=3/2",
    subplot = 2,
)
scatter!(s_P15.n, g_P15, xlabel = "ν", ylabel = "g factor", subplot = 4)
hline!([g_F(1/2, 1/2, 1, 1, 0)], l = :dash, label = "¹P₁", subplot = 3)
hline!([g_F(1/2, 1/2, 1, 1, 1)], l = :dash, label = "³P₁", subplot = 3)
hline!([g_F(1/2, 1/2, 0, 1, 1)], l = :dash, label = "³P₀", subplot = 3)
hline!([g_F(3/2, 1/2, 1, 1, 0)], l = :dash, label = "¹P₁", subplot = 4, ylim = (0, 2))
hline!([g_F(3/2, 1/2, 1, 1, 1)], l = :dash, label = "³P₁", subplot = 4)
hline!([g_F(3/2, 1/2, 2, 1, 1)], l = :dash, label = "³P₂", subplot = 4)

savefig("Yb171_P_series.pdf")
