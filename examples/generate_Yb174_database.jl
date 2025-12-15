using MQDT

# load Yb174 data
parameters = MQDT.Yb174.PARA
low_n_models = MQDT.Yb174.FMODEL_LOWN_P1
low_l_models = [
    MQDT.Yb174.FMODEL_HIGHN_S0,
    MQDT.Yb174.FMODEL_HIGHN_S1,
    MQDT.Yb174.FMODEL_HIGHN_P0,
    MQDT.Yb174.FMODEL_HIGHN_P1,
    MQDT.Yb174.FMODEL_HIGHN_P2,
    MQDT.Yb174.FMODEL_HIGHN_D1,
    MQDT.Yb174.FMODEL_HIGHN_D2,
    MQDT.Yb174.FMODEL_HIGHN_D3,
    MQDT.Yb174.FMODEL_HIGHN_F2,
    MQDT.Yb174.FMODEL_HIGHN_F3,
    MQDT.Yb174.FMODEL_HIGHN_F4,
    MQDT.Yb174.FMODEL_HIGHN_G3,
    MQDT.Yb174.FMODEL_HIGHN_G4,
    MQDT.Yb174.FMODEL_HIGHN_G5,
]

# bounds
n_min = [3, 27, 6, 6, 5, 26, 5, 18, 25, 7, 25, 25, 25, 25]
n_max = 30

# calculate low \nu MQDT states
low_n_states = eigenstates(2, 2, low_n_models, parameters)

# calculate high \nu, low \ell MQDT states
low_l_states = [eigenstates(n_min[i], n_max, low_l_models[i], parameters) for i in eachindex(n_min)]

# calculate high \ell SQDT states
l_max = n_max - 1
MQDT.wigner_init_float(n_max, "Jmax", 9) # initialize Wigner symbol calculation
high_l_models = single_channel_jj_models(:Yb174, 5:l_max)
high_l_states = [eigenstates(25, n_max, M, parameters) for M in high_l_models]

# generate basis and calculate matrix elements
basis = basisarray(vcat(low_n_states, low_l_states, high_l_states), vcat(low_n_models, low_l_models, high_l_models))
@time me = matrix_elements(basis, parameters)

# prepare tables
using DataFrames
col_names = [:id_initial, :id_final, :value]

e1 = DataFrame(me["dipole"], col_names)
e2 = DataFrame(me["quadrupole"], col_names)
m1 = DataFrame(me["paramagnetic"], col_names)
m2 = DataFrame(me["diamagnetic"], col_names)

states_table = DataFrame(;
    id=collect(1:size(basis)),
    energy=MQDT.get_e(basis, parameters),
    parity=MQDT.get_p(basis),
    f=MQDT.get_f(basis),
    nu=MQDT.get_nu(basis),
    term=MQDT.get_term(basis),
    lead=MQDT.get_lead(basis),
)
sort!(states_table, [:nu])

# store tables as csv files
using CSV
CSV.write("Yb174_mqdt_states.csv", states_table)
CSV.write("Yb174_mqdt_matrix_elements_d.csv", e1)
CSV.write("Yb174_mqdt_matrix_elements_q.csv", e2)
CSV.write("Yb174_mqdt_matrix_elements_mu.csv", m1)
CSV.write("Yb174_mqdt_matrix_elements_q0.csv", m2)

# store tables as parquet files
using Parquet2
Parquet2.writefile("Yb174_mqdt_states.parquet", states_table)
Parquet2.writefile("Yb174_mqdt_matrix_elements_d.parquet", e1)
Parquet2.writefile("Yb174_mqdt_matrix_elements_q.parquet", e2)
Parquet2.writefile("Yb174_mqdt_matrix_elements_mu.parquet", m1)
Parquet2.writefile("Yb174_mqdt_matrix_elements_q0.parquet", m2)
