using MQDT

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

# calculate low \nu MQDT states
low_n_states = eigenstates(2, 2, low_n_models, parameters)

# calculate high \nu, low \ell MQDT states
n_min, n_max = 25, 30
low_l_states = [eigenstates(n_min, n_max, M, parameters) for M in low_l_models]

# calculate high \ell SQDT states
l_max = n_max - 1
MQDT.wigner_init_float(n_max, "Jmax", 9) # initialize Wigner symbol caluclation
high_l_models = single_channel_models(5:l_max, parameters)
high_l_states = [eigenstates(n_min, n_max, M, parameters) for M in high_l_models]

# generate state table
basis = basisarray(
    vcat(low_n_states, low_l_states, high_l_states),
    vcat(low_n_models, low_l_models, high_l_models),
)
state_table = state_data(basis, parameters)

# calculate matrix elements
@time d1 = matrix_element(1, basis) # dipole
@time d2 = matrix_element(2, basis) # quadrupole
@time dm = matrix_element(parameters, basis) # Zeeman
@time dd = matrix_element(basis) # diamagnetic

# generate matrix element table
m1 = matrix_data(d1)
m2 = matrix_data(d2)
mm = matrix_data(dm)
md = matrix_data(dd)

# prepare PAIRINTERACTION state table
db = databasearray(
    vcat(low_n_states, low_l_states, high_l_states),
    vcat(low_n_models, low_l_models, high_l_models),
)
st = state_data(db, parameters)

# store tables as csv files
using CSV
CSV.write("Yb174_mqdt_states.csv", st)
CSV.write("Yb174_mqdt_matrix_elements_d.csv", m1)
CSV.write("Yb174_mqdt_matrix_elements_q.csv", m2)
CSV.write("Yb174_mqdt_matrix_elements_mu.csv", mm)
CSV.write("Yb174_mqdt_matrix_elements_q0.csv", md)

# store tables as parquet files
using Parquet2
Parquet2.writefile("Yb174_mqdt_states.parquet", st)
Parquet2.writefile("Yb174_mqdt_matrix_elements_d.parquet", m1)
Parquet2.writefile("Yb174_mqdt_matrix_elements_q.parquet", m2)
Parquet2.writefile("Yb174_mqdt_matrix_elements_mu.parquet", mm)
Parquet2.writefile("Yb174_mqdt_matrix_elements_q0.parquet", md)
