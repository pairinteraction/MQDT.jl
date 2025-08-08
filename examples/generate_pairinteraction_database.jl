using MQDT

const MODELS_TABLE = Dict(
    :Sr87 => [
        MQDT.Sr87.FMODEL_HIGHN_S35,
        MQDT.Sr87.FMODEL_HIGHN_S45,
        MQDT.Sr87.FMODEL_HIGHN_S55,
        MQDT.Sr87.FMODEL_HIGHN_P25,
        MQDT.Sr87.FMODEL_HIGHN_P35,
        MQDT.Sr87.FMODEL_HIGHN_P45,
        MQDT.Sr87.FMODEL_HIGHN_P55,
        MQDT.Sr87.FMODEL_HIGHN_P65,
        MQDT.Sr87.FMODEL_HIGHN_D15,
        MQDT.Sr87.FMODEL_HIGHN_D25,
        MQDT.Sr87.FMODEL_HIGHN_D35,
        MQDT.Sr87.FMODEL_HIGHN_D45,
        MQDT.Sr87.FMODEL_HIGHN_D55,
        MQDT.Sr87.FMODEL_HIGHN_D65,
        MQDT.Sr87.FMODEL_HIGHN_D75,
    ],
    :Sr88 => [
        MQDT.Sr88.FMODEL_HIGHN_S0,
        MQDT.Sr88.FMODEL_HIGHN_S1,
        MQDT.Sr88.FMODEL_HIGHN_P0,
        MQDT.Sr88.FMODEL_HIGHN_P1,
        MQDT.Sr88.FMODEL_HIGHN_P2,
        MQDT.Sr88.FMODEL_HIGHN_D1,
        MQDT.Sr88.FMODEL_HIGHN_D2,
        MQDT.Sr88.FMODEL_HIGHN_D3,
    ],
    :Yb171 => [
        MQDT.Yb171.FMODEL_HIGHN_S05,
        MQDT.Yb171.FMODEL_HIGHN_S15,
        MQDT.Yb171.FMODEL_HIGHN_P05,
        MQDT.Yb171.FMODEL_HIGHN_P15,
        MQDT.Yb171.FMODEL_HIGHN_D05,
        MQDT.Yb171.FMODEL_HIGHN_D15,
        MQDT.Yb171.FMODEL_HIGHN_D25,
        MQDT.Yb171.FMODEL_HIGHN_D35,
        MQDT.Yb171.FMODEL_HIGHN_F25,
        MQDT.Yb171.FMODEL_HIGHN_F35,
        MQDT.Yb171.FMODEL_HIGHN_F45,
        MQDT.Yb171.FMODEL_HIGHN_G25,
        MQDT.Yb171.FMODEL_HIGHN_G35,
        MQDT.Yb171.FMODEL_HIGHN_G45,
        MQDT.Yb171.FMODEL_HIGHN_G55,
    ],
    :Yb173 => [
        MQDT.Yb173.FMODEL_HIGHN_S15,
        MQDT.Yb173.FMODEL_HIGHN_S25,
        MQDT.Yb173.FMODEL_HIGHN_S35,
        MQDT.Yb173.FMODEL_HIGHN_P05,
        MQDT.Yb173.FMODEL_HIGHN_P15,
        MQDT.Yb173.FMODEL_HIGHN_P25,
        MQDT.Yb173.FMODEL_HIGHN_P35,
        MQDT.Yb173.FMODEL_HIGHN_P45,
    ],
    :Yb174 => [
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
    ],
)
const PARA_TABLE = Dict(
    :Sr87 => MQDT.Sr87.PARA,
    :Sr88 => MQDT.Sr88.PARA,
    :Yb171 => MQDT.Yb171.PARA,
    :Yb173 => MQDT.Yb173.PARA,
    :Yb174 => MQDT.Yb174.PARA,
)

# choose species
species = :Yb174

# calculate low \ell MQDT states
n_min, n_max = 35, 40
low_l_models = MODELS_TABLE[species]
parameters = PARA_TABLE[species]
low_l_states = [eigenstates(n_min, n_max, M, parameters) for M in low_l_models]

# calculate high \ell SQDT states
l_max = n_max - 1
MQDT.wigner_init_float(n_max, "Jmax", 9) # initialize Wigner symbol caluclation
high_l_models = single_channel_models(5:l_max, parameters)
high_l_states = [eigenstates(n_min, n_max, M, parameters) for M in high_l_models]

# generate state table
basis = basisarray(vcat(low_l_states, high_l_states), vcat(low_l_models, high_l_models))
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
db = databasearray(vcat(low_l_states, high_l_states), vcat(low_l_models, high_l_models))
st = state_data(db, parameters)

# store tables as csv files
using CSV
CSV.write("$(species)_mqdt_states.csv", st)
CSV.write("$(species)_mqdt_matrix_elements_d.csv", m1)
CSV.write("$(species)_mqdt_matrix_elements_q.csv", m2)
CSV.write("$(species)_mqdt_matrix_elements_mu.csv", mm)
CSV.write("$(species)_mqdt_matrix_elements_q0.csv", md)

# store tables as parquet files
using Parquet2
Parquet2.writefile("$(species)_mqdt_states.parquet", st)
Parquet2.writefile("$(species)_mqdt_matrix_elements_d.parquet", m1)
Parquet2.writefile("$(species)_mqdt_matrix_elements_q.parquet", m2)
Parquet2.writefile("$(species)_mqdt_matrix_elements_mu.parquet", mm)
Parquet2.writefile("$(species)_mqdt_matrix_elements_q0.parquet", md)
