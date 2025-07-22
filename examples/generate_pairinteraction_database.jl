using MQDT

# choose species
species = :Yb174

# generate dispatch tables # TODO check whether the MODELS_TABLE is not missing any models
const MODELS_TABLE = Dict(
    :Sr87 => [
        MQDT.Sr87.FMODEL_HIGHN_S35,
        MQDT.Sr87.FMODEL_HIGHN_S45,
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
        MQDT.Yb171.FMODEL_HIGHN_P25,
        MQDT.Yb171.FMODEL_HIGHN_D05,
        MQDT.Yb171.FMODEL_HIGHN_D15,
        MQDT.Yb171.FMODEL_HIGHN_D25,
        MQDT.Yb171.FMODEL_HIGHN_D35,
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
    ],
)
const PARA_TABLE = Dict(
    :Sr87 => MQDT.Sr87.PARA,
    :Sr88 => MQDT.Sr88.PARA,
    :Yb171 => MQDT.Yb171.PARA,
    :Yb173 => MQDT.Yb173.PARA,
    :Yb174 => MQDT.Yb174.PARA,
)

# calculate bound states
n_min, n_max = 35, 40
models = MODELS_TABLE[species]
parameters = PARA_TABLE[species]
clock_states = MQDT.eigenstates(1.5, 2.5, MQDT.Yb174.FMODEL_LOWN_P0, parameters)  # TODO  make it work for species != Yb174
rydberg_states = [MQDT.eigenstates(n_min, n_max, M, parameters) for M in models]

# generate state table
clock_basis = MQDT.basisarray(clock_states, MQDT.Yb174.FMODEL_LOWN_P0) # TODO make it work for species != Yb174
rydberg_basis = MQDT.basisarray(rydberg_states, models)
basis = union(clock_basis, rydberg_basis)
state_table = MQDT.state_data(basis, parameters)

# calculate matrix elements
@time d1 = MQDT.matrix_element(1, basis) # dipole
@time d2 = MQDT.matrix_element(2, basis) # quadrupole
@time dm = MQDT.matrix_element(parameters, basis) # Zeeman
@time dd = MQDT.matrix_element(basis) # diamagnetic

# generate matrix element table
m1 = MQDT.matrix_data(d1)
m2 = MQDT.matrix_data(d2)
mm = MQDT.matrix_data(dm)
md = MQDT.matrix_data(dd)

# prepare PAIRINTERACTION output
c_db = MQDT.databasearray(clock_states, MQDT.Yb174.FMODEL_LOWN_P0) # TODO make it work for species != Yb174
r_db = MQDT.databasearray(rydberg_states, models)
b_db = union(c_db, r_db)
ST = MQDT.state_data(b_db, parameters)

# store full matrix for PAIRINTERACTION (as opposed to upper triangle)
M1 = MQDT.tri_to_full(m1, ST)
M2 = MQDT.tri_to_full(m2, ST)
MM = MQDT.tri_to_full(mm, ST)
MD = MQDT.tri_to_full(md, ST)

# store tables as csv files
using CSV
CSV.write("$(species)_mqdt_states.csv", ST)
CSV.write("$(species)_mqdt_matrix_elements_d.csv", M1)
CSV.write("$(species)_mqdt_matrix_elements_q.csv", M2)
CSV.write("$(species)_mqdt_matrix_elements_mu.csv", MM)
CSV.write("$(species)_mqdt_matrix_elements_q0.csv", MD)

# store tables as parquet files
using Parquet2
Parquet2.writefile("$(species)_mqdt_states.parquet", ST)
Parquet2.writefile("$(species)_mqdt_matrix_elements_d.parquet", M1)
Parquet2.writefile("$(species)_mqdt_matrix_elements_q.parquet", M2)
Parquet2.writefile("$(species)_mqdt_matrix_elements_mu.parquet", MM)
Parquet2.writefile("$(species)_mqdt_matrix_elements_q0.parquet", MD)
