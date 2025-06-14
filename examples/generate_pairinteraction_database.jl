using MQDT

# choose species
SP = :Yb174

# generate dispatch tables # TODO check whether the MODELS_TABLE is not missing any models
const MODELS_TABLE = Dict(
    :Sr87 => [MQDT.Sr87.MODEL_S35, MQDT.Sr87.MODEL_S45, 
              MQDT.Sr87.MODEL_P25, MQDT.Sr87.MODEL_P35, MQDT.Sr87.MODEL_P45, MQDT.Sr87.MODEL_P55, MQDT.Sr87.MODEL_P65,
              MQDT.Sr87.MODEL_D15, MQDT.Sr87.MODEL_D25, MQDT.Sr87.MODEL_D35, MQDT.Sr87.MODEL_D45, MQDT.Sr87.MODEL_D55, MQDT.Sr87.MODEL_D65, MQDT.Sr87.MODEL_D75],
    :Sr88 => [MQDT.Sr88.FMODEL_S0, MQDT.Sr88.FMODEL_S1,
              MQDT.Sr88.FMODEL_P0, MQDT.Sr88.FMODEL_P1, MQDT.Sr88.FMODEL_P2, 
              MQDT.Sr88.FMODEL_D1, MQDT.Sr88.FMODEL_D2, MQDT.Sr88.FMODEL_D3],
    :Yb171 => [MQDT.Yb171.RYDBERG_S05, MQDT.Yb171.RYDBERG_S15,
               MQDT.Yb171.RYDBERG_P05, MQDT.Yb171.RYDBERG_P15, MQDT.Yb171.RYDBERG_P25,
               MQDT.Yb171.RYDBERG_D05, MQDT.Yb171.RYDBERG_D15, MQDT.Yb171.RYDBERG_D25, MQDT.Yb171.RYDBERG_D35],
    :Yb173 => [MQDT.Yb173.MODEL_S15, MQDT.Yb173.MODEL_S25, MQDT.Yb173.MODEL_S35,
               MQDT.Yb173.MODEL_P05, MQDT.Yb173.MODEL_P15, MQDT.Yb173.MODEL_P25, MQDT.Yb173.MODEL_P35, MQDT.Yb173.MODEL_P45],
    :Yb174 => [MQDT.Yb174.RYDBERG_S0, MQDT.Yb174.RYDBERG_S1,
               MQDT.Yb174.RYDBERG_P0, MQDT.Yb174.RYDBERG_P1, MQDT.Yb174.RYDBERG_P2,
               MQDT.Yb174.RYDBERG_D1, MQDT.Yb174.RYDBERG_D2, MQDT.Yb174.RYDBERG_D3],
)
const PARA_TABLE = Dict(
    :Sr87 => MQDT.Sr87.PARA,
    :Sr88 => MQDT.Sr88.PARA,
    :Yb171 => MQDT.Yb171.PARA,
    :Yb173 => MQDT.Yb173.PARA,
    :Yb174 => MQDT.Yb174.PARA,
)

# calculate bound states
n1, n2 = 35, 40 # above n ~ 105 will cause problems with the current integrator for matrix elements
models = MODELS_TABLE[SP]
para = PARA_TABLE[SP]
c_states = MQDT.eigenstates(1.5, 2.5, MQDT.Yb174.CORE_P0, para) # clock states # TODO  make it work for SP != Yb174
r_states = [MQDT.eigenstates(n1, n2, M, para) for M in models] # Rydberg states

# generate state table
c_basis = MQDT.basisarray(c_states, MQDT.Yb174.CORE_P0) # TODO make it work for SP != Yb174
r_basis = MQDT.basisarray(r_states, models)
basis = union(c_basis, r_basis)
state_table = MQDT.state_data(basis, para)

# calculate matrix elements
@time d1 = MQDT.matrix_element(1, basis) # dipole
@time d2 = MQDT.matrix_element(2, basis) # quadrupole
@time dm = MQDT.matrix_element(para, basis) # Zeeman
@time dd = MQDT.matrix_element(basis) # diamagnetic

# generate matrix element table
m1 = MQDT.matrix_data(d1)
m2 = MQDT.matrix_data(d2)
mm = MQDT.matrix_data(dm)
md = MQDT.matrix_data(dd)

# prepare PAIRINTERACTION output
c_db = MQDT.databasearray(c_states, MQDT.Yb174.CORE_P0) # TODO make it work for SP != Yb174
r_db = MQDT.databasearray(r_states, models)
b_db = union(c_db, r_db)
ST = MQDT.state_data(b_db, para)

# store full matrix for PAIRINTERACTION (as opposed to upper triangle)
M1 = MQDT.tri_to_full(m1, ST)
M2 = MQDT.tri_to_full(m2, ST)
MM = MQDT.tri_to_full(mm, ST)
MD = MQDT.tri_to_full(md, ST)

# store tables as csv files
using CSV
CSV.write("$(SP)_mqdt_states.csv", ST)
CSV.write("$(SP)_mqdt_matrix_elements_d.csv", M1)
CSV.write("$(SP)_mqdt_matrix_elements_q.csv", M2)
CSV.write("$(SP)_mqdt_matrix_elements_mu.csv", MM)
CSV.write("$(SP)_mqdt_matrix_elements_q0.csv", MD)

# store tables as parquet files
using Parquet2
Parquet2.writefile("$(SP)_mqdt_states.parquet", ST)
Parquet2.writefile("$(SP)_mqdt_matrix_elements_d.parquet", M1)
Parquet2.writefile("$(SP)_mqdt_matrix_elements_q.parquet", M2)
Parquet2.writefile("$(SP)_mqdt_matrix_elements_mu.parquet", MM)
Parquet2.writefile("$(SP)_mqdt_matrix_elements_q0.parquet", MD)
