using MQDT

# choose species
SP = :Yb174

# generate dispatch tables
const MODELS_TABLE = Dict(
    :Sr87 => [MQDT.Sr87.MODEL_S0, MQDT.Sr87.MODEL_S1, 
              MQDT.Sr87.MODEL_P0, MQDT.Sr87.MODEL_P1, MQDT.Sr87.MODEL_P2, 
              MQDT.Sr87.MODEL_D1, MQDT.Sr87.MODEL_D2, MQDT.Sr87.MODEL_D3],
    :Sr88 => [MQDT.Sr88.MODEL_S0, MQDT.Sr88.MODEL_S1,
              MQDT.Sr88.MODEL_P0, MQDT.Sr88.MODEL_P1, MQDT.Sr88.MODEL_P2, 
              MQDT.Sr88.MODEL_D1, MQDT.Sr88.MODEL_D2, MQDT.Sr88.MODEL_D3],
    :Yb171 => [MQDT.Yb171.MODEL_S0, MQDT.Yb171.MODEL_S1,
               MQDT.Yb171.MODEL_P0, MQDT.Yb171.MODEL_P1, MQDT.Yb171.MODEL_P2,
               MQDT.Yb171.MODEL_D1, MQDT.Yb171.MODEL_D2, MQDT.Yb171.MODEL_D3],
    :Yb173 => [MQDT.Yb173.MODEL_S0, MQDT.Yb173.MODEL_S1,
               MQDT.Yb173.MODEL_P0, MQDT.Yb173.MODEL_P1, MQDT.Yb173.MODEL_P2,
               MQDT.Yb173.MODEL_D1, MQDT.Yb173.MODEL_D2, MQDT.Yb173.MODEL_D3],
    :Yb174 => [MQDT.Yb174.MODEL_S0, MQDT.Yb174.MODEL_S1,
               MQDT.Yb174.MODEL_P0, MQDT.Yb174.MODEL_P1, MQDT.Yb174.MODEL_P2,
               MQDT.Yb174.MODEL_D1, MQDT.Yb174.MODEL_D2, MQDT.Yb174.MODEL_D3],
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
c_states = MQDT.eigenstates(1.5, 2.5, MODEL_G1, para) # clock states
r_states = [MQDT.eigenstates(n1, n2, M, para) for M in models] # Rydberg states

# generate state table
c_basis = MQDT.basisarray(c_states, MODEL_G1)
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
c_db = MQDT.databasearray(c_states, MODEL_G1)
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
CSV.write("data/$(SP)_mqdt_states.csv", ST)
CSV.write("data/$(SP)_mqdt_matrix_elements_d.csv", M1)
CSV.write("data/$(SP)_mqdt_matrix_elements_q.csv", M2)
CSV.write("data/$(SP)_mqdt_matrix_elements_mu.csv", MM)
CSV.write("data/$(SP)_mqdt_matrix_elements_q0.csv", MD)

# store tables as parquet files
using Parquet2
Parquet2.writefile("data/$(SP)_mqdt_states.parquet", ST)
Parquet2.writefile("data/$(SP)_mqdt_matrix_elements_d.parquet", M1)
Parquet2.writefile("data/$(SP)_mqdt_matrix_elements_q.parquet", M2)
Parquet2.writefile("data/$(SP)_mqdt_matrix_elements_mu.parquet", MM)
Parquet2.writefile("data/$(SP)_mqdt_matrix_elements_q0.parquet", MD)
