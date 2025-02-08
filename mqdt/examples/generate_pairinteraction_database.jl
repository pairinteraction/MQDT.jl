# setup environment
using Pkg
Pkg.activate(".")

# load mqdt module
include("$(pwd())/mqdt/src/mqdt.jl")
using .mqdt

# choose species
SP = :Yb174
include("$(pwd())/mqdt/parameters/$(SP).jl")

# calculate bound states
n1 = 35 
n2 = 40 # above n ~ 105 will have cause problems with the current integrator for matrix elements
models = [
    MODEL_S0, MODEL_S1, 
    MODEL_P0, MODEL_P1, MODEL_P2, 
    MODEL_D1, MODEL_D2, MODEL_D3
    ]
c_states = eigenstates(1.5, 2.5, MODEL_G1, PARA) # clock states
r_states = Vector{EigenStates}() # Rydberg states
for i in eachindex(models)
    push!(r_states, eigenstates(n1, n2, models[i], PARA)) 
end

# generate state table
c_basis = basisarray(c_states, MODEL_G1)
r_basis = basisarray(r_states, models)
basis = union(c_basis, r_basis)
state_table = state_data(basis, PARA)

# calculate matrix elements
@time d1 = matrix_element(1, basis) # dipole
@time d2 = matrix_element(2, basis) # quadrupole
@time dm = matrix_element(PARA, basis) # Zeeman
@time dd = matrix_element(basis) # diamagnetic

# generate matrix element table
m1 = matrix_data(d1)
m2 = matrix_data(d2)
mm = matrix_data(dm)
md = matrix_data(dd)

# prepare PAIRINTERACTION output
c_db = databasearray(c_states, MODEL_G1)
r_db = databasearray(r_states, models)
b_db = union(c_db, r_db)
ST = state_data(b_db, PARA)

# store full matrix for PAIRINTERACTION (as opposed to upper triangle)
M1 = tri_to_full(m1, ST)
M2 = tri_to_full(m2, ST)
MM = tri_to_full(mm, ST)
MD = tri_to_full(md, ST)

# store tables as csv files
using CSV
CSV.write("$(pwd())/data/$(SP)_mqdt_states.csv", ST)
CSV.write("$(pwd())/data/$(SP)_mqdt_matrix_elements_d.csv", M1)
CSV.write("$(pwd())/data/$(SP)_mqdt_matrix_elements_q.csv", M2)
CSV.write("$(pwd())/data/$(SP)_mqdt_matrix_elements_mu.csv", MM)
CSV.write("$(pwd())/data/$(SP)_mqdt_matrix_elements_q0.csv", MD)

# store tables as parquet files
using Parquet2
Parquet2.writefile("$(pwd())/data/$(SP)_mqdt_states.parquet", ST)
Parquet2.writefile("$(pwd())/data/$(SP)_mqdt_matrix_elements_d.parquet", M1)
Parquet2.writefile("$(pwd())/data/$(SP)_mqdt_matrix_elements_q.parquet", M2)
Parquet2.writefile("$(pwd())/data/$(SP)_mqdt_matrix_elements_mu.parquet", MM)
Parquet2.writefile("$(pwd())/data/$(SP)_mqdt_matrix_elements_q0.parquet", MD)
