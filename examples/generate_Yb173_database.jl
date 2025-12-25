using MQDT

# load Yb173 data
parameters = MQDT.Yb173.PARA
low_l_models = [
    MQDT.Yb173.FMODEL_HIGHN_S15,
    MQDT.Yb173.FMODEL_HIGHN_S25,
    MQDT.Yb173.FMODEL_HIGHN_S35,
    MQDT.Yb173.FMODEL_HIGHN_P05,
    MQDT.Yb173.FMODEL_HIGHN_P15,
    MQDT.Yb173.FMODEL_HIGHN_P25,
    MQDT.Yb173.FMODEL_HIGHN_P35,
    MQDT.Yb173.FMODEL_HIGHN_P45,
]

# bounds
n_min = 20
n_max = 30

# calculate high \nu, low \ell MQDT states
low_l_states = [eigenstates(n_min, n_max, low_l_models[i], parameters) for i in eachindex(low_l_models)]

# calculate high \ell SQDT states
l_max = n_max - 1
MQDT.wigner_init_float(n_max, "Jmax", 9) # initialize Wigner symbol caluclation
high_l_models = single_channel_fj_models(:Yb173, 5:l_max, parameters)
high_l_states = [eigenstates(25, n_max, M, parameters) for M in high_l_models]

# generate basis and calculate matrix elements
basis = basisarray(vcat(low_l_states, high_l_states), vcat(low_l_models, high_l_models))
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
output_dir = "data/Yb173_mqdt/"
mkpath(output_dir)

using CSV
CSV.write("$(output_dir)states.csv", states_table)
CSV.write("$(output_dir)matrix_elements_d.csv", e1)
CSV.write("$(output_dir)matrix_elements_q.csv", e2)
CSV.write("$(output_dir)matrix_elements_mu.csv", m1)
CSV.write("$(output_dir)matrix_elements_q0.csv", m2)

# store tables as parquet files
using Parquet2
Parquet2.writefile("$(output_dir)states.parquet", states_table)
Parquet2.writefile("$(output_dir)matrix_elements_d.parquet", e1)
Parquet2.writefile("$(output_dir)matrix_elements_q.parquet", e2)
Parquet2.writefile("$(output_dir)matrix_elements_mu.parquet", m1)
Parquet2.writefile("$(output_dir)matrix_elements_q0.parquet", m2)
