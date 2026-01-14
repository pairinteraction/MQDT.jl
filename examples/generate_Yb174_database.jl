using MQDT

species = :Yb174
parameters = MQDT.get_species_parameters(species)
n_max = 30
n_min_high_l = 25 # minimum n for high-l states

all_models = Vector{MQDT.fModel}()

sr = 1 / 2
Jc = 1 / 2
ic = parameters.spin
for lr in 0:(n_max - 1)
    for Jr in abs(lr - sr):1:(lr + sr)
        for Fc in abs(Jc - ic):1:(Jc + ic)
            for F in abs(Fc - Jr):1:(Fc + Jr)
                models = MQDT.get_fmodels(species, lr, Jr, Fc, F, parameters)
                for model in models
                    if !any(m -> m.name == model.name, all_models)
                        push!(all_models, model)
                    end
                end
            end
        end
    end
end

@info "Calculating MQDT states..."
states = Vector{MQDT.EigenStates}(undef, length(all_models))
for (i, M) in enumerate(all_models)
    n_min = NaN
    if startswith(M.name, "SQDT")
        n_min = n_min_high_l
    end
    @info "$(M.name)"
    states[i] = MQDT.eigenstates(n_min, n_max, M, parameters)
    @info "  found nu_min=$(minimum(states[i].n)), nu_max=$(maximum(states[i].n)), total states=$(length(states[i].n))"
end

# generate basis
basis = basisarray(states, all_models)

# calculate matrix elements
MQDT.wigner_init_float(n_max, "Jmax", 9) # initialize Wigner symbol calculation
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
