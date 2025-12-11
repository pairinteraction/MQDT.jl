# --------------------------------------------------------
# get basis state properties
# --------------------------------------------------------

function get_e(T::BasisArray, P::Parameters)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = P.threshold - P.rydberg ./ T.states[i].energy .^ 2
    end
    return t
end

function get_p(T::BasisArray)
    t = Vector{Int}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].parity
    end
    return t
end

function get_f(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].f
    end
    return t
end

function get_nu(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].energy
    end
    return t
end

function get_term(T::BasisArray)
    t = Vector{String}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].term
    end
    return t
end

function get_lead(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].lead
    end
    return t
end
