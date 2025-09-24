# --------------------------------------------------------
# data frames
# --------------------------------------------------------

"""
    state_data(T::BasisArray, P::Parameters)
    state_data(T::DataBaseArray, P::Parameters)

Store a bound state data as a data frame as used by the PAIRINTERACTION software.
"""
function state_data(T::BasisArray, P::Parameters)
    df = DataFrame(
        id=collect(1:size(T)),
        energy=get_e(T, P),
        parity=get_p(T),
        f=get_f(T),
        nu=get_nu(T),
        l=exp_l(T),
        term=get_term(T),
        lead=get_lead(T),
        L=get_L(T),
        S=get_S(T),
    )
    return sort!(df, [:nu, :l])
end

function state_data(T::DataBaseArray, P::Parameters)
    df = DataFrame(
        id=collect(1:size(T)),
        energy=get_e(T, P),
        parity=get_p(T),
        n=get_n(T, P),
        nu=get_nu(T),
        f=get_f(T),
        exp_nui=exp_nui(T),
        exp_l=exp_L(T),
        exp_j=exp_J(T),
        exp_s=exp_S(T),
        exp_l_ryd=exp_lr(T),
        exp_j_ryd=exp_Jr(T),
        std_nui=std_nui(T),
        std_l=std_L(T),
        std_j=std_J(T),
        std_s=std_S(T),
        std_l_ryd=std_lr(T),
        std_j_ryd=std_Jr(T),
        is_j_total_momentum=is_J(T, P),
        is_calculated_with_mqdt=is_mqdt(T),
        underspecified_channel_contribution=get_neg(T),
        leading_term=get_term(T),
        leading_percentage=get_percentage(T),
    )
    return sort!(df, [:nu, :exp_l_ryd])
end

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

function get_L(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].L
    end
    return t
end

function get_S(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].S
    end
    return t
end

function exp_q(q::Vector, n::Vector)
    if allequal(q)
        return Float64(q[1])
    else
        m = n .^ 2
        M = sum(m)
        if M > 1
            m /= M
        end
        return dot(q, m)
    end
end

function std_q(q::Vector, n::Vector)
    if allequal(q)
        return 0.0
    else
        m = n .^ 2
        M = sum(m)
        if M > 1
            m /= M
        end
        e1 = dot(q, m)^2
        e2 = dot(q .^ 2, m)
        if abs(e1 - e2) < 1e-11
            return 0.0
        else
            return sqrt(e2 - e1)
        end
    end
end

function exp_l(T::BasisState)
    l = T.lr
    n = T.coeff
    return exp_q(l, n)
end

function exp_l(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_l(T.states[i])
    end
    return t
end

# --------------------------------------------------------
# Pairinteraction output utility
# --------------------------------------------------------

function get_e(T::DataBaseArray, P::Parameters)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = P.threshold - P.rydberg ./ T.states[i].nu .^ 2
    end
    return t / 219474.6313632 # convert 1/cm to atomic units
end

function get_p(T::DataBaseArray)
    t = Vector{Int}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].parity
    end
    return t
end

function get_n(T::DataBaseArray, P::Parameters)
    nu = get_nu(T)
    l = round.(Int, exp_lr(T))
    return get_n(nu, l, P.species)
end

function get_n(nu::Vector{Float64}, l::Vector{Int}, species::Symbol)
    i0 = findall(iszero, l)
    i1 = findall(iszero, l .- 1)
    i2 = findall(iszero, l .- 2)
    i3 = findall(iszero, l .- 3)
    j0 = findall(x->x<2, nu)
    nu[j0] .+= 1
    if occursin("Yb", String(species))
        nu[i0] .+= 4
        nu[i1] .+= 3
        nu[i2] .+= 2
        nu[i3] .+= 1
    else
        nu[i0] .+= 3
        nu[i1] .+= 2
        nu[i2] .+= 2
    end
    return ceil.(Int, nu)
end

function get_nu(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].nu
    end
    return t
end

function get_f(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].f
    end
    return t
end

function get_neg(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].neg
    end
    return t
end

function exp_nui(T::DataBaseState)
    ν = T.nui_all
    n = T.ai_all
    return exp_q(ν, n)
end

function exp_nui(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_nui(T.states[i])
    end
    return t
end

function std_nui(T::DataBaseState)
    ν = T.nui_all
    n = T.ai_all
    return std_q(ν, n)
end

function std_nui(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = std_nui(T.states[i])
    end
    return t
end

function exp_L(T::DataBaseState)
    l = T.L
    n = T.transform' * T.ai
    return exp_q(l, n)
end

function exp_L(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_L(T.states[i])
    end
    return t
end

function std_L(T::DataBaseState)
    l = T.L
    n = T.transform' * T.ai
    return std_q(l, n)
end

function std_L(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = std_L(T.states[i])
    end
    return t
end

function exp_J(T::DataBaseState)
    j = T.J
    n = T.transform' * T.ai
    return exp_q(j, n)
end

function exp_J(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_J(T.states[i])
    end
    return t
end

function std_J(T::DataBaseState)
    j = T.J
    n = T.transform' * T.ai
    return std_q(j, n)
end

function std_J(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = std_J(T.states[i])
    end
    return t
end

function exp_S(T::DataBaseState)
    s = T.S
    n = T.transform' * T.ai
    return exp_q(s, n)
end

function exp_S(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_S(T.states[i])
    end
    return t
end

function std_S(T::DataBaseState)
    s = T.S
    n = T.transform' * T.ai
    return std_q(s, n)
end

function std_S(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = std_S(T.states[i])
    end
    return t
end

function exp_lr(T::DataBaseState)
    l = T.lr
    n = T.ai
    return exp_q(l, n)
end

function exp_lr(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_lr(T.states[i])
    end
    return t
end

function std_lr(T::DataBaseState)
    l = T.lr
    n = T.ai
    return std_q(l, n)
end

function std_lr(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = std_lr(T.states[i])
    end
    return t
end

function exp_Jr(T::DataBaseState)
    j = T.Jr
    n = T.ai
    return exp_q(j, n)
end

function exp_Jr(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_Jr(T.states[i])
    end
    return t
end

function std_Jr(T::DataBaseState)
    j = T.Jr
    n = T.ai
    return std_q(j, n)
end

function std_Jr(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = std_Jr(T.states[i])
    end
    return t
end

function is_J(T::DataBaseArray, P::Parameters)
    return repeat([iszero(P.spin)], size(T))
end

function is_mqdt(T::DataBaseState)
    return !isone(length(T.nui_all))
end

function is_mqdt(T::DataBaseArray)
    t = Vector{Bool}(undef, size(T))
    for i in eachindex(t)
        t[i] = is_mqdt(T.states[i])
    end
    return t
end

function get_term(T::DataBaseArray)
    t = Vector{String}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].term
    end
    return t
end

function get_percentage(T::DataBaseArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = T.states[i].lead
    end
    return t
end
