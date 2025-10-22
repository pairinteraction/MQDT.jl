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

# --------------------------------------------------------
# Pairinteraction output utility
# --------------------------------------------------------

function get_n(T::BasisArray, P::Parameters)
    nu = get_nu(T)
    l = round.(Int, calc_exp_qn(T, "l_r"))
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

function get_neg(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for (i, state) in enumerate(T.states)
        neglected_inds = findall(iszero, state.core)
        ai_neglected = state.ai_all[neglected_inds]
        neglected = sum(ai_neglected .^ 2)
        t[i] = neglected
    end
    return t
end

function exp_nui(T::BasisState)
    nu = T.nui_all
    n = T.ai_all
    return exp_q(nu, n)
end

function exp_nui(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_nui(T.states[i])
    end
    return t
end

function std_nui(T::BasisState)
    nu = T.nui_all
    n = T.ai_all
    return std_q(nu, n)
end

function std_nui(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = std_nui(T.states[i])
    end
    return t
end

function is_J(T::BasisArray, P::Parameters)
    return repeat([iszero(P.spin)], size(T))
end

function is_mqdt(T::BasisState)
    return !isone(length(T.nui_all))
end

function is_mqdt(T::BasisArray)
    t = Vector{Bool}(undef, size(T))
    for i in eachindex(t)
        t[i] = is_mqdt(T.states[i])
    end
    return t
end
