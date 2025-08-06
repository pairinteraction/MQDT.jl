# --------------------------------------------------------
# data frames
# --------------------------------------------------------

"""
See also [`matrix_data`](@ref)

    state_data(T::BasisArray, P::Parameters)
    state_data(T::DataBaseArray, P::Parameters)

Store a bound state data as a data frame as used by the PAIRINTERACTION software.
"""
function state_data(T::BasisArray, P::Parameters)
    df = DataFrame(
        id = collect(1:size(T)),
        energy = get_e(T, P),
        parity = get_p(T),
        f = get_f(T),
        nu = get_nu(T),
        l = exp_l(T),
        S = exp_S(T)
        )
    S = sortperm(get_e(T, P))
    return df[S,:]
end

function state_data(T::DataBaseArray, P::Parameters)
    df = DataFrame(
        id = collect(1:size(T)),
        energy = get_e(T, P),
        parity = get_p(T),
        n = get_n(T, P),
        nu = get_nu(T),
        f = get_f(T),
        exp_nui = exp_nui(T),
        exp_l = exp_L(T),
        exp_j = exp_J(T),
        exp_s = exp_S(T),
        exp_l_ryd = exp_lr(T),
        exp_j_ryd = exp_Jr(T),
        std_nui = std_nui(T),
        std_l = std_L(T),
        std_j = std_J(T),
        std_s = std_S(T),
        std_l_ryd = std_lr(T),
        std_j_ryd = std_Jr(T),
        is_j_total_momentum = is_J(T, P), 
        is_calculated_with_mqdt = is_mqdt(T),
        underspecified_channel_contribution = get_neg(T)
        )
    S = sortperm(get_e(T, P))
    return df[S,:]
end

"""
See also [`state_data`](@ref)

    matrix_data(T::SparseMatrixCSC{Float64, Int64})

Store a reduced matrix elements as a data frame as used by the PAIRINTERACTION software.
"""
function matrix_data(T::SparseMatrixCSC{Float64, Int64})
    m = findnz(T) # non-zero elements
    t = findall(m[1] .<= m[2]) # upper triangle
    df = DataFrame(
        id_initial = m[1][t], 
        id_final = m[2][t], 
        value = m[3][t])
    return df
end

"""
See also [`matrix_data`](@ref)

    tri_to_full(M::DataFrame, S::DataFrame)

Converts a reduced matrix elements data frame storing the upper triangle to a data frame 
storing the full matrix including the wigner phase convention (-1)^(f_final-f_initial).
"""
function tri_to_full(M::DataFrame, S::DataFrame)
    s = size(S, 1)
    f = S.f
    i1 = M.id_initial
    i2 = M.id_final
    val = M.value
    D = sparse(i1, i2, val, s, s)
    for i in axes(D, 1)
        for j in i:s
            d = D[i, j]
            if !iszero(d)
                D[j, i] = d * (-1)^(f[i]-f[j])
            end
        end
    end
    nz = findnz(D)
    df = DataFrame(:id_initial => nz[1], :id_final => nz[2], :val => nz[3])
    return df
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

function exp_q(q::Vector, n::Vector)
    if allequal(q)
        return Float64(q[1])
    else
        return dot(q, n .^ 2)
    end
end

function std_q(q::Vector, n::Vector)
    if allequal(q)
        return 0.0
    else
        e1 = dot(q, n .^ 2)^2
        e2 = dot(q .^ 2, n .^ 2)
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

function exp_S(T::BasisState)
    S = T.S
    n = T.coeff
    return exp_q(S, n)
end

function exp_S(T::BasisArray)
    t = Vector{Float64}(undef, size(T))
    for i in eachindex(t)
        t[i] = exp_S(T.states[i])
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
    n = get_nu(T)
    l = round.(Int, exp_lr(T))
    i0 = findall(iszero, l)
    i1 = findall(iszero, l .- 1)
    i2 = findall(iszero, l .- 2)
    if occursin("Yb", String(P.species))
        n[i0] .+= 4.5
        n[i1] .+= 3.9
        n[i2] .+= 2.8
    else
        n[i0] .+= 3.3
        n[i1] .+= 2.8
        n[i2] .+= 2.5
    end
    return round.(Int, n)
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
    n = T.ai
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
    n = T.ai
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
    n = T.ai
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
    n = T.ai
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
    n = T.ai
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
    n = T.ai
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

function is_mqdt(T::DataBaseArray)
    return repeat([true], size(T))
end
