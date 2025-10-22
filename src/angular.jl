
function calc_exp_qn(T::BasisState, qn::String)
    state = get_state(T)
    value = state.calc_exp_qn(qn)
    return pyconvert(Float64, value)
end

function calc_exp_qn(T::BasisArray, qn::String)
    exp_list = Vector{Float64}(undef, size(T))
    for (i, state) in enumerate(T.states)
        exp_list[i] = calc_exp_qn(state, qn)
    end
    return exp_list
end

function calc_std_qn(T::BasisState, qn::String)
    state = get_state(T)
    value = state.calc_std_qn(qn)
    return pyconvert(Float64, value)
end

function calc_std_qn(T::BasisArray, qn::String)
    std_list = Vector{Float64}(undef, size(T))
    for (i, state) in enumerate(T.states)
        std_list[i] = calc_std_qn(state, qn)
    end
    return std_list
end

function get_state(T::BasisState)
    ryd_numerov = pyimport("ryd_numerov")
    channels = T.channels.i
    kets = Vector{Any}(undef, size(channels))
    for i in eachindex(channels)
        qn = channels[i]
        kets[i] = get_ket(qn, T.species)
    end
    state = ryd_numerov.angular.AngularState(T.coeff, kets; warn_if_not_normalized=false)
    return state
end

function get_ket(qn::lsQuantumNumbers, species::Symbol)
    ryd_numerov = pyimport("ryd_numerov")
    ket = ryd_numerov.angular.AngularKetLS(;
        s_c=qn.sc,
        s_tot=qn.S,
        l_c=qn.lc,
        l_r=qn.lr,
        l_tot=qn.L,
        j_tot=qn.J,
        f_tot=qn.F,
        species=String(species),
    )
    return ket
end

function get_ket(qn::jjQuantumNumbers, species::Symbol)
    ryd_numerov = pyimport("ryd_numerov")
    ket = ryd_numerov.angular.AngularKetJJ(;
        s_c=qn.sc,
        l_c=qn.lc,
        l_r=qn.lr,
        j_c=qn.Jc,
        j_r=qn.Jr,
        j_tot=qn.J,
        f_tot=qn.F,
        species=String(species),
    )
    return ket
end

function get_ket(qn::fjQuantumNumbers, species::Symbol)
    ryd_numerov = pyimport("ryd_numerov")
    ket = ryd_numerov.angular.AngularKetFJ(;
        s_c=qn.sc,
        l_c=qn.lc,
        l_r=qn.lr,
        j_c=qn.Jc,
        f_c=qn.Fc,
        j_r=qn.Jr,
        f_tot=qn.F,
        species=String(species),
    )
    return ket
end
