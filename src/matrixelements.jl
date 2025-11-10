# --------------------------------------------------------
# radial
# --------------------------------------------------------
const lru_get_rydberg_state = LRU{Tuple{Symbol,Float64,Int64},Any}(; maxsize=20_000)

"""
    get_radial_state_cached(species::Symbol, nu::Float64, l::Int64)

Get a rydberg state from the rydstate package and calculate its wavefunction.
This function is cached using the `LRUCache` package.
"""
function get_radial_state_cached(species::Symbol, nu::Float64, l::Int64)
    get!(lru_get_rydberg_state, (species, nu, l)) do
        rydstate = pyimport("rydstate")

        state = rydstate.radial.RadialState(String(species); nu=nu, l_r=l)
        state.create_model(; potential_type="model_potential_fei_2009")
        state.create_wavefunction("numerov"; sign_convention="positive_at_outer_bound")
        return state
    end
end

const lru_radial_moment = LRU{Tuple{Symbol,Float64,Int,Float64,Int,Int},Float64}(; maxsize=500_000)

"""
    radial_moment_cached(species::Symbol, nu1::Float64, l1::Int, nu2::Float64, l2::Int, order::Int)

Calculate the radial matrix element using the python `rydstate` package.
This function is cached using the `LRUCache` package.

# Examples

```jldoctest
radial_moment = MQDT.radial_moment_cached(:Yb174, 70.2667, 1, 71.2638, 2, 1)
isapprox(radial_moment, 1712.16; rtol=1e-5)

# output

true
```
"""
function radial_moment_cached(species::Symbol, nu1::Float64, l1::Int, nu2::Float64, l2::Int, order::Int)
    get!(lru_radial_moment, (species, nu1, l1, nu2, l2, order)) do
        state_i = get_radial_state_cached(species, nu1, l1)
        state_f = get_radial_state_cached(species, nu2, l2)
        radial = state_f.calc_matrix_element(state_i, order; unit="a.u.")
        return pyconvert(Float64, radial)
    end
end

"""
See also [`radial_moment_cached`](@ref)

    function radial_matrix(s1::BasisState, s2::BasisState, order::Int)

Evaluates radial_moment over two BasisStates, returns a matrix,
where each entry corresponds to a pair of sqdt states from the two BasisStates.
"""
function radial_matrix(s1::BasisState, s2::BasisState, order::Int)
    @assert s1.species == s2.species "Species mismatch: s1.species ($(s1.species)) != s2.species ($(s2.species))"
    R = zeros(length(s1.nu), length(s2.nu))

    for i in eachindex(s1.nu)
        nui = s1.nu[i]
        li = s1.lr[i]
        for j in eachindex(s2.nu)
            nuj = s2.nu[j]
            lj = s2.lr[j]

            if abs(li-lj) <= order
                if max(nui, nuj) < 25 || abs(nui-nuj) < 11 # cut off calculation of matrix elements for F states and higher \ell
                    R[i, j] = radial_moment_cached(s1.species, nui, li, nuj, lj, order)
                end
            end
        end
    end

    return R
end

# ---------------------------------
# angular
# ---------------------------------

"""
See also [`angular_matrix_cached`](@ref)

    angular_moment(k, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    angular_moment(k, q1::jjQuantumNumbers, q2::jjQuantumNumbers)
    angular_moment(k, q1::lsQuantumNumbers, q2::lsQuantumNumbers)

Returns the angular matrix elements (i.e. an analytically evaluated integral of spherical harmonics).
Formula is found in Robicheaux2018 Eq. 20 and in Vaillant2014 Eq. C3
"""
function angular_moment(k, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    a = 0.0
    if (q1.sc, q1.lc, q1.Jc, q1.Fc) == (q2.sc, q2.lc, q2.Jc, q2.Fc)
        if iseven(q1.lr+q2.lr+k) && abs(q1.F-q2.F) <= k && abs(q1.lr-q2.lr) <= k && abs(q1.Jr-q2.Jr) <= k
            Λ = q1.F + q1.Fc + q1.Jr + q2.Jr + q1.lr + q2.lr + 0.5
            sq = square_brakets([q1.F, q2.F, q1.Jr, q2.Jr, q1.lr, q2.lr])
            qn1 = Vector{Int64}(2*[q1.lr, k, q2.lr, 0, 0, 0])
            qn2 = Vector{Int64}(2*[q1.Jr, q1.F, q1.Fc, q2.F, q2.Jr, k])
            qn3 = Vector{Int64}(2*[q1.lr, q1.Jr, 0.5, q2.Jr, q2.lr, k])
            a = (-1.0)^Λ * sq * f3j(qn1...) * f6j(qn2...) * f6j(qn3...)
        end
    end
    return a
end

function angular_moment(k, q1::jjQuantumNumbers, q2::jjQuantumNumbers)
    c1 = fj_quantum_numbers(q1)
    c2 = fj_quantum_numbers(q2)
    return angular_moment(k, c1, c2)
end

# function angular_moment(k, q1::jjQuantumNumbers, q2::jjQuantumNumbers)
#     a = 0.
#     if (q1.sc, q1.lc, q1.Jc) == (q2.sc, q2.lc, q2.Jc)
#         if iseven(q1.lr+q2.lr+k) && abs(q1.J-q2.J) <= k && abs(q1.lr-q2.lr) <= k && abs(q1.Jr-q2.Jr) <= k
#             Λ = q1.J + q1.Jc + q1.Jr + q2.Jr + q1.lr + q2.lr + 0.5
#             sq = square_brakets([q1.J, q2.J, q1.Jr, q2.Jr, q1.lr, q2.lr])
#             qn1 = Vector{Int64}(2*[q1.lr, k, q2.lr, 0, 0, 0])
#             qn2 = Vector{Int64}(2*[q1.Jr, q1.J, q1.Jc, q2.J, q2.Jr, k])
#             qn3 = Vector{Int64}(2*[q1.lr, q1.Jr, 0.5, q2.Jr, q2.lr, k])
#             a = (-1.)^Λ * sq * f3j(qn1...) * f6j(qn2...) * f6j(qn3...)
#         end
#     end
#     return a
# end

function angular_moment(k, q1::lsQuantumNumbers, q2::lsQuantumNumbers)
    a = 0.0
    if (q1.lc, q1.S) == (q2.lc, q2.S)
        if iseven(q1.lr+q2.lr+k) && abs(q1.J-q2.J) <= k && abs(q1.lr-q2.lr) <= k
            Λ = q1.lc + q1.S
            sq = square_brakets([q1.lr, q2.lr, q1.L, q2.L, q1.J, q2.J])
            qn1 = Vector{Int64}(2*[q1.lr, k, q2.lr, 0, 0, 0])
            qn2 = Vector{Int64}(2*[q1.J, k, q2.J, q1.L, q1.S, q2.L])
            qn3 = Vector{Int64}(2*[q1.L, k, q2.L, q1.lr, q1.lc, q2.lr])
            a = (-1.0)^Λ * sq * f3j(qn1...) * f6j(qn2...) * f6j(qn3...)
        end
    end
    return a
end

const lru_angular_matrix = LRU{Tuple{Int,Channels,Channels},Any}(; maxsize=1_000)

"""
See also [`angular_moment`](@ref)

    angular_matrix_cached(k::Int, k1::Channels, k2::Channels)

Iterates `angular_moment` over Channels, returns a matrix.
This function is cached using the `LRUCache` package.
"""
function angular_matrix_cached(k::Int, k1::Channels, k2::Channels)
    get!(lru_angular_matrix, (k, k1, k2)) do
        c1 = k1.i
        c2 = k2.i
        A = zeros(length(c1), length(c2))
        for i in eachindex(c1)
            for j in eachindex(c2)
                A[i, j] = angular_moment(k, c1[i], c2[j])
            end
        end
        return A
    end
end

# --------------------------------------------------------
# magnetic field
# --------------------------------------------------------

const lru_magnetic_matrix = LRU{Tuple{Any,Any,Any,Channels,Channels},Any}(; maxsize=1_000)

"""
See also [`magneton`](@ref)

    magnetic_matrix_cached(nd, mp, ic, k1::Channels, k2::Channels)

Iterates `magneton` over Channels, returns a matrix.
This function is cached using the `LRUCache` package.
"""
function magnetic_matrix_cached(nd, mp, ic, k1::Channels, k2::Channels)
    get!(lru_magnetic_matrix, (nd, mp, ic, k1, k2)) do
        c1 = k1.i
        c2 = k2.i
        A = zeros(length(c1), length(c2))
        for i in eachindex(c1)
            for j in eachindex(c2)
                A[i, j] = magneton(nd, mp, ic, c1[i], c2[j])
            end
        end
        return A
    end
end

"""
    magneton(nd, mp, ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    magneton(nd, mp, ic, q1::jjQuantumNumbers, q2::jjQuantumNumbers)

Returns the reduced magnetic dipole moment (excluding radial contributions from the Rydberg electron)
"""
function magneton(nd, mp, ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    μ_B = 0.5
    μ_N = nd/mp
    g_s = 2.0023192
    m_lr = element_lr(q1, q2)
    m_sr = element_sr(q1, q2)
    m_ic = element_ic(ic, q1, q2)
    m_lc = element_lc(ic, q1, q2)
    m_sc = element_sc(ic, q1, q2)
    return - μ_B*(m_lr + m_lc + g_s*(m_sr + m_sc)) + μ_N*m_ic
end

function magneton(nd, mp, ic, q1::jjQuantumNumbers, q2::jjQuantumNumbers)
    c1 = fj_quantum_numbers(q1)
    c2 = fj_quantum_numbers(q2)
    return magneton(nd, mp, ic, c1, c2)
end

"""
    element_lr(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of Rydberg orbital anugular momentum (Robicheaux2018 Eq. 24)
"""
function element_lr(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(q1.lr) * G1(q1, q2) * G2(q1, q2)
    else
        return 0.0
    end
end

"""
    element_sr(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of Rydberg spin (Robicheaux2018 Eq. 24)
"""
function element_sr(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(0.5) * G1(q1, q2) * G3(q1, q2)
    else
        return 0.0
    end
end

"""
    element_ic(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of core nuclear spin (Robicheaux2018 Eq. 24)
"""
function element_ic(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(ic) * G4(q1, q2) * G5(q1, q2, ic)
    else
        return 0.0
    end
end

"""
    element_lc(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of core orbital anugular momentum (Robicheaux2018 Eq. 24)
"""
function element_lc(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(q1.lc) * G4(q1, q2) * G6(q1, q2, ic) * G7(q1, q2)
    else
        return 0.0
    end
end

"""
    element_sc(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of core spin (Robicheaux2018 Eq. 24)
"""
function element_sc(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(0.5) * G4(q1, q2) * G6(q1, q2, ic) * G8(q1, q2)
    else
        return 0.0
    end
end

"""
    Λ(x) = sqrt((2x+1) * (x+1) * x)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)

# Examples

```jldoctest
MQDT.Λ(5)

# output

18.16590212458495
```
"""
function Λ(x)
    return sqrt((2x+1) * (x+1) * x)
end

"""
    square_brakets(a::Vector) = sqrt(prod(2a.+1))

Auxiliary function to calculate angular overlaps (Robicheaux2018 Eq. 11)

# Examples

```jldoctest
MQDT.square_brakets([1, 2, 3])

# output

10.246950765959598
```
"""
function square_brakets(a::Vector)
    return sqrt(prod(2a .+ 1))
end

"""
    G1(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G1(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.0
    if q1.Jc == q2.Jc && q1.Fc == q2.Fc
        Λ = q1.Fc + q2.Jr + q1.F + 1
        sq = square_brakets([q1.F, q2.F])
        qn = Vector{Int64}(2*[q1.Jr, q1.F, q1.Fc, q2.F, q2.Jr, 1])
        g = (-1.0)^Λ * sq * f6j(qn...)
    end
    return g
end

"""
    G2(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G2(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.0
    if q1.lr == q2.lr
        Λ = 0.5 + q1.lr + q1.Jr + 1
        sq = square_brakets([q1.Jr, q2.Jr])
        qn = Vector{Int64}(2*[q1.lr, q1.Jr, 0.5, q2.Jr, q1.lr, 1])
        g = (-1.0)^Λ * sq * f6j(qn...)
    end
    return g
end

"""
    G3(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G3(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.0
    if q1.lr == q2.lr
        Λ = 0.5 + q1.lr + q2.Jr + 1
        sq = square_brakets([q1.Jr, q2.Jr])
        qn = Vector{Int64}(2*[0.5, q1.Jr, q1.lr, q2.Jr, 0.5, 1])
        g = (-1.0)^Λ * sq * f6j(qn...)
    end
    return g
end

"""
    G4(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G4(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.0
    if q1.Jr == q2.Jr
        Λ = q1.Fc + q1.Jr + q2.F + 1
        sq = square_brakets([q1.F, q2.F])
        qn = Vector{Int64}(2*[q1.Fc, q1.F, q1.Jr, q2.F, q2.Fc, 1])
        g = (-1.0)^Λ * sq * f6j(qn...)
    end
    return g
end

"""
    G5(q1::fjQuantumNumbers, q2::fjQuantumNumbers, ic)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G5(q1::fjQuantumNumbers, q2::fjQuantumNumbers, ic)
    g = 0.0
    if q1.Jc == q2.Jc
        Λ = q1.Jc + ic + q1.Fc + 1
        sq = square_brakets([q1.Fc, q2.Fc])
        qn = Vector{Int64}(2*[ic, q1.Fc, q1.Jc, q2.Fc, ic, 1])
        g = (-1.0)^Λ * sq * f6j(qn...)
    end
    return g
end

"""
    G6(q1::fjQuantumNumbers, q2::fjQuantumNumbers, ic)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G6(q1::fjQuantumNumbers, q2::fjQuantumNumbers, ic)
    Λ = q1.Jc + ic + q2.Fc + 1
    sq = square_brakets([q1.Fc, q2.Fc])
    qn = Vector{Int64}(2*[q1.Jc, q1.Fc, ic, q2.Fc, q2.Jc, 1])
    g = (-1.0)^Λ * sq * f6j(qn...)
    return g
end

"""
    G7(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G7(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.0
    if q1.sc == q2.sc && q1.lc == q2.lc
        Λ = q1.sc + q1.lc + q1.Jc + 1
        sq = square_brakets([q1.Jc, q2.Jc])
        qn = Vector{Int64}(2*[q1.lc, q1.Jc, q1.sc, q2.Jc, q1.lc, 1])
        g = (-1.0)^Λ * sq * f6j(qn...)
    end
    return g
end

"""
    G8(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G8(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.0
    if q1.sc == q2.sc && q1.lc == q2.lc
        Λ = q1.sc + q1.lc + q2.Jc + 1
        sq = square_brakets([q1.Jc, q2.Jc])
        qn = Vector{Int64}(2*[q1.sc, q1.Jc, q1.lc, q2.Jc, q1.sc, 1])
        g = (-1.0)^Λ * sq * f6j(qn...)
    end
    return g
end

# ---------------------------------
# multipole moments
# ---------------------------------

"""
    multipole_moment(s1::BasisState, s2::BasisState, P::Paramaters)

Returns the reduced matrix elements for electric dipole, electric quadrupole, magnetic dipole, and electric dipole squared.
"""
function multipole_moments(s1::BasisState, s2::BasisState, P::Parameters)
    M = zeros(4)
    if s1.parity != s2.parity
        # electric dipole
        R = radial_matrix(s1, s2, 1)
        Y = angular_matrix_cached(1, s1.channels, s2.channels)
        M[1] = s1.coeff' * (Y .* R) * s2.coeff
    end
    if s1.parity == s2.parity
        # electric quadrupole
        R = radial_matrix(s1, s2, 2)
        Y = angular_matrix_cached(2, s1.channels, s2.channels)
        M[2] = s1.coeff' * (Y .* R) * s2.coeff
        # electric dipole squared
        Y = angular_matrix_cached(0, s1.channels, s2.channels)
        M[3] = s1.coeff' * (Y .* R) * s2.coeff
        # magnetic dipole
        R = radial_matrix(s1, s2, 0)
        Y = magnetic_matrix_cached(P.dipole, P.mass, P.spin, s1.channels, s2.channels)
        M[4] = s1.coeff' * (Y .* R) * s2.coeff
    end
    return M
end

function matrix_elements(B::BasisArray, P::Parameters)
    M = Dict(
        "dipole" => Tuple{Int64,Int64,Float64}[],
        "quadrupole" => Tuple{Int64,Int64,Float64}[],
        "diamagnetic" => Tuple{Int64,Int64,Float64}[],
        "paramagnetic" => Tuple{Int64,Int64,Float64}[],
    )
    st = B.states
    for i in eachindex(st)
        b1 = st[i]
        for j in i:length(st)
            b2 = st[j]
            m = multipole_moments(b1, b2, P)
            if !iszero(m[1])
                push!(M["dipole"], (i, j, m[1]))
            end
            if !iszero(m[2])
                push!(M["quadrupole"], (i, j, m[2]))
            end
            if !iszero(m[3])
                push!(M["diamagnetic"], (i, j, m[3]))
            end
            if !iszero(m[4])
                push!(M["paramagnetic"], (i, j, m[4]))
            end
        end
    end
    return M
end

# --------------------------------------------------------
# frame transformations
# --------------------------------------------------------
#   Additional code to span spin frames
#   and calculate transformations
# --------------------------------------------------------

struct AngularMomenta
    lc::Vector{Number}
    lr::Vector{Number}
    ic::Number
end

function ls_channels(Q::AngularMomenta)
    out = Vector{lsQuantumNumbers}()
    for S in 0:1
        for lc in Q.lc
            for lr in Q.lr
                for L in abs(lc - lr):(lc + lr)
                    for J in abs(L - S):(L + S)
                        for F in abs(J - Q.ic):(J + Q.ic)
                            push!(out, lsQuantumNumbers(0.5, S, lc, lr, L, J, F))
                        end
                    end
                end
            end
        end
    end
    C = lsChannels(out)
    F = get_F(C)
    channels = Vector{lsChannels}()
    for j in sort(unique(F))
        c = F .== f
        push!(channels, lsChannels(C.i[c]))
    end
    return channels
end

function jj_channels(Q::AngularMomenta)
    out = Vector{jjQuantumNumbers}()
    for lc in Q.lc
        for jc in abs(lc - 0.5):(lc + 0.5)
            for lr in Q.lr
                for jr in abs(lr - 0.5):(lr + 0.5)
                    for J in abs(jc - jr):(jc + jr)
                        for F in abs(J - Q.ic):(J + Q.ic)
                            push!(out, jjQuantumNumbers(0.5, lc, jc, lr, jr, J, F))
                        end
                    end
                end
            end
        end
    end
    C = jjChannels(out)
    F = get_F(C)
    channels = Vector{jjChannels}()
    for f in sort(unique(F))
        c = F .== f
        push!(channels, jjChannels(C.i[c]))
    end
    return channels
end

function fj_channels(Q::AngularMomenta)
    out = Vector{fjQuantumNumbers}()
    for lc in Q.lc
        for jc in abs(lc - 0.5):(lc + 0.5)
            for fc in abs(Q.ic - jc):(Q.ic + jc)
                for lr in Q.lr
                    for jr in abs(lr - 0.5):(lr + 0.5)
                        for F in abs(fc - jr):(fc + jr)
                            push!(out, fjQuantumNumbers(0.5, lc, jc, fc, lr, jr, F))
                        end
                    end
                end
            end
        end
    end
    C = fjChannels(out)
    F = get_F(C)
    channels = Vector{fjChannels}()
    for f in sort(unique(F))
        c = F .== f
        push!(channels, fjChannels(C.i[c]))
    end
    return channels
end

function ls_to_jj(q1::lsQuantumNumbers, q2::jjQuantumNumbers)
    res = 0.0
    if (q1.sc, q1.lc, q1.lr, q1.J) == (q2.sc, q2.lc, q2.lr, q2.J)
        sq = square_brakets([q1.S, q1.L, q2.Jc, q2.Jr])
        qn = Vector{Int64}(2*[q1.sc, 0.5, q1.S, q1.lc, q1.lr, q1.L, q2.Jc, q2.Jr, q1.J])
        res = sq*f9j(qn...)
        if isapprox(round(2res), 2res; rtol=1e-14)
            res = round(res; digits=14)
        end
    end
    return res
end

function jj_to_fj(q1::jjQuantumNumbers, q2::fjQuantumNumbers, ic)
    res = 0.0
    if (q1.sc, q1.lc, q1.lr, q1.Jc, q1.Jr) == (q2.sc, q2.lc, q2.lr, q2.Jc, q2.Jr)
        Λ_io = q1.Jr + q2.Fc + ic + q1.J
        sq = square_brakets([q1.J, q2.Fc])
        qn = Vector{Int64}(2*[q1.Jr, q1.Jc, q1.J, ic, q2.F, q2.Fc])
        res = (-1.0)^Λ_io*sq*f6j(qn...)
        if isapprox(round(2res), 2res; rtol=1e-14)
            res = round(res; digits=14)
        end
    end
    return res
end

function matrix_ls_to_jj(in::lsChannels, out::jjChannels)
    t = Matrix{Float64}(undef, size(in), size(out))
    for i in axes(t, 1)
        for j in axes(t, 2)
            t[i, j] = ls_to_jj(in.i[i], out.i[j])
        end
    end
    return t
end

function matrix_jj_to_fj(in::jjChannels, out::fjChannels, ic::Number)
    t = Matrix{Float64}(undef, size(in), size(out))
    for i in axes(t, 1)
        for j in axes(t, 2)
            t[i, j] = jj_to_fj(in.i[i], out.i[j], ic)
        end
    end
    return t
end

function jj_frame_transformation(l::Integer)
    am = AngularMomenta([0], [l], 0)
    ls = ls_channels(am)
    jj = jj_channels(am)
    t = Vector{Matrix{Float64}}(undef, length(ls))
    for k in eachindex(ls)
        t[k] = matrix_ls_to_jj(ls[k], jj[k])'
    end
    return ls, jj, t
end

function fj_frame_transformation(l::Integer, ic::Number)
    am = AngularMomenta([0], [l], ic)
    ls = ls_channels(am)
    jj = jj_channels(am)
    fj = fj_channels(am)
    t = Vector{Matrix{Float64}}(undef, length(ls))
    for k in eachindex(ls)
        t1 = matrix_ls_to_jj(ls[k], jj[k])
        t2 = matrix_jj_to_fj(jj[k], fj[k], ic)
        t[k] = t2' * t1'
    end
    return ls, fj, t
end
