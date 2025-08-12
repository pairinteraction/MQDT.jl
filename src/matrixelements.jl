# --------------------------------------------------------
# radial
# --------------------------------------------------------
const lru_get_rydberg_state = LRU{Tuple{String,Float64,Int64},Any}(maxsize = 20_000)

"""
    get_rydberg_state_cached(species::String, nu::Float64, l::Int64)

Get a rydberg state from the ryd-numerov package and calculate its wavefunction.
This function is cached using the `LRUCache` package.
"""
function get_rydberg_state_cached(species::String, nu::Float64, l::Int64)
    get!(lru_get_rydberg_state, (species, nu, l)) do
        ryd_numerov = pyimport("ryd_numerov")

        # disable warnings from ryd_numerov for now
        logging = pyimport("logging")
        logging.getLogger("ryd_numerov").setLevel(logging.ERROR)

        energy_au = -0.5 / nu^2  # simple hydrogenic energy with effective principal quantum number nu
        n = ceil(Int, max(nu, l + 1))  # FIXME, n is just used for sanity checks of the wavefunction, not for calculating the wavefunction

        state = ryd_numerov.RydbergState(species, n, l, j = l + 0.5)
        state.set_energy(energy_au)
        state.create_model(potential_type = "coulomb")
        state.create_wavefunction("numerov", sign_convention = "positive_at_outer_bound")
        state
    end
end

const lru_radial_moment = LRU{Tuple{Int,Float64,Float64,Int,Int},Float64}(maxsize = 500_000)

"""
    radial_moment_cached(order::Int, n1, n2, l1, l2)

Calculate the radial matrix element using the python `ryd-numerov` package.
This function is cached using the `LRUCache` package.

# Examples

```julia-repl
MQDT.radial_moment_cached(1, 30, 31, 1, 2)
329.78054480806406
```
"""
function radial_moment_cached(order::Int, n1, n2, l1, l2)
    get!(lru_radial_moment, (order, n1, n2, l1, l2)) do
        state_i = get_rydberg_state_cached("H_textbook", n1, l1)
        state_f = get_rydberg_state_cached("H_textbook", n2, l2)
        radial = state_i.calc_radial_matrix_element(state_f, order, unit = "a.u.")
        pyconvert(Float64, radial)
    end
end

"""
See also [`radial_moment_cached`](@ref)

    function radial_matrix(k::Int, n1, n2, l1, l2)

Evaluates radial_moment over lists.

# Examples

```julia-repl
MQDT.radial_matrix(1, [30, 30], [31, 31], [1, 2], [2, 1])
2×2 Matrix{Float64}:
  329.781  -302.282
 -302.278   276.085
```
"""
function radial_matrix(k::Int, n1, n2, l1, l2)
    R = zeros(length(n1), length(n2))
    for i in eachindex(n1)
        ni = n1[i]
        li = l1[i]
        for j in eachindex(n2)
            nj = n2[j]
            lj = l2[j]
            if abs(li-lj) <= k
                if max(ni, nj) < 25 || abs(ni-nj) < 11 # cut off calculation of matrix elements for F states and higher \ell
                    R[i, j] = radial_moment_cached(k, ni, nj, li, lj)
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
        if iseven(q1.lr+q2.lr+k) &&
           abs(q1.F-q2.F) <= k &&
           abs(q1.lr-q2.lr) <= k &&
           abs(q1.Jr-q2.Jr) <= k
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

const lru_angular_matrix = LRU{Tuple{Int,Channels,Channels},Any}(maxsize = 1_000)

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
        A
    end
end



# ---------------------------------
# multipole moments
# ---------------------------------

"""
See also [`magnetic_dipole_moment`](@ref), [`special_quadrupole_moment`](@ref)

    multipole_moment(k::Int, s1::BasisState, s2::BasisState)

Returns the reduced electric multipole matrix elements of order `k` of multi-channel wave functions.
"""
function multipole_moment(k::Int, s1::BasisState, s2::BasisState)
    p1 = s1.parity
    n1 = s1.nu
    l1 = s1.lr
    a1 = s1.coeff
    k1 = s1.channels
    p2 = s2.parity
    n2 = s2.nu
    l2 = s2.lr
    a2 = s2.coeff
    k2 = s2.channels
    M = 0.0
    if (p1 == p2 && iseven(k)) || (p1 != p2 && isodd(k))
        R = radial_matrix(k, n1, n2, l1, l2)
        Y = angular_matrix_cached(k, k1, k2)
        M = a1' * (Y .* R) * a2
    end
    return M
end

# --------------------------------------------------------
# magnetic field
# --------------------------------------------------------

"""
See also [`multipole_moment`](@ref), [`special_quadrupole_moment`](@ref)

    magnetic_dipole_moment(nd, mp, ic, s1::BasisState, s2::BasisState)

Returns the reduced magentic dipole moment matrix elements of multi-channel wave functions, relevant for the Zeeman Hamiltonian.
"""
function magnetic_dipole_moment(nd, mp, ic, s1::BasisState, s2::BasisState)
    p1 = s1.parity
    n1 = s1.nu
    l1 = s1.lr
    a1 = s1.coeff
    k1 = s1.channels
    p2 = s2.parity
    n2 = s2.nu
    l2 = s2.lr
    a2 = s2.coeff
    k2 = s2.channels
    M = 0.0
    """
    this combination of radial and angular terms is for paramagnetic interaction (Zeeman term)
    """
    if p1 == p2
        R = radial_matrix(0, n1, n2, l1, l2)
        Y = magnetic_matrix_cached(nd, mp, ic, k1, k2)
        M = a1' * (Y .* R) * a2
    end
    return M
end

"""
See also [`multipole_moment`](@ref), [`magnetic_dipole_moment`](@ref)

    special_quadrupole_moment(s1::BasisState, s2::BasisState)

Returns the special electric quadrupole moment matrix elements (i.e. r^2Y_{00}) of multi-channel wave functions, relevant for the diamagnetic Hamiltonian.
"""
function special_quadrupole_moment(s1::BasisState, s2::BasisState)
    p1 = s1.parity
    n1 = s1.nu
    l1 = s1.lr
    a1 = s1.coeff
    k1 = s1.channels
    p2 = s2.parity
    n2 = s2.nu
    l2 = s2.lr
    a2 = s2.coeff
    k2 = s2.channels
    M = 0.0
    """
    this combination of radial and angular terms is for diamagentic interaction
    """
    if p1 == p2
        R = radial_matrix(2, n1, n2, l1, l2)
        Y = angular_matrix_cached(0, k1, k2)
        M = a1' * (Y .* R) * a2
    end
    return M
end

const lru_magnetic_matrix = LRU{Tuple{Any,Any,Any,Channels,Channels},Any}(maxsize = 1_000)

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
        A
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

```julia-repl
MQDT.Λ(5)
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

```julia-repl
MQDT.square_brakets([1, 2, 3])
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

# --------------------------------------------------------
# matrix elements for multipole and magnetic interactions
# --------------------------------------------------------

"""
See also [`multipole_moment`](@ref), [`magnetic_dipole_moment`](@ref), [`special_quadrupole_moment`](@ref)

    matrix_element(k::Int, B::BasisArray)
    matrix_element(B::BasisArray)
    matrix_element(A::Parameters, B::BasisArray)

Iterates `multipole_moment`, `magnetic_dipole_moment`, and `special_quadrupole_moment` over a list of states passed as `BasisArray`
in order to calculate reduced elements.
"""
function matrix_element(k::Int, B::BasisArray)
    """
    for general multipole couplings of order k
    """
    st = B.states
    M = spzeros(length(st), length(st))
    for i in eachindex(st)
        b1 = st[i]
        f1 = b1.f
        for j = i:length(st)
            b2 = st[j]
            f2 = b2.f
            m = multipole_moment(k, b1, b2)
            if !iszero(m)
                M[i, j] = m
                if j != i
                    M[j, i] = (-1)^(f2-f1) * m
                end
            end
        end
    end
    return M
end

function matrix_element(B::BasisArray)
    """
    for diamagnetic couplings (quadratic in the magnetic field)
    """
    st = B.states
    M = spzeros(length(st), length(st))
    for i in eachindex(st)
        b1 = st[i]
        f1 = b1.f
        for j = i:length(st)
            b2 = st[j]
            f2 = b2.f
            m = special_quadrupole_moment(b1, b2)
            if !iszero(m)
                M[i, j] = m
                if j != i
                    M[j, i] = (-1)^(f2-f1) * m
                end
            end
        end
    end
    return M
end

function matrix_element(A::Parameters, B::BasisArray)
    """
    for paramagnetic couplings (linear in the magnetic field)
    """
    nd = A.dipole
    mp = A.mass
    ic = A.spin
    st = B.states
    M = spzeros(size(st, 1), size(st, 1))
    for i in eachindex(st)
        b1 = st[i]
        f1 = b1.f
        for j = i:length(st)
            b2 = st[j]
            f2 = b2.f
            m = magnetic_dipole_moment(nd, mp, ic, b1, b2)
            if !iszero(m)
                M[i, j] = m
                if j != i
                    M[j, i] = (-1)^(f2-f1) * m
                end
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
    for S = 0:1
        for lc in Q.lc
            for lr in Q.lr
                for L = abs(lc-lr):(lc+lr)
                    for J = abs(L-S):(L+S)
                        for F = abs(J-Q.ic):(J+Q.ic)
                            push!(out, lsQuantumNumbers(0.5, S, lc, lr, L, J, F))
                        end
                    end
                end
            end
        end
    end
    return lsChannels(out)
end

function jj_channels(Q::AngularMomenta)
    out = Vector{jjQuantumNumbers}()
    for lc in Q.lc
        for jc = abs(lc-0.5):(lc+0.5)
            for lr in Q.lr
                for jr = abs(lr-0.5):(lr+0.5)
                    for J = abs(jc-jr):(jc+jr)
                        for F = abs(J-Q.ic):(J+Q.ic)
                            push!(out, jjQuantumNumbers(0.5, lc, jc, lr, jr, J, F))
                        end
                    end
                end
            end
        end
    end
    return jjChannels(out)
end

function fj_channels(Q::AngularMomenta)
    out = Vector{fjQuantumNumbers}()
    for lc in Q.lc
        for jc = abs(lc-0.5):(lc+0.5)
            for fc = abs(Q.ic-jc):(Q.ic+jc)
                for lr in Q.lr
                    for jr = abs(lr-0.5):(lr+0.5)
                        for F = abs(fc-jr):(fc+jr)
                            push!(out, fjQuantumNumbers(0.5, lc, jc, fc, lr, jr, F))
                        end
                    end
                end
            end
        end
    end
    return fjChannels(out)
end

function ls_channels(C::lsChannels)
    J = get_J(C)
    values = sort(unique(J))
    channels = Vector{lsChannels}()
    for j in values
        c = J .== j
        push!(channels, lsChannels(C.i[c]))
    end
    return values, channels
end

function jj_channels(C::jjChannels)
    J = get_J(C)
    values = sort(unique(J))
    channels = Vector{jjChannels}()
    for j in values
        c = J .== j
        push!(channels, jjChannels(C.i[c]))
    end
    return values, channels
end

function fj_channels(C::fjChannels)
    F = get_F(C)
    values = sort(unique(F))
    channels = Vector{fjChannels}()
    for f in values
        c = F .== f
        push!(channels, fjChannels(C.i[c]))
    end
    return values, channels
end

function ls_to_jj(q1::lsQuantumNumbers, q2::jjQuantumNumbers)
    res = 0.0
    if (q1.sc, q1.lc, q1.lr, q1.J) == (q2.sc, q2.lc, q2.lr, q2.J)
        sq = square_brakets([q1.S, q1.L, q2.Jc, q2.Jr])
        qn = Vector{Int64}(2*[q1.sc, 0.5, q1.S, q1.lc, q1.lr, q1.L, q2.Jc, q2.Jr, q1.J])
        res = sq * f9j(qn...)
    end
    return res
end

function jj_to_fj(q1::jjQuantumNumbers, q2::fjQuantumNumbers, ic)
    res = 0.0
    if (q1.sc, q1.lc, q1.lr, q1.Jc, q1.Jr) == (q2.sc, q2.lc, q2.lr, q2.Jc, q2.Jr)
        Λ = q1.Jr + q2.Fc - ic - q1.J
        sq = square_brakets([q1.J, q2.Fc])
        qn = Vector{Int64}(2*[q1.Jr, q1.Jc, q1.J, ic, q2.F, q2.Fc])
        res = (-1.0)^Λ * sq * f6j(qn...)
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
