# --------------------------------------------------------
# radial
# --------------------------------------------------------

"""
    radialwavefunction(r::Number, n::Int, l::Int) = r* 2/n^2* sqrt(sf_gamma(n-l)/sf_gamma(n+l+1))* exp(-r/n)* (2r/n)^l* sf_laguerre_n(n-l-1, 2l+1, 2r/n)
    radialwavefunction(r::Number, n::Float64, l::Int) = 1/n* 1/sqrt(sf_gamma(n+l+1)* sf_gamma(n-l))* exp(-r/n)* (2r/n)^(l+1)* sf_hyperg_U(l+1-n, 2(l+1), 2r/n)

Evaluate the radial wave function for hydrogenic states using Laguerre polynomials and for quantum defect states using Whitaker functions.
Calls function from the Gnu Scientific Library package `GSL`.

# Examples

```julia-repl
MQDT.radialwavefunction(1000, 40, 1)
-0.014041583633649013

MQDT.radialwavefunction(1000, 40., 1)
-0.014041583633649006
```
"""
function radialwavefunction(r::Number, n::Int, l::Int)
    r * 2/n^2 * sqrt(sf_gamma(n-l)/sf_gamma(n+l+1)) * exp(-r/n) * (2r/n)^l * sf_laguerre_n(n-l-1, 2l+1, 2r/n)
end

function radialwavefunction(r::Number, n::Float64, l::Int)
    n = round(n, digits=8)
    1/n * 1/sqrt(sf_gamma(n+l+1) * sf_gamma(n-l)) * exp(-r/n) * (2r/n)^(l+1) * sf_hyperg_U(l+1-n, 2(l+1), 2r/n)
end

"""
    integrator_start(n, l)

Heuristic for determination of the inner bound for radial integration depending on angular momentum.
"""
function integrator_start(n, l)
    if l < 1
        return 1e-3
    elseif l < 2
        return 1/4
    elseif l < 4
        return (2l-1)/2
    else
        return n^2 - n*sqrt(n^2 - l^2)
    end
end

"""
    integrate(f::Function, g::Function, order::Int, start::Number, stop::Number)

Radial integration using the package `QuadGK`.
"""
function integrate(f::Function, g::Function, order::Int, start::Number, stop::Number)
    func(x) = f(x) * g(x) * x^order
    return quadgk(func, start, stop, rtol=1e-8)
end

"""
See also [`radial_moment`](@ref)

    radial_overlap(n1, n2, l1, l2)

Analytic overlap of two radial wave functions.

# Examples

```julia-repl
MQDT.radial_overlap(30, 31, 1, 2)
0.999865618490289
```
"""
function radial_overlap(n1, n2, l1, l2)
    b1 = beta(n1, l1)
    b2 = beta(n2, l2)
    if n1 == n2 && l1 == l2
        return 1.
    elseif b1 != b2
        return 2sqrt(n1*n2)/(n1+n2) * sin(b1-b2)/(b1-b2)
    else
        return 2sqrt(n1*n2)/(n1+n2)
    end
end

"""
See also [`radial_overlap`](@ref)

    radial_moment(order::Int, n1, n2, l1, l2)

Numeric integral of two radial wave functions.
Returns the result and the accuracy as a tuple.

# Examples

```julia-repl
MQDT.radial_moment(1, 30, 31, 1, 2)
(329.78054480806827, 2.2900528483368755e-6)
```
"""
@memoize function radial_moment(order::Int, n1, n2, l1, l2)
    w1(r) = radialwavefunction(r, n1, l1)
    w2(r) = radialwavefunction(r, n2, l2)
    if l1 > l2
        s1 = integrator_start(n1, l1)
    else
        s1 = integrator_start(n2, l2)
    end
    s2 = 5min(n1, n2)^2
    return integrate(w1, w2, order, s1, s2)
end

"""
See also [`radial_overlap`](@ref), [`radial_moment`](@ref), 

    radial_integral(order::Int, n1, n2, l1, l2)

Combines the functions `radial_overlap` and `radial_moment`.

# Examples

```julia-repl
MQDT.radial_integral(1, 30, 31, 1, 2)
(329.78054480806827, 2.2900528483368755e-6)
```
"""
function radial_integral(order::Int, n1, n2, l1, l2)
    if iszero(order)
        return (radial_overlap(n1, n2, l1, l2), 0.)
    else
        return radial_moment(order, n1, n2, l1, l2)
    end
end

"""
See also [`radial_integral`](@ref)

    function radial_matrix(k::Int, n1, n2, l1, l2)

Evaluates `radial_integral` over lists.

# Examples

```julia-repl
MQDT.radial_matrix(1, [30, 30], [31, 31], [1, 2], [2, 1])
2×2 Matrix{Float64}:
  329.781  -302.282
 -302.278   276.085
```
"""
@memoize function radial_matrix(k::Int, n1, n2, l1, l2)
    R = zeros(length(n1), length(n2))
    for i in eachindex(n1)
        for j in eachindex(n2)
            R[i,j] = radial_integral(k, n1[i], n2[j], l1[i], l2[j])[1]
        end
    end
    return R
end

# ---------------------------------
# angular
# ---------------------------------

"""
See also [`angular_matrix`](@ref)

    angular_moment(k, q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the angular matrix elements (i.e. an analytically evaluated integral of spherical harmonics).
Formula is found in Robicheaux2018 Eq. 20
"""
function angular_moment(k, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    a = 0.
    if (q1.sc, q1.lc, q1.Jc, q1.Fc) == (q2.sc, q2.lc, q2.Jc, q2.Fc)
        q1.sr == q2.sr && iseven(q1.lr+q2.lr+k) && abs(q1.F-q2.F) <= k && abs(q1.lr-q2.lr) <= k && abs(q1.Jr-q2.Jr) <= k
        Λ = q1.F + q1.Fc + q1.Jr + q2.Jr + q1.lr + q2.lr + q1.sr
        sq = square_brakets([q1.F, q2.F, q1.Jr, q2.Jr, q1.lr, q2.lr])
        qn1 = Vector{Int64}(2*[q1.lr, k, q2.lr, 0, 0, 0])
        qn2 = Vector{Int64}(2*[q1.Jr, q1.F, q1.Fc, q2.F, q2.Jr, k])
        qn3 = Vector{Int64}(2*[q1.lr, q1.Jr, q1.sr, q2.Jr, q2.lr, k])
        a = (-1.)^Λ * sq * sf_coupling_3j(qn1...) * sf_coupling_6j(qn2...) * sf_coupling_6j(qn3...)
    end
    return a
end

"""
See also [`fjChannels`](@ref)

    fj_channels(C::Channels)
    fj_channels(C::jjChannels)

Type checking function for `Channels`. 
If the input is an `fjChannels`, returns the input. 
If the input is a `jjChannels`, constructs and returns the corresponding `fjChannel`.
Does nothing for other inputs.
"""
function fj_channels(C::Channels)
    if typeof(C) == fjChannels
        return C
    elseif typeof(C) == jjChannels
        return fj_channels(C)
    end
end

function fj_channels(C::jjChannels)
    return fjChannels(C.sc, C.lc, C.Jc, C.Jc, C.sr, C.lr, C.Jr, C.J)
end

"""
See also [`fjQuantumNumbers`](@ref)

    fj_quantum_numbers(C::Vector)

Creates a `fjQuantumNumbers` type required for evaluating functions that calculate matrix elements.
Note that `fjQuantumNumbers` requires a list of arguments, while `fj_quantum_numbers` supports all arguments being passed as a single vector.

# Examples

```julia-repl
MQDT.fj_quantum_numbers([0.5, 0, 0.5, 1., 0.5, 0, 0.5, 1.5])
Main.MQDT.fjQuantumNumbers(0.5, 0, 0.5, 1.0, 0.5, 0, 0.5, 1.5)
```
"""
function fj_quantum_numbers(C::Vector)
    return fjQuantumNumbers(C...)
end

"""
See also [`angular_moment`](@ref)

    angular_matrix(k::Int, k1::Channels, k2::Channels)

Iterates `angular_moment` over Channels, returns a matrix.
This function is 'memoized' using the `Memoize` package.
"""
@memoize function angular_matrix(k::Int, k1::Channels, k2::Channels)
    c1 = fj_channels(k1)
    c2 = fj_channels(k2)
    v1 = cat(c1)
    v2 = cat(c2)
    A = zeros(size(v1, 2), size(v2, 2))
    for i in axes(v1, 2)
        q1 = fj_quantum_numbers(v1[:,i])
        for j in axes(v2, 2)
            q2 = fj_quantum_numbers(v2[:,j])
            A[i,j] = angular_moment(k, q1, q2)
        end
    end
    return A
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
    M = 0.
    if p1 == p2 && iseven(k) || p1 != p2 && isodd(k)
        R = radial_matrix(k, n1, n2, l1, l2)
        Y = angular_matrix(k, k1, k2)
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
    M = 0.
    """
    this combination of radial and angular terms is for paramagnetic interaction (Zeeman term)
    """
    if p1 == p2
        R = radial_matrix(0, n1, n2, l1, l2)
        Y = magnetic_matrix(nd, mp, ic, k1, k2)
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
    M = 0.
    """
    this combination of radial and angular terms is for diamagentic interaction
    """
    if p1 == p2
        R = radial_matrix(2, n1, n2, l1, l2)
        Y = angular_matrix(0, k1, k2)
        M = a1' * (Y .* R) * a2
    end 
    return M
end

"""
See also [`magneton`](@ref)

    magnetic_matrix(nd, mp, ic, k1::Channels, k2::Channels)

Iterates `magneton` over Channels, returns a matrix.
This function is 'memoized' using the `Memoize` package.
"""
@memoize function magnetic_matrix(nd, mp, ic, k1::Channels, k2::Channels)
    c1 = fj_channels(k1)
    c2 = fj_channels(k2)
    v1 = cat(c1)
    v2 = cat(c2)
    A = zeros(size(v1, 2), size(v2, 2))
    for i in axes(v1, 2)
        q1 = fj_quantum_numbers(v1[:,i])
        for j in axes(v2, 2)
            q2 = fj_quantum_numbers(v2[:,j])
            A[i,j] = magneton(nd, mp, ic, q1, q2)
        end
    end
    return A
end

"""
    magneton(nd, mp, ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)

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

"""
    element_lr(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of Rydberg orbital anugular momentum (Robicheaux2018 Eq. 24)
"""
function element_lr(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.sr == q2.sr && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(q1.lr) * G1(q1, q2) * G2(q1, q2)
    else
        return 0.
    end
end

"""
    element_sr(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of Rydberg spin (Robicheaux2018 Eq. 24)
"""
function element_sr(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.sr == q2.sr && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(0.5) * G1(q1, q2) * G3(q1, q2)
    else
        return 0.
    end
end

"""
    element_ic(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of core nuclear spin (Robicheaux2018 Eq. 24)
"""
function element_ic(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.sr == q2.sr && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(ic) * G4(q1, q2) * G5(q1, q2, ic)
    else
        return 0.
    end
end

"""
    element_lc(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of core orbital anugular momentum (Robicheaux2018 Eq. 24)
"""
function element_lc(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.sr == q2.sr && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(q1.lc) * G4(q1, q2) * G6(q1, q2, ic) * G7(q1, q2)
    else
        return 0.
    end
end

"""
    element_sc(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Returns the reduced matrix element of core spin (Robicheaux2018 Eq. 24)
"""
function element_sc(ic, q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    if q1.sc == q2.sc && q1.sr == q2.sr && q1.lc == q2.lc && q1.lr == q2.lr
        return Λ(0.5) * G4(q1, q2) * G6(q1, q2, ic) * G8(q1, q2)
    else
        return 0.
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
    return sqrt(prod(2a.+1))
end

"""
    G1(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G1(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.
    if q1.Jc == q2.Jc && q1.Fc == q2.Fc
        Λ = q1.Fc + q2.Jr + q1.F + 1
        sq = square_brakets([q1.F, q2.F])
        qn = Vector{Int64}(2*[q1.Jr, q1.F, q1.Fc, q2.F, q2.Jr, 1])
        g = (-1.)^Λ * sq * sf_coupling_6j(qn...)
    end
    return g
end

"""
    G2(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G2(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.
    if q1.sr == q2.sr && q1.lr == q2.lr
        Λ = q1.sr + q1.lr + q1.Jr + 1
        sq = square_brakets([q1.Jr, q2.Jr])
        qn = Vector{Int64}(2*[q1.lr, q1.Jr, q1.sr, q2.Jr, q1.lr, 1])
        g = (-1.)^Λ * sq * sf_coupling_6j(qn...)
    end
    return g
end

"""
    G3(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G3(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.
    if q1.sr == q2.sr && q1.lr == q2.lr
        Λ = q1.sr + q1.lr + q2.Jr + 1
        sq = square_brakets([q1.Jr, q2.Jr])
        qn = Vector{Int64}(2*[q1.sr, q1.Jr, q1.lr, q2.Jr, q1.sr, 1])
        g = (-1.)^Λ * sq * sf_coupling_6j(qn...)
    end
    return g
end

"""
    G4(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G4(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.
    if q1.Jr == q2.Jr
        Λ = q1.Fc + q1.Jr + q2.F + 1
        sq = square_brakets([q1.F, q2.F])
        qn = Vector{Int64}(2*[q1.Fc, q1.F, q1.Jr, q2.F, q2.Fc, 1])
        g = (-1.)^Λ * sq * sf_coupling_6j(qn...)
    end
    return g
end

"""
    G5(q1::fjQuantumNumbers, q2::fjQuantumNumbers, ic)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G5(q1::fjQuantumNumbers, q2::fjQuantumNumbers, ic)
    g = 0.
    if q1.Jc == q2.Jc
        Λ = q1.Jc + ic + q1.Fc + 1
        sq = square_brakets([q1.Fc, q2.Fc])
        qn = Vector{Int64}(2*[ic, q1.Fc, q1.Jc, q2.Fc, ic, 1])
        g = (-1.)^Λ * sq * sf_coupling_6j(qn...)
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
    g = (-1.)^Λ * sq * sf_coupling_6j(qn...)
    return g
end

"""
    G7(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G7(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.
    if q1.sc == q2.sc && q1.lc == q2.lc
        Λ = q1.sc + q1.lc + q1.Jc + 1
        sq = square_brakets([q1.Jc, q2.Jc])
        qn = Vector{Int64}(2*[q1.lc, q1.Jc, q1.sc, q2.Jc, q1.lc, 1])
        g = (-1.)^Λ * sq * sf_coupling_6j(qn...)
    end
    return g
end

"""
    G8(q1::fjQuantumNumbers, q2::fjQuantumNumbers)

Auxiliary function to calculate reduced matrix elements for the magnetic moment (Robicheaux2018 Eq. 25)
"""
function G8(q1::fjQuantumNumbers, q2::fjQuantumNumbers)
    g = 0.
    if q1.sc == q2.sc && q1.lc == q2.lc
        Λ = q1.sc + q1.lc + q2.Jc + 1
        sq = square_brakets([q1.Jc, q2.Jc])
        qn = Vector{Int64}(2*[q1.sc, q1.Jc, q1.lc, q2.Jc, q1.sc, 1])
        g = (-1.)^Λ * sq * sf_coupling_6j(qn...)
    end
    return g
end

# --------------------------------------------------------
# matrix elements for multipole and magnetic interactions
# --------------------------------------------------------

"""
See also [`multipole_moment`](@ref), [`magnetic_dipole_moment`](@ref), [`special_quadrupol_moment`](@ref)

    matrix_element(k::Int, B::BasisArray)
    matrix_element(B::BasisArray)
    matrix_element(A::Parameters, B::BasisArray)

Iterates `multipole_moment`, `magnetic_dipole_moment`, and `special_quadrupol_moment` over a list of states passed as `BasisArray`
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
        for j in i:length(st)
            b2 = st[j]
            f2 = b2.f
            m = multipole_moment(k, b1, b2)
            if !iszero(m)
                M[i,j] = m
                """
                # no need to save the lower triangle
                if j != i
                    M[j,i] = (-1)^(f2-f1) * m
                end
                """
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
        for j in i:length(st)
            b2 = st[j]
            f2 = b2.f
            m = special_quadrupole_moment(b1, b2)
            if !iszero(m)
                M[i,j] = m
                """
                # no need to save the lower triangle
                if j != i
                    M[j,i] = (-1)^(f2-f1) * m
                end
                """
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
        for j in i:length(st)
            b2 = st[j]
            f2 = b2.f
            m = magnetic_dipole_moment(nd, mp, ic, b1, b2)
            if !iszero(m)
                M[i,j] = m
                """
                # no need to save the lower triangle
                if j != i
                    M[j,i] = (-1)^(f2-f1) * m
                end
                """
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

struct QuantumNumbers
    sc::Vector{Number}
    lc::Vector{Number}
    sr::Vector{Number}
    lr::Vector{Number}
    ic::Number
end

function ls_channels(Q::QuantumNumbers)
    sc_list = Q.sc
    lc_list = Q.lc
    sr_list = Q.sr
    lr_list = Q.lr
    Sc = Float64[]
    Sr = Float64[]
    St = Float64[]
    Lc = Int[]
    Lr = Int[]
    Lt = Int[]
    Jt = Float64[]
    for sc in sc_list
        for sr in sr_list
            for st in abs(sc-sr):sc+sr
                for lc in lc_list
                    for lr in lr_list
                        for lt in abs(lc-lr):lc+lr
                            for jt in abs(lt-st):lt+st
                                push!(Sc, sc)
                                push!(Sr, sr)
                                push!(St, st)
                                push!(Lc, lc)
                                push!(Lr, lr)
                                push!(Lt, lt)
                                push!(Jt, jt)
                            end
                        end
                    end
                end
            end
        end
    end
    return lsChannels(Sc, Sr, St, Lc, Lr, Lt, Jt)
end

function jj_channels(Q::QuantumNumbers)
    sc_list = Q.sc
    lc_list = Q.lc
    sr_list = Q.sr
    lr_list = Q.lr
    Sc = Float64[]
    Lc = Int[]
    Jc = Float64[]
    Sr = Float64[]
    Lr = Int[]
    Jr = Float64[]
    Jt = Float64[]
    for sc in sc_list
        for lc in lc_list
            for jc in abs(sc-lc):sc+lc
                for sr in sr_list
                    for lr in lr_list
                        for jr in abs(sr-lr):sr+lr
                            for jt in abs(jc-jr):jc+jr
                                push!(Sc, sc)
                                push!(Lc, lc)
                                push!(Jc, jc)
                                push!(Sr, sr)
                                push!(Lr, lr)
                                push!(Jr, jr)
                                push!(Jt, jt)
                            end
                        end
                    end
                end
            end
        end
    end
    return jjChannels(Sc, Lc, Jc, Sr, Lr, Jr, Jt)
end

function ff_channels(Q::QuantumNumbers)
    sc_list = Q.sc
    lc_list = Q.lc
    sr_list = Q.sr
    lr_list = Q.lr
    ic = Q.ic
    Sc = Float64[]
    Lc = Int[]
    Jc = Float64[]
    Fc = Float64[]
    Sr = Float64[]
    Lr = Int[]
    Jr = Float64[]
    Ft = Float64[]
    for sc in sc_list
        for lc in lc_list
            for jc in abs(sc-lc):sc+lc
                for fc in abs(ic-jc):ic+jc
                    for sr in sr_list
                        for lr in lr_list
                            for jr in abs(sr-lr):sr+lr
                                for ft in abs(fc-jr):fc+jr
                                    push!(Sc, sc)
                                    push!(Lc, lc)
                                    push!(Jc, jc)
                                    push!(Fc, fc)
                                    push!(Sr, sr)
                                    push!(Lr, lr)
                                    push!(Jr, jr)
                                    push!(Ft, ft)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return fjChannels(Sc, Lc, Jc, Fc, Sr, Lr, Jr, Ft)
end

function ff_channels(C::fjChannels)
    sc = C.sc
    lc = C.lc
    Jc = C.Jc
    Fc = C.Fc
    sr = C.sr
    lr = C.lr
    Jr = C.Jr
    Ft = C.F
    values = sort(unique(Ft))
    C = Vector{Channels}()
    for f in values
        c = Ft .== f
        push!(C, fjChannels(sc[c], lc[c], Jc[c], Fc[c], sr[c], lr[c], Jr[c], Ft[c]))
    end
    return values, C
end

function ls_to_jj(sc1, sr1, S, lc1, lr1, L, J1, sc2, lc2, Jc, sr2, lr2, Jr, J2)
    res = 0.
    if (sc1, lc1, sr1, lr1, J1) == (sc2, lc2, sr2, lr2, J2)
        sq = square_brakets([S, L, Jc, Jr])
        qn = Vector{Int64}(2*[sc1, sr1, S, lc1, lr1, L, Jc, Jr, J1])
        res = sq*f9j(qn...)
    end
    return res 
end

function jj_to_fj(sc1, lc1, Jc1, sr1, lr1, Jr1, J, sc2, lc2, Jc2, Fc, sr2, lr2, Jr2, Ft, ic)
    res = 0.
    if (sc1, lc1, sr1, lr1, Jc1, Jr1) == (sc2, lc2, sr2, lr2, Jc2, Jr2)
        Λ_io = Jr1 + Fc - ic - J
        sq = square_brakets([J, Fc])
        qn = Vector{Int64}(2*[Jr1, Jc1, J, ic, Ft, Fc])
        res = (-1.)^Λ_io*sq*sf_coupling_6j(qn...)
    end
    return res 
end

function matrix_ls_to_jj(in::lsChannels, out::jjChannels)
    t = Matrix{Float64}(undef, size(in), size(out))
    for i in axes(t, 1)
        for j in axes(t, 2)
            t[i,j] = ls_to_jj(cat(in)[:,i]..., cat(out)[:,j]...)
        end
    end
    return t
end

function matrix_jj_to_fj(in::jjChannels, out::fjChannels, ic::Number)
    t = Matrix{Float64}(undef, size(in), size(out))
    for i in axes(t, 1)
        for j in axes(t, 2)
            t[i,j] = jj_to_fj(cat(in)[:,i]..., cat(out)[:,j]..., ic)
        end
    end
    return t
end
