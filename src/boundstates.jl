# --------------------------------------------------------
# energy dependence
# --------------------------------------------------------

"""
See also [`theta_rr`](@ref)

    theta(N::Number, M::Matrix)
    theta(N::Vector, M::Matrix)

Given `M` contains quantum defect parameters, where the k-th column represents the (k-1)st order of energy dependence,
this function returns the specific quantum defect at principal quantum number `N`, where `N` can also be a list.

# Examples

```jldoctest
MQDT.theta(60, [0.4 10 100; 0.4 1 0])

# output

2-element Vector{Float64}:
 0.4027854938271605
 0.4002777777777778
```
"""
function theta(N::Number, M::Matrix)
    l, w = size(M)
    t = zeros(l)
    for i in 1:w
        t .+= M[:, i] * (1/N^2)^(i-1)
    end
    return t
end

function theta(N::Vector, M::Matrix)
    l, w = size(M)
    t = zeros(l)
    for i in 1:w
        t .+= M[:, i] .* (1 ./ N .^ 2) .^ (i-1)
    end
    return t
end

"""
See also [`theta`](@ref)

    theta_rr(N::Number, M::Matrix, index::Int)

Similar to `theta`, but here the energy dependence is assumed to be according to the Rydberg-Ritz formula.
Given `M` contains quantum defect parameters, where the k-th column represents the (k-1)st order of energy dependence,
this function returns the specific quantum defect at principal quantum number `N`, where `N` can also be a list.

# Examples

```jldoctest
MQDT.theta_rr(60, [0.4 10 100; 0.4 1 0], 1)

# output

0.4028231137913411
```
"""
function theta_rr(N::Number, M::Matrix, index::Int)
    if !isempty(index)
        m = M[index, :]
        t = 0.0
        for i in eachindex(m)
            t += m[i] * (1/(N-m[1])^2)^(i-1)
        end
        return t
    end
end

"""
    beta(N, L) = π .* (N .- L)
"""
function beta(N, L)
    return π .* (N .- L)
end

"""
    nu(N::Number, T::Vector, I::Number, R::Number)
    nu(N::Number, M::Model, P::Parameters)
    nu(N::Vector, M::Model, P::Parameters)

Given a principal quantum number, return the channel-dependent principal quantum numbers, depending on the channel thresholds and the Rydberg constant.

# Examples

```jldoctest
MQDT.nu(60, [50442.795744, 83967.7, 80835.39, 77504.98, 50443.217463], 50443.217463, 109736.9635066)

# output

5-element Vector{Float64}:
 60.41940064523729
  1.80841488964728
  1.8992315361002043
  2.0125837973909384
 60.0
```
"""
function nu(N::Number, T::Vector, I::Number, R::Number)
    return 1 ./ sqrt.((T .- I) ./ R .+ 1/N^2)
end

function nu(N::Number, M::Model, P::Parameters)
    thresholds = get_thresholds(M, P)
    i = P.threshold
    r = P.rydberg
    return nu(N, thresholds, i, r)
end

function nu(N::Vector, M::Model, P::Parameters)
    n = Matrix{Float64}(undef, M.size, length(N))
    for i in eachindex(N)
        n[:, i] = nu(N[i], M, P)
    end
    return n
end

"""
    epsilon(N, T, R)
    epsilon(N, P::Parameters)

Given a principal quantum number, return the energy, depending on the ionization thresholds and the Rydberg constant.

# Examples

```jldoctest
ϵ = MQDT.epsilon(60, 50443.217463, 109736.9635066)

# output

50412.734973137056
```
"""
function epsilon(N, T, R)
    return T .- R ./ N .^ 2
end

function epsilon(N, P::Parameters)
    t = P.threshold
    r = P.rydberg
    return epsilon(N, t, r)
end

# --------------------------------------------------------
# k matrix
# --------------------------------------------------------

"""
    kbar(M::Number) = tan(π * M)
    kbar(M::Vector{Float64}) = diagm(kbar.(M))
    kbar(N::Number, M::fModel, P::Parameters)

Function that returns the K matrix in the close-coupling frame, i.e. in its diagonal representation
"""
function kbar(M::Number)
    return tan(π * M)
end

function kbar(M::Vector{Float64})
    return diagm(kbar.(M))
end

function kbar(N::Number, M::fModel, P::Parameters)
    n = nu(N, M, P)
    m = theta(n, M.defects)
    return kbar(m)
end

# function kbar(N::Number, M::fModel, P::Parameters)
#     n = nu(N, M, P)
#     m = theta(n, M.defects)
#     i = findall(isone, M.rydbergritz)
#     if !isempty(i)
#         for j in i
#             m[j] = theta_rr(n[j], M.defects, j)
#         end
#     end
#     return kbar(m)
# end

"""
    couple(C::String)

From a string of type "ij", return a touple of indices (i, j).
With this method, maximum for j is 10.

# Examples

```jldoctest
MQDT.couple("710")

# output

(7, 10)
```
"""
function couple(C::String)
    if isempty(C)
        C = "11"
    end
    c = parse(Int, C)
    i = div(c, 10)
    j = mod(c, 10)
    if j == 0
        i = div(c, 100)
        j = mod(c, 100)
    end
    return i, j
end

"""
    rot(T::Number, C::String, D::Int)
    rot(T::Vector, C::Vector{String}, D::Int)
    rot(N::Number, M::fModel, P::Parameters)

Function that returns a rotation matrix of dimensions `D` with rotation angle `T`.
Works for subsequent rotations passed as vectors.

# Examples

```jldoctest
println(MQDT.rot(0.4, "23", 3))
println(MQDT.rot([0.4, 0.2], ["23", "12"], 3))

# output

[1.0 0.0 0.0; 0.0 0.9210609940028851 -0.3894183423086505; 0.0 0.3894183423086505 0.9210609940028851]
[0.9800665778412416 -0.19866933079506122 0.0; 0.18298657129998708 0.90270109637546 -0.3894183423086505; 0.07736548146578168 0.3816559020950483 0.9210609940028851]
```
"""
function rot(T::Number, C::String, D::Int)
    r = [cos(T) -sin(T); sin(T) cos(T)]
    out = diagm(ones(D))
    i, j = couple(C)
    out[[i, j], [i, j]] = r
    return out
end

function rot(T::Vector, C::Vector{String}, D::Int)
    r = diagm(ones(D))
    for i in eachindex(T)
        r = r * rot(T[i], C[i], D)
    end
    return r
end

function rot(N::Number, M::fModel, P::Parameters)
    if iszero(M.angles) # check presence of channel mixing
        return diagm(ones(M.size))
    else
        t = M.angles[:, 1]
        i = findall(!iszero, M.angles[:, 2]) # check energy dependence
        if !isempty(i)
            j = couple(M.mixing[i[1]])[1] # find involved channel
            n = nu(N, M, P)[j] # find nu relative to threshold
            t = theta(n, M.angles)
        end
        return rot(t, M.mixing, M.size)
    end
end

"""
    transform(N::Number, M::fModel, P::Parameters)

Function that returns the frame transformation matrix from close-coupling to fragmentation frame.
"""
function transform(N::Number, M::fModel, P::Parameters)
    r = rot(N, M, P)
    u = M.unitary
    return u * r
end

"""
    kmat(N::Number, M::fModel, P::Parameters)
    kmat(N::Number, M::kModel, P::Parameters)

Function that returns the K matrix in the fragmentation frame.
"""
function kmat(N::Number, M::fModel, P::Parameters)
    k = kbar(N, M, P)
    t = transform(N, M, P)
    return t * k * t'
end

function kmat(N::Number, M::kModel, P::Parameters)
    k = M.K0 + diagm(M.K1 * (1 - epsilon(N, P)/P.threshold))
    return k
end

"""
    mmat(N::Number, M::Model, P::Parameters)

Function that returns the M matrix.
"""
function mmat(N::Number, M::Model, P::Parameters)
    k = kmat(N, M, P)
    n = nu(N, M, P)
    l = get_lr(M)
    b = beta(n, l)
    return diagm(sin.(b)) .+ cos.(b)' .* k
end

# --------------------------------------------------------
# root and nullspace
# --------------------------------------------------------

"""
    mroots(N::Number, M::Model, P::Parameters)
    mroots(N1::Number, N2::Number, M::Model, P::Parameters)

Function that returns the roots of the M matrix in the range from `N` to `N+1` or, if provided, from `N1` to `N2`.
Root finding algorithm used from the package `Roots`.
"""
function mroots(N::Number, M::Model, P::Parameters)
    m(n) = det(mmat(n, M, P))
    z = find_zeros(m, (N-0.5, N+0.5))
    c = Int[]
    for i in eachindex(z)
        mz = m(z[i])
        if abs(mz) > 1e-10
            println("Warning: skipped a root that seems inaccurate. Value was $(mz) at n=$(z[i]) for $(M.name).")
        else
            push!(c, i)
        end
    end
    return z[c]
end

function mroots(N1::Number, N2::Number, M::Model, P::Parameters)
    z = Float64[]
    for i in N1:N2
        append!(z, mroots(i, M, P))
    end
    return z
end

"""
See also [`mroots`](@ref)

    eigenenergies(N1::Number, N2::Number, M::Model, P::Parameters)

Function that returns both the roots of the M matrix as well as the channel-dependent principal quantum numbers.
"""
function eigenenergies(N1::Number, N2::Number, M::Model, P::Parameters)
    z = mroots(N1, N2, M, P)
    n = nu(z, M, P)
    return z, n
end

"""
    nullspace(n::Vector{Float64}, M::Matrix)

Returns the nullspace of the M matrix, which corresponds to the channel coefficents for an MQDT bound state.
"""
function LinearAlgebra.nullspace(n::Vector{Float64}, M::Matrix)
    m = M ./ n' .^ (3/2)
    e = eigen(m)
    s = sortperm(abs.(e.values))
    i = s[1]
    f = e.values[i]
    if abs(f) > 1e-11
        println("Warning: nullspace may not be accurate. Smallest eigenvalue is $f for nu=$(n[1]).")
    end
    if length(s) > 1
        g = e.values[s[2]]
        if abs(g) < 1e-14
            println("Warning: nullspace may not be unique. Second smallest eigenvalue is $g for nu=$(n[1]).")
        end
    end
    return e.vectors[:, i]
end

"""
See also [`eigenenergies`](@ref)

    eigenstates(N1::Number, N2::Number, M::Model, P::Parameters)

Function that returns the bound states corresponding to an MQDT model in the form of an instance of `EigenStates`,
which contains the global reference principal quantum number, the channel-dependent principal quantum numbers, as well as the channel coefficients.
"""
function eigenstates(N1::Number, N2::Number, M::Model, P::Parameters)
    z, n = eigenenergies(N1, N2, M, P)
    mfunc(e) = mmat(e, M, P)
    a = similar(n)
    for i in eachindex(z)
        m = mfunc(z[i])
        t = nullspace(n[:, i], m)
        a[:, i] = t
    end
    return EigenStates(z, n, a)
end

# --------------------------------------------------------
# basis array
# --------------------------------------------------------

function exp_LS(A::Matrix{Float64}, T::Matrix{Float64}, V)
    if allequal(V)
        return repeat([V[1]], size(A, 2))
    else
        t = T' * A
        return sum(t .^ 2 .* V; dims=1)[:]
    end
end

function find_leading_term(A::Matrix{Float64}, T::Matrix{Float64}, L::Vector{String})
    term = Vector{String}(undef, size(A, 2))
    lead = Vector{Float64}(undef, size(A, 2))
    for i in axes(A, 2)
        val, ind = findmax((T' * A[:, i]) .^ 2)
        term[i] = L[ind]
        lead[i] = val
    end
    return term, lead
end

"""
    basisarray(T::EigenStates, M::fModel)
    basisarray(T::Vector{EigenStates}, M::Vector{fModel})

Function that generates all relevant bound-state data.
Returns a `BasisArray`, which is a list of `BasisState` instances.
"""
function basisarray(T::EigenStates, M::fModel)
    F = findall(M.core)
    B = Vector{BasisState}()
    c = M.outer_channels
    e = T.n
    p = unique_parity(c)
    f = good_quantum_number(c)
    n = T.nu[F, :]
    l = get_lr(c)
    a = T.a[F, :]
    term, lead = find_leading_term(T.a, M.unitary, M.terms)
    for i in eachindex(e)
        if !isone(M.size) || l[1] < n[1, i]
            ei = e[i]
            if abs(ei - round(Int, ei)) < 1e2eps()
                ei = round(ei)
            end
            push!(B, BasisState(M.species, ei, p, f, n[:, i], l, a[:, i], c, term[i], lead[i]))
        end
    end
    return BasisArray(B)
end

function basisarray(T::Vector{EigenStates}, M::Vector{fModel})
    B = Vector{BasisState}()
    for i in eachindex(T)
        append!(B, basisarray(T[i], M[i]).states)
    end
    ordered = sortperm(get_nu(BasisArray(B)))
    return BasisArray(B[ordered])
end
