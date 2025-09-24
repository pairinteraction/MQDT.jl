# --------------------------------------------------------
# Store Parameters
# --------------------------------------------------------

"""
    Parameters(species::Symbol, mass::Float64, spin::Float64, rydberg::Float64, threshold::Float64, hyperfine::Float64, dipole::Float64)

Type to store relevant parameters for each atomic element. A corresponding
file in the parameters folder of the mqdt module initializes this type.

# Examples

```jldoctest
julia> PARA = Parameters(
           :Yb174,
           1822.88848628*173.9388621, # nuclear mass
           0, # nuclear spin
           109736.9695858, # Rydberg constant in 1/cm
           50443.070393, # lowest ionization threshold in 1/cm
           0, # hyperfine constant in 1/cm
           2.1, # nuclear dipole
       );

```
"""
struct Parameters
    species::Symbol
    mass::Float64
    spin::Float64
    rydberg::Float64
    threshold::Float64
    hyperfine::Float64
    dipole::Float64
end

# --------------------------------------------------------
# Quantum Number structs
# --------------------------------------------------------

"""
See also [`lsQuantumNumbers`](@ref), [`jjQuantumNumbers`](@ref), [`fjQuantumNumbers`](@ref)

Abstract type to store the quantum numbers of an MQDT channel representation.
"""
abstract type QuantumNumbers end

"""
See also [`lsQuantumNumbers`](@ref), [`jjQuantumNumbers`](@ref)

    fjQuantumNumbers(sc::Float64, lc::Int, Jc::Float64, Fc::Float64, lr::Int, Jr::Float64, F::Float64)

Type to store the quantum numbers of a single MQDT channel in the "fj" (sometimes called "fragmentation")
coupling scheme, i.e. core angular momentum quantum numbers are coupled to the nuclear spin to give a
core hyperfine configuration (Fc) which is coupled to the total Rydberg angular momentum (Jr).
"""
struct fjQuantumNumbers <: QuantumNumbers
    sc::Float64
    lc::Int
    Jc::Float64
    Fc::Float64
    lr::Int
    Jr::Float64
    F::Float64
end

"""
See also [`lsQuantumNumbers`](@ref), [`fjQuantumNumbers`](@ref)

    jjQuantumNumbers(lc::Int, Jc::Float64, lr::Int, Jr::Float64, J::Float64)

Type to store the quantum numbers of a single MQDT channel in the "jj" coupling scheme,
i.e. core and Rydberg angular momentum (AM) quantum numbers are coupled separately to give a total core AM (Jc) and total Rydberg AM (Jr).
"""
struct jjQuantumNumbers <: QuantumNumbers
    sc::Float64
    lc::Int
    Jc::Float64
    lr::Int
    Jr::Float64
    J::Float64
    F::Float64
end

function jjQuantumNumbers(sc, lc, Jc, lr, Jr, J)
    return jjQuantumNumbers(sc, lc, Jc, lr, Jr, J, J)
end

"""
See also [`jjQuantumNumbers`](@ref), [`fjQuantumNumbers`](@ref)

    lsQuantumNumbers(S::Float64, lc::Int, lr::Int, L::Int, J::Float64)

Type to store the quantum numbers of a single MQDT channel in the "LS" coupling scheme,
i.e. core and Rydberg angular momentum quantum numbers are coupled individually to give a total spin (S) and total orbital angular momentum (L).
"""
struct lsQuantumNumbers <: QuantumNumbers
    sc::Float64
    S::Float64
    lc::Int
    lr::Int
    L::Int
    J::Float64
    F::Float64
end

function lsQuantumNumbers(sc, S, lc, lr, L, J)
    return lsQuantumNumbers(sc, S, lc, lr, L, J, J)
end

# --------------------------------------------------------
# Base functions for QuantumNumbers structs
# --------------------------------------------------------

function Base.cat(T::lsQuantumNumbers)
    return Vector{Float64}([T.sc, T.S, T.lc, T.lr, T.L, T.J, T.F])
end

function Base.cat(T::jjQuantumNumbers)
    return Vector{Float64}([T.sc, T.lc, T.Jc, T.lr, T.Jr, T.J, T.F])
end

function Base.cat(T::fjQuantumNumbers)
    return Vector{Float64}([T.sc, T.lc, T.Jc, T.Fc, T.lr, T.Jr, T.F])
end

function fj_quantum_numbers(T::jjQuantumNumbers)
    return fjQuantumNumbers(T.sc, T.lc, T.Jc, T.Jc, T.lr, T.Jr, T.J)
end

# --------------------------------------------------------
# Channel structs
# --------------------------------------------------------

"""
See also [`lsChannels`](@ref), [`jjChannels`](@ref), [`fjChannels`](@ref)

Abstract type to store the quantum numbers of an MQDT channel representation.
"""
abstract type Channels end

"""
See also [`jjChannels`](@ref), [`fjChannels`](@ref)

    lsChannels(i::Vector{lsQuantumNumbers})

Type to store the quantum numbers of an MQDT channel representation, here in the "LS" coupling scheme,
i.e. core and Rydberg angular momentum quantum numbers are coupled individually to give a total spin (S) and total orbital angular momentum (L).
"""
struct lsChannels <: Channels
    i::Vector{lsQuantumNumbers}
end

"""
See also [`lsChannels`](@ref), [`fjChannels`](@ref)

    jjChannels(i::Vector{jjQuantumNumbers})

Type to store the quantum numbers of an MQDT channel representation, here in the
"jj" coupling scheme, i.e. core and Rydberg angular momentum (AM) quantum numbers
are coupled separately to give a total core AM (Jc) and total Rydberg AM (Jr).
"""
struct jjChannels <: Channels
    i::Vector{jjQuantumNumbers}
end

"""
See also [`lsChannels`](@ref), [`jjChannels`](@ref)

    fjChannels(i::Vector{fjQuantumNumbers})

Type to store the quantum numbers of an MQDT channel representation, here in the "fj" (sometimes called
"fragmentation") coupling scheme
"""
struct fjChannels <: Channels
    i::Vector{fjQuantumNumbers}
end

# --------------------------------------------------------
# Base functions for Channel structs
# --------------------------------------------------------

function Base.size(T::Channels)
    return length(T.i)
end

function get_S(T::lsChannels)
    t = T.i
    r = Vector{Int}(undef, length(t))
    for i in eachindex(t)
        r[i] = t[i].S
    end
    return r
end

function get_S(T::jjChannels)
    t = T.i
    r = Vector{Float64}(undef, length(t))
    for i in eachindex(t)
        lr = t[i].lr
        jr = t[i].Jr
        jt = t[i].J
        theta = atan(sqrt(lr/(lr+1)))
        if jt == lr
            if jr > lr
                r[i] = cos(theta)^2
            else
                r[i] = sin(theta)^2
            end
        else
            r[i] = 1.0
        end
    end
    return r
end

function get_lc(T::Channels)
    t = T.i
    r = Vector{Int}(undef, length(t))
    for i in eachindex(t)
        r[i] = t[i].lc
    end
    return r
end

function get_lr(T::Channels)
    t = T.i
    r = Vector{Int}(undef, length(t))
    for i in eachindex(t)
        r[i] = t[i].lr
    end
    return r
end

function get_L(T::lsChannels)
    t = T.i
    r = Vector{Int}(undef, length(t))
    for i in eachindex(t)
        r[i] = t[i].L
    end
    return r
end

function get_L(T::jjChannels)
    return get_lr(T)
end

function get_Jr(T::Channels)
    t = T.i
    r = Vector{Float64}(undef, length(t))
    for i in eachindex(t)
        r[i] = t[i].Jr
    end
    return r
end

function get_J(T::Channels)
    t = T.i
    r = Vector{Float64}(undef, length(t))
    for i in eachindex(t)
        r[i] = t[i].J
    end
    return r
end

function get_F(T::Channels)
    t = T.i
    r = Vector{Float64}(undef, length(t))
    for i in eachindex(t)
        r[i] = t[i].F
    end
    return r
end

"""
See also [`unique_parity`](@ref)

    parity(T::Channels)

Given a channel representation, this function returns the parity of a each
channel based of the orbital angular momenta of the core and the Rydberg electron.
"""
function parity(T::Channels)
    pc = 2iseven.(get_lc(T)) .- 1
    pr = 2iseven.(get_lr(T)) .- 1
    return pc .* pr
end

"""
See also [`parity`](@ref)

    unique_parity(T::Channels)

Given a channel representation, this function checks the parity of
each channel and returns it in case each channel has the same parity.
"""
function unique_parity(T::Channels)
    f = unique(parity(T))
    if length(f) == 1
        return f[1]
    else
        error("unique_parity: channels have different parity")
    end
end

"""
    good_quantum_number(T::lsChannels)
    good_quantum_number(T::jjChannels)
    good_quantum_number(T::fjChannels)

Given a ls or jj (fj) channel representation, this function checks whether each channel
has the same total (hyperfine) angular momentum J (F) and, if true, returns it.
"""
function good_quantum_number(T::lsChannels)
    f = unique(get_J(T))
    if length(f) == 1
        return f[1]
    else
        error("good_quantum_number: channels have different J")
    end
end

function good_quantum_number(T::jjChannels)
    f = unique(get_J(T))
    if length(f) == 1
        return f[1]
    else
        error("good_quantum_number: channels have different J")
    end
end

function good_quantum_number(T::fjChannels)
    f = unique(get_F(T))
    if length(f) == 1
        return f[1]
    else
        error("good_quantum_number: channels have different F")
    end
end

# --------------------------------------------------------
# K matrix model structs
# --------------------------------------------------------

"""
See also [`fModel`](@ref), [`kModel`](@ref)

Abstract type to store MQDT models.
"""
abstract type Model end

"""
See also [`kModel`](@ref)

    fModel(
        species::Symbol,
        name::String,
        size::Int,
        terms::Vector{String},
        core::Vector{Bool},
        thresholds::Vector{Float64},
        defects::Matrix{Float64},
        mixing::Vector{String},
        angles::Matrix{Float64},
        inner_channels::Channels
        outer_channels::Channels
        unitary::Matrix{Float64}
    )

Type to store MQDT models inspired by [PRX 15, 011009 (2025)] using sparse K matriced and frame transformation.
Contains all relevant parameters to form the K matrix and calculate the spectrum for a specific Rydberg series.

# Examples

```jldoctest
julia> RYDBERG_S0 = fModel(
           :Yb174,
           "S J=0, ν > 2", # fit for states 6s7s upward [Phys. Rev. X 15, 011009 (2025)]
           6,
           ["6sns 1S0", "4f13 5d 6snl a", "6pnp 1S0", "4f13 5d 6snl b", "6pnp 3P0", "4f13 5d 6snl c"],
           Bool[1, 0, 1, 0, 1, 0],
           [50443.070393, 83967.7, 80835.39, 83967.7, 77504.98, 83967.7],
           [
               0.355097325 0.278368431;
               0.204537279 0;
               0.116394359 0;
               0.295432196 0;
               0.25765161 0;
               0.155807042 0
           ],
           ["12", "13", "14", "34", "35", "16"],
           [0.12654859 0; 0.30010744 0; 0.05703381 0; 0.11439805 0; 0.09864375 0; 0.14248210 0],
           lsChannels([
               lsQuantumNumbers(0.5, 0, 0, 0, 0, 0),
               lsQuantumNumbers(0.5, 0, 1, 1, 0, 0),
               lsQuantumNumbers(0.5, 1, 1, 1, 1, 0),
           ]),
           jjChannels([
               jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 0),
               jjQuantumNumbers(0.5, 1, 0.5, 1, 1.5, 0),
               jjQuantumNumbers(0.5, 1, 0.5, 1, 0.5, 0),
           ]),
           [
               1 0 0 0 0 0;
               0 1 0 0 0 0;
               0 0 -sqrt(2/3) 0 sqrt(1/3) 0;
               0 0 0 1 0 0;
               0 0 sqrt(1/3) 0 sqrt(2/3) 0;
               0 0 0 0 0 1
           ],
       );

```
"""
struct fModel <: Model
    species::Symbol
    name::String
    size::Int
    terms::Vector{String}
    core::Vector{Bool}
    thresholds::Vector{Float64}
    defects::Matrix{Float64}
    mixing::Vector{String}
    angles::Matrix{Float64}
    inner_channels::Channels
    outer_channels::Channels
    unitary::Matrix{Float64}
end

"""
See also [`fModel`](@ref)

    kModel(
        species::Symbol
        name::String
        size::Int
        terms::Vector{String}
        jjscheme::Vector{Bool}
        lschannels::lsChannels
        jjchannels::jjChannels
        thresholds::Vector{Float64}
        K0::Matrix{Float64}
        K1::Vector{Float64}
    )

Type to store MQDT models inspired by [JPB 47, 155001 (2014)] using dense, energy-dependent K matrices.
Contains all relevant parameters to form the K matrix and calculate the spectrum for a specific Rydberg series.

# Examples

```jldoctest
julia> KMODEL_S0 = kModel(
           :Sr88,
           "1S0",
           3,
           ["(5s1/2)(ns1/2)", "(4d5/2)(nd5/2)", "(4d3/2)(nd3/2)"],
           Bool[1, 1, 1],
           lsChannels([
               lsQuantumNumbers(0.5, 0, 0, 0, 0, 0),
               lsQuantumNumbers(0.5, 0, 2, 2, 0, 0),
               lsQuantumNumbers(0.5, 0, 2, 2, 0, 0),
           ]),
           jjChannels([
               jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 0),
               jjQuantumNumbers(0.5, 2, 2.5, 2, 2.5, 0),
               jjQuantumNumbers(0.5, 2, 1.5, 2, 1.5, 0),
           ]),
           [45932.2002, 60768.43, 60488.09],
           [
               1.051261 0.3759864 -0.02365485;
               0.3759864 -0.6400925 -0.0002063825;
               -0.02365485 -0.0002063825 3.009087
           ],
           [0.8763911, 0.4042584, 17.22631],
       );

```
"""
struct kModel <: Model
    species::Symbol
    name::String
    size::Int
    terms::Vector{String}
    jjscheme::Vector{Bool}
    lschannels::lsChannels
    jjchannels::jjChannels
    thresholds::Vector{Float64}
    K0::Matrix{Float64}
    K1::Vector{Float64}
end

# --------------------------------------------------------
# Base functions for model structs
# --------------------------------------------------------

function Base.length(T::Model)
    return T.size
end

function Base.size(T::Model)
    return (T.size,)
end

function get_lr(T::fModel)
    l = zeros(Int, T.size)
    l[findall(T.core)] = get_lr(T.outer_channels)
    return l
end

function get_lr(T::kModel)
    return get_lr(T.jjchannels)
end

function get_S(T::fModel)
    return get_S(T.inner_channels)
end

function get_L(T::fModel)
    return get_L(T.inner_channels)
end

function get_J(T::fModel)
    return get_J(T.inner_channels)
end

function single_channel_models(species::Symbol, l::Integer, p::Parameters)
    @assert l > 0 "l must be positive and nonzero for this function"
    jr = [l-1/2, l-1/2, l+1/2, l+1/2]
    jt = [l-1, l, l, l+1]
    m = Vector{fModel}(undef, 4)
    for i in eachindex(jt)
        m[i] = fModel(
            species,
            "L=$l, J=$(jt[i]), Jr=$(jr[i])",
            1,
            [""],
            Bool[1],
            [p.threshold],
            [0;;],
            [""],
            [0;;],
            jjChannels([jjQuantumNumbers(0.5, 0, 0.5, l, jr[i], jt[i])]),
            jjChannels([jjQuantumNumbers(0.5, 0, 0.5, l, jr[i], jt[i])]),
            [1;;],
        )
    end
    return m
end

function single_channel_models(species::Symbol, l_list::UnitRange{Int64}, p::Parameters)
    m = Vector{fModel}()
    for l in l_list
        append!(m, single_channel_models(species, l, p))
    end
    return m
end

function test_model(T::fModel)
    if T.size != length(T.terms) != length(T.core) != length(T.thresholds) != size(T.defects, 1)
        return println("Model size does not correspond to provided parameters.")
    elseif size(T.unitary) != (T.size, T.size)
        return println("Frame transformation matrix does not have the correct dimensions.")
    elseif length(findall(T.core)) != size(T.outer_channels) != size(T.inner_channels)
        return println("Wrong number of channels.")
    elseif length(T.mixing) != size(T.angles, 1)
        return println("Wrong number of close-coupling rotations.")
    else
        return println("Model $(T.name) passed.")
    end
end

function test_model(T::kModel)
    if T.size !=
       length(T.terms) !=
       length(jjscheme) !=
       length(T.lschannels) !=
       length(T.jjchannels) !=
       length(T.thresholds) !=
       length(T.K1)
        return println("Model size does not correspond to provided parameters.")
    elseif size(T.K0) != (T.size, T.size)
        return println("Frame transformation matrix does not have the correct dimensions.")
    else
        return println("Model $(T.name) passed.")
    end
end

function test_model(T::Union{Vector{fModel},Vector{kModel}})
    for i in eachindex(T)
        test_model(T[i])
    end
end

function test_unitary(T::fModel, P::Parameters)
    c = findall(T.core)
    t0 = T.unitary[c, c]
    if typeof(T.outer_channels) == jjChannels
        if typeof(T.inner_channels) == lsChannels
            t1 = matrix_ls_to_jj(T.inner_channels, T.outer_channels)'
        else
            s = size(T.inner_channels)
            t1 = diagm(repeat([1.0], s))
        end
    elseif typeof(T.outer_channels) == fjChannels
        if typeof(T.inner_channels) == jjChannels
            t1 = matrix_jj_to_fj(T.inner_channels, T.outer_channels, P.spin)'
        else
            lc = unique(get_lc(T.inner_channels))
            lr = unique(get_lr(T.inner_channels))
            ft = unique(get_F(T.outer_channels))[1]
            qn = mqdt.AngularMomenta(lc, lr, P.spin)
            jj = jj_channels(qn)
            fj = [unique.(get_F.(jj))[i][1] for i in 1:length(jj)]
            i = findfirst(isequal(ft), fj)
            tjj = matrix_ls_to_jj(T.inner_channels, jj[i])
            tfj = matrix_jj_to_fj(jj[i], T.outer_channels, P.spin)
            t1 = (tjj * tfj)'
        end
    end
    return isapprox(t0, t1), t0, t1
end

function test_unitary(T::kModel)
    t1 = eigen(T.K0).vectors
    t2 = matrix_ls_to_jj(T.lschannels, T.jjchannels)'
    return t1, t2
end

# --------------------------------------------------------
# Bound state structs
# --------------------------------------------------------

"""
See also [`eigenstates`](@ref), [`BasisState`](@ref)

    EigenStates(n::Vector{Float64}, nu::Matrix{Float64}, a::Matrix{Float64})

Type to store multi-channel bound states.
`n` is a global energy reference to a global threshold specified in the `Parameters` type.
`nu` is the principal quantum number in each channel with respect to the channel-specific threshold.
`a` is the coefficient of each channel contributing to the bound state.
`EigenStates` is populated by the `eigenstates` function.
"""
struct EigenStates
    n::Vector{Float64}
    nu::Matrix{Float64}
    a::Matrix{Float64}
end

"""
See also [`basisarray`](@ref), [`BasisArray`](@ref), [`EigenStates`](@ref)

    BasisState(
        species::Symbol,
        energy::Float64,
        parity::Int,
        f::Float64,
        nu::Vector{Float64},
        lr::Vector{Int},
        coeff::Vector{Float64},
        channels::Channels,
        term::String,
        lead::Float64,
        L::Float64,
        S::Float64,
    )

Type to store all relevant information of multi-channel bound states for a given Rydberg series.
`BasisState` is generated by the `basisarray` function using a `Model` and `EigenStates`.
"""
struct BasisState
    species::Symbol
    energy::Float64
    parity::Int
    f::Float64
    nu::Vector{Float64}
    lr::Vector{Int}
    coeff::Vector{Float64}
    channels::Channels
    term::String
    lead::Float64
    L::Float64
    S::Float64
end

"""
See also [`basisarray`](@ref), [`BasisState`](@ref)

    BasisArray(states::Vector{BasisState})

Type to store a basis. A basis is a list of `BasisState` types, typically corresponding to different Rydberg series.
`BasisArray` is generated by the `basisarray` function using lists of `Model` and `EigenStates`.
"""
struct BasisArray
    states::Vector{BasisState}
end

"""
See also [`databasearray`](@ref), [`DataBaseArray`](@ref), [`EigenStates`](@ref), [`BasisState`](@ref)

    DataBaseState(
        nu::Float64,
        parity::Int,
        f::Float64,
        nui::Vector{Float64},
        nui_all::Vector{Float64},
        ai::Vector{Float64},
        ai_all::Vector{Float64},
        S::Vector{Float64},
        L::Vector{Float64},
        J::Vector{Float64},
        lr::Vector{Int},
        Jr::Vector{Float64},
        transform::Matrix{Float64},
        term::String,
        lead::Float64,
        neg::Float64,
    )

Type to store all relevant to the PAIRINTERACTION database for a given Rydberg series.
`DataBaseState` is generated by the `databasearray` function using a `Model` and `EigenStates`.
"""
struct DataBaseState
    nu::Float64
    parity::Int
    f::Float64
    nui::Vector{Float64}
    nui_all::Vector{Float64}
    ai::Vector{Float64}
    ai_all::Vector{Float64}
    S::Vector{Float64}
    L::Vector{Int}
    J::Vector{Float64}
    lr::Vector{Int}
    Jr::Vector{Float64}
    transform::Matrix{Float64}
    term::String
    lead::Float64
    neg::Float64
end

"""
See also [`DataBaseState`](@ref), [`databasearray`](@ref), [`basisarray`](@ref), [`BasisArray`](@ref)

    DataBaseArray(states::Vector{DataBaseState})

Type to store a bound state data as a list of `DataBaseState` types as used by the PAIRINTERACTION software.
`DataBaseArray` is generated by the `databasearray` function using lists of `Model` and `EigenStates`.
"""
struct DataBaseArray
    states::Vector{DataBaseState}
end

# --------------------------------------------------------
# Base functions for state structs
# --------------------------------------------------------

function Base.length(T::EigenStates)
    return length(T.n)
end

function Base.size(T::EigenStates)
    return size(T.nu)
end

function Base.size(T::BasisArray)
    return length(T.states)
end

function Base.union(T1::BasisArray, T2::BasisArray)
    return BasisArray(vcat(T1.states, T2.states))
end

function Base.union(T1::BasisArray, T2::BasisArray, T3::BasisArray)
    return BasisArray(vcat(T1.states, T2.states, T3.states))
end

function Base.size(T::DataBaseArray)
    return length(T.states)
end

function Base.union(T1::DataBaseArray, T2::DataBaseArray)
    return DataBaseArray(vcat(T1.states, T2.states))
end

function Base.union(T1::DataBaseArray, T2::DataBaseArray, T3::DataBaseArray)
    return DataBaseArray(vcat(T1.states, T2.states, T3.states))
end
