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

"""
    coreQuantumNumbers(lc::Int, Jc::Float64, Fc::Float64)

Type to store the quantum numbers of the core excitations.
If Fc is `NaN`, hyperfine structure is neglected for this core excitation.

# Examples

```jldoctest
julia> core_exc = coreQuantumNumbers(0, 0.5, NaN);

```
"""
struct coreQuantumNumbers <: QuantumNumbers
    lc::Int
    Jc::Float64
    Fc::Float64
end

function coreQuantumNumbers(lc, Jc)
    return coreQuantumNumbers(lc, Jc, NaN)
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
           Dict(coreQuantumNumbers(0, 0.5) => 50443.070393),
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
    thresholds_dict::Dict{Union{coreQuantumNumbers,String},Float64}
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
            if jr < lr
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
        F::Float64,
        range::Tuple{Float64,Float64},
        terms::Vector{String},
        core::Vector{Bool},
        defects::Matrix{Float64},
        mixing::Vector{String},
        angles::Matrix{Float64},
        inner_channels::Channels,
        outer_channels::Channels,
        unitary::Matrix{Float64},
        thresholds_dict::Union{Dict{Union{coreQuantumNumbers,String},Float64},Nothing}
    )

Type to store MQDT models inspired by [PRX 15, 011009 (2025)] using sparse K matriced and frame transformation.
Contains all relevant parameters to form the K matrix and calculate the spectrum for a specific Rydberg series.

# Examples

```jldoctest
julia> RYDBERG_S0 = fModel(
           :Yb174,
           "S J=0, ν > 2", # fit for states 6s7s upward [Phys. Rev. X 15, 011009 (2025)]
           0,
           (2, Inf),
           ["6sns 1S0", "4f13 5d 6snl a", "6pnp 1S0", "4f13 5d 6snl b", "6pnp 3P0", "4f13 5d 6snl c"],
           Bool[1, 0, 1, 0, 1, 0],
           [
               0.355097325 0.278368431;
               0.204537279 0;
               0.116394359 0;
               0.295432196 0;
               0.25765161 0;
               0.155807042 0
           ],
           ["1.2", "1.3", "1.4", "3.4", "3.5", "1.6"],
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
    F::Float64
    range::Tuple{Float64,Float64}
    terms::Vector{String}
    core::Vector{Bool}
    defects::Matrix{Float64}
    mixing::Vector{String}
    angles::Matrix{Float64}
    inner_channels::Channels
    outer_channels::Channels
    unitary::Matrix{Float64}
    thresholds_dict::Union{Dict{Union{coreQuantumNumbers,String},Float64},Nothing}
end

function fModel(species, name, F, range, terms, core, defects, mixing, angles, inner_channels, outer_channels, unitary)
    return fModel(
        species,
        name,
        F,
        range,
        terms,
        core,
        defects,
        mixing,
        angles,
        inner_channels,
        outer_channels,
        unitary,
        nothing,
    )
end

"""
See also [`fModel`](@ref)

    kModel(
        species::Symbol,
        name::String,
        size::Int,
        terms::Vector{String},
        jjscheme::Vector{Bool},
        lschannels::lsChannels,
        jjchannels::jjChannels,
        K0::Matrix{Float64},
        K1::Vector{Float64},
        thresholds_dict::Union{Dict{Union{coreQuantumNumbers,String},Float64},Nothing}
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
    K0::Matrix{Float64}
    K1::Vector{Float64}
    thresholds_dict::Union{Dict{Union{coreQuantumNumbers,String},Float64},Nothing}
end

function kModel(species, name, size, terms, jjscheme, lschannels, jjchannels, K0, K1)
    return kModel(species, name, size, terms, jjscheme, lschannels, jjchannels, K0, K1, nothing)
end

# --------------------------------------------------------
# Base functions for model structs
# --------------------------------------------------------

function Base.length(T::Model)
    return length(T.terms)
end

function get_lr(T::fModel)
    l = zeros(Int, length(T.terms))
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

"""
Get the minimum and maximum effective principal quantum number (nu_min, nu_max) from a given Model.

The name of a Model usually is something like "S J=0, ν > 2" or "S F=1/2, ν > 26",
from this we can extract nu_min = 2 and nu_max Inf.

Some Model names also contain a upper limit, e.g. "P J=1, 1.7 < ν < 2.7",
from which we extract nu_min = 1.7 and nu_max = 2.7.
"""
function get_nu_limits_from_model(model::Model)
    m = match(r"([0-9/\.]+)\s*<\s*ν\s*<\s*([0-9/\.]+)", model.name)
    if m !== nothing
        nu_min = parse(Float64, m.captures[1])
        nu_max = parse(Float64, m.captures[2])
        return nu_min, nu_max
    end

    m = match(r"ν\s*>\s*([0-9/\.]+)", model.name)
    if m !== nothing
        nu_min = parse(Float64, m.captures[1])
        return nu_min, Inf
    end

    throw(ArgumentError("No match found for 'ν > ...' or '... < ν < ...' in model: $(model.name)"))
end

function get_thresholds(M::fModel, P::Parameters)
    thresholds_dict = P.thresholds_dict
    if !isnothing(M.thresholds_dict)
        thresholds_dict = merge(thresholds_dict, M.thresholds_dict)
    end

    thresholds = Vector{Float64}(undef, length(M.terms))

    i_core = 0
    for (i, term) in enumerate(M.terms)
        if !M.core[i]
            threshold_keys = [key for key in keys(thresholds_dict) if isa(key, String) && occursin(key, term)]
            if length(threshold_keys) != 1
                error("get_thresholds: no or multiple threshold found for term $term")
            end
            thresholds[i] = thresholds_dict[threshold_keys[1]]
            continue
        end

        # is_core
        i_core += 1
        outer_channel = M.outer_channels.i[i_core]
        lc = outer_channel.lc
        Jc = outer_channel.Jc

        if isa(outer_channel, fjQuantumNumbers)
            Fc = outer_channel.Fc
            Fc = Fc % 1 == 0 ? Int(Fc) : Fc
            key = coreQuantumNumbers(lc, Jc, Fc)
            if haskey(thresholds_dict, key)
                thresholds[i] = thresholds_dict[key]
                continue
            end
        end

        key = coreQuantumNumbers(lc, Jc)
        if haskey(thresholds_dict, key)
            thresholds[i] = thresholds_dict[key]
            continue
        end

        error("get_thresholds: missing threshold for key='$key' (term=$(term))")
    end

    return thresholds
end

function get_thresholds(M::kModel, P::Parameters)
    thresholds_dict = P.thresholds_dict
    if !isnothing(M.thresholds_dict)
        thresholds_dict = merge(thresholds_dict, M.thresholds_dict)
    end

    thresholds = Vector{Float64}(undef, length(M.terms))

    for (i, term) in enumerate(M.terms)
        outer_channel = M.jjchannels.i[i]
        lc = outer_channel.lc
        Jc = outer_channel.Jc

        key = coreQuantumNumbers(lc, Jc)
        if haskey(thresholds_dict, key)
            thresholds[i] = thresholds_dict[key]
            continue
        end

        key = coreQuantumNumbers(lc, NaN)
        if haskey(thresholds_dict, key)
            thresholds[i] = thresholds_dict[key]
            continue
        end

        error("get_thresholds: missing threshold for key='$key' (term=$(term))")
    end

    return thresholds
end

function get_species_parameters(species::Symbol)
    species_module = getfield(MQDT, species)
    for nm in names(species_module; all=true)
        obj = getfield(species_module, nm)
        if isa(obj, Parameters) && obj.species == species
            return obj
        end
    end
    return error("No Parameters found for species: $species")
end

"""
    get_fmodels(species::Symbol, lr::Int, Jr::Float64, Fc::Float64, F::Float64, param::Parameters)

Given atomic `species` and quantum numbers of a Rydberg channel (`lr`, `Jr`, `Fc`, `F`),
this function searches for all existing `fModel` in the species module, which includes the specified channel.
If no matching `fModel` is found, a single-channel `fModel` is generated and returned.
"""
function get_fmodels(species::Symbol, lr::Int, Jr::Float64, Fc::Float64, F::Float64, param::Parameters)
    ic = param.spin

    if ic == 0
        qn = jjQuantumNumbers(0.5, 0, 0.5, lr, Jr, F)
    else
        qn = fjQuantumNumbers(0.5, 0, 0.5, Fc, lr, Jr, F)
    end

    species_module = getfield(MQDT, species)
    models = Vector{fModel}()
    for nm in names(species_module; all=true)
        obj = getfield(species_module, nm)
        if !isa(obj, fModel) || (obj.species != species) || (obj.F != F)
            continue
        end
        for channel in obj.outer_channels.i
            if channel == qn
                push!(models, obj)
            end
        end
    end
    if !isempty(models)
        return models
    end

    # else generate single-channel model
    return [single_channel_model(species, lr, Jr, Fc, F, param)]
end

function single_channel_model(species::Symbol, lr::Int, Jr::Float64, Fc::Float64, F::Float64, param::Parameters)
    sr = 1 / 2
    Jc = 1 / 2
    ic = param.spin
    if !(abs(lr - sr) <= Jr <= (lr + sr)) || !(abs(Jc - ic) <= Fc <= (Jc + ic)) || !(abs(Jr - Fc) <= F <= (Jr + Fc))
        error("Quantum numbers not allowed: lr=$lr, Jr=$Jr, Fc=$Fc, F=$F")
    end

    if lr <= 0
        error("expected only lr > 0: lr=$lr, Jr=$Jr, Fc=$Fc, F=$F")
    end

    if param.spin == 0
        channel = jjChannels([jjQuantumNumbers(0.5, 0, 0.5, lr, Jr, F, F)])
        name = "SQDT L=$lr, Jr=$Jr, J=$F, ν > $(lr+1)"
        term = "JJ: L=$lr, Jr=$Jr, J=$F"
    else
        channel = fjChannels([fjQuantumNumbers(0.5, 0, 0.5, Fc, lr, Jr, F)])
        name = "SQDT L=$lr, Jr=$Jr, F=$F, ν > $(lr+1)"
        term = "FJ: I=$(param.spin), L=$lr, Jr=$Jr, F=$F"
    end

    return fModel(species, name, F, (lr+1, Inf), [term], Bool[1], [0;;], [""], [0;;], channel, channel, [1;;])
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
        nu_list::Vector{Float64},
        lr_list::Vector{Union{Int, Nothing}},
        coefficients::Vector{Float64},
        channels::Channels,
        term::String,
        lead::Float64,
        model::fModel,
    )

Type to store all relevant information of multi-channel bound states for a given Rydberg series.
`BasisState` is generated by the `basisarray` function using a `Model` and `EigenStates`.
"""
struct BasisState
    species::Symbol
    energy::Float64
    parity::Int
    f::Float64
    nu_list::Vector{Float64}
    lr_list::Vector{Union{Int,Nothing}}
    coefficients::Vector{Float64}
    channels::Channels
    term::String
    lead::Float64
    model::fModel
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
