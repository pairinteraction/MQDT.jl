# --------------------------------------------------------
# Store Parameters
# --------------------------------------------------------

"""
    Parameters(species::Symbol, mass::Float64, spin::Float64, rydberg::Float64, threshold::Float64, hyperfine::Float64, dipole::Float64)

Type to store relevant parameters for each atomic element. A corresponding 
file in the parameters folder of the mqdt module initializes this type.

# Examples

```julia-repl
PARA = Parameters(
    :Yb174,
    1822.88848628*173.9388621, # nuclear mass
    0, # nuclear spin
    109736.9695858, # Rydberg constant in 1/cm
    50443.070393, # lowest ionization threshold in 1/cm
    0, # hyperfine constant in 1/cm
    2.1 # nuclear dipole
)
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
# Channel structs
# --------------------------------------------------------

"""
See also [`lsChannels`](@ref), [`jjChannels`](@ref), [`fjChannels`](@ref)

Abstract type to store the quantum numbers of an MQDT channel representation. 
"""
abstract type Channels end

"""
See also [`jjChannels`](@ref), [`fjChannels`](@ref)

    lsChannels(sc::Vector{Float64}, sr::Vector{Float64}, S::Vector{Float64}, lc::Vector{Int}, lr::Vector{Int}, L::Vector{Int}, J::Vector{Float64})

Type to store the quantum numbers of an MQDT channel representation, here in the "LS" coupling scheme, 
i.e. core and Rydberg angular momentum quantum numbers are coupled individually to give a total spin (S) and total orbital angular momentum (L).

# Examples

```julia-repl
mqdt.lsChannels(
    [0.5, 0.5],
    [0.5, 0.5],
    [0., 1.],
    [0, 0],
    [0, 0],
    [0, 0],
    [0., 1.]
)
```
"""
struct lsChannels <: Channels
    sc::Vector{Float64}
    sr::Vector{Float64}
    S::Vector{Float64}
    lc::Vector{Int}
    lr::Vector{Int}
    L::Vector{Int}
    J::Vector{Float64}
end

"""
See also [`lsChannels`](@ref), [`fjChannels`](@ref)

    jjChannels(sc::Vector{Float64}, lc::Vector{Int}, Jc::Vector{Float64}, sr::Vector{Float64}, lr::Vector{Int}, Jr::Vector{Float64}, J::Vector{Float64})

Type to store the quantum numbers of an MQDT channel representation, here in the 
"jj" coupling scheme, i.e. core and Rydberg angular momentum (AM) quantum numbers 
are coupled separately to give a total core AM (Jc) and total Rydberg AM (Jr).

# Examples

```julia-repl
mqdt.jjChannels(
    [0.5, 0.5],
    [0, 0],
    [0.5, 0.5],
    [0.5, 0.5],
    [0, 0],
    [0.5, 0.5],
    [0., 1.]
)
```
"""
struct jjChannels <: Channels
    sc::Vector{Float64}
    lc::Vector{Int}
    Jc::Vector{Float64}
    sr::Vector{Float64}
    lr::Vector{Int}
    Jr::Vector{Float64}
    J::Vector{Float64}
end

"""
See also [`lsChannels`](@ref), [`jjChannels`](@ref)

    fjChannels(sc::Vector{Float64}, lc::Vector{Int}, Jc::Vector{Float64}, Fc::Vector{Float64}, sr::Vector{Float64}, lr::Vector{Int}, Jr::Vector{Float64}, F::Vector{Float64})

Type to store the quantum numbers of an MQDT channel representation, here in the "fj" (sometimes called 
"fragmentation") coupling scheme, i.e. core angular momentum quantum numbers are coupled to the nuclear spin 
to give a core hyperfine configuration (Fc) which is coupled to the total Rydberg angular momentum (Jr).

# Examples

```julia-repl
mqdt.fjChannels(
    [0.5, 0.5, 0.5],
    [0, 0, 0],
    [0.5, 0.5, 0.5],
    [0., 1., 1.],
    [0.5, 0.5, 0.5],
    [0, 0, 0],
    [0.5, 0.5, 0.5],
    [0.5, 0.5, 1.5],
)
```
"""
struct fjChannels <: Channels
    sc::Vector{Float64}
    lc::Vector{Int}
    Jc::Vector{Float64}
    Fc::Vector{Float64}
    sr::Vector{Float64}
    lr::Vector{Int}
    Jr::Vector{Float64}
    F::Vector{Float64}
end

"""
See also [`fjChannels`](@ref), [`fj_quantum_numbers`](@ref)

    fjQuantumNumbers(sc::Float64, lc::Int, Jc::Float64, Fc::Float64, sr::Float64, lr::Int, Jr::Float64, F::Float64)

Type to store the quantum numbers of a single MQDT channel in the "fj" (sometimes called "fragmentation") coupling scheme.
Can conveniently be constructed using from [`fjChannels`](@ref) using the [`fj_quantum_numbers`](@ref) function.

# Examples

```julia-repl
mqdt.fjQuantumNumbers(
    0.5,
    0,
    0.5,
    1.,
    0.5,
    0,
    0.5,
    1.5,
)
```
"""
struct fjQuantumNumbers <: Channels
    sc::Float64
    lc::Int
    Jc::Float64
    Fc::Float64
    sr::Float64
    lr::Int
    Jr::Float64
    F::Float64
end

# --------------------------------------------------------
# Base functions for Channel structs
# --------------------------------------------------------

function Base.size(T::Channels)
    return length(T.lr)
end

function Base.cat(T::lsChannels)
    return vcat(T.sc', T.sr', T.S', T.lc', T.lr', T.L', T.J')
end

function Base.cat(T::jjChannels)
    return vcat(T.sc', T.lc', T.Jc', T.sr', T.lr', T.Jr', T.J')
end

function Base.cat(T::fjChannels)
    return vcat(T.sc', T.lc', T.Jc', T.Fc', T.sr', T.lr', T.Jr', T.F')
end

"""
See also [`unique_parity`](@ref)

    parity(T::Channels)

Given a channel representation, this function returns the parity of a each 
channel based of the orbital angular momenta of the core and the Rydberg electron.
"""
function parity(T::Channels)
    pc = 2iseven.(T.lc) .- 1
    pr = 2iseven.(T.lr) .- 1
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
    good_quantum_number(T::jjChannels)
    good_quantum_number(T::fjChannels)

Given a jj (fj) channel representation, this function checks whether each channel 
has the same total (hyperfine) angular momentum J (F) and, if true, returns it.
"""
function good_quantum_number(T::jjChannels)
    f = unique(T.J)
    if length(f) == 1
        return f[1]
    else
        error("good_quantum_number: channels have different J")
    end
end

function good_quantum_number(T::fjChannels)
    f = unique(T.F)
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
See also [`Model`](@ref)

Abstract type to store MQDT models. 
"""
abstract type Models end

"""
    Model(
        name::String,
        size::Int,
        lsterms::Vector{String},
        # jjterms::Vector{String},
        core::Vector{Bool},
        rydbergritz::Vector{Bool},
        thresholds::Vector{Float64},
        defects::Matrix{Float64},
        mixing::Vector{String},
        angles::Matrix{Float64},
        SLJ::Matrix{Int},
        channels::Channels,
        unitary::Matrix{Float64}
    )

Type to store MQDT models. Contains all relevant parameters to form the K matrix and calculate the spectrum for a specific Rydberg series.

# Examples

```julia-repl
MODEL_S0 = Model(
    "S0",
    6, 
    ["6sns 1S0", "4f13 5d 6snp a", "6pnl 1S0", "4f13 5d 6snp b", "6pnl 3P0", "4f13 5d 6snp c"],
    # ["(6s1/2)(ns1/2)", "5d a", "(6p3/2)(np3/2)", "5d b", "(6p1/2)(np1/2)", "5d c"],
    Bool[1, 0, 1, 0, 1, 0], 
    Bool[0, 0, 0, 0, 0, 0], 
    [50443.070393, 83967.7, 80835.39, 83967.7, 77504.98, 83967.7],
    [0.355097325 0.278368431; 0.204537279 0; 0.116394359 0; 0.295432196 0; 0.25765161 0; 0.155807042 0],
    ["12", "13", "14", "34", "35", "16"],
    [0.12654859 0; 0.30010744 0; 0.05703381 0; 0.11439805 0; 0.09864375 0; 0.14248210 0],
    [0 0 0; 0 0 0; 1 1 0],
    mqdt.jjChannels([0.5, 0.5, 0.5], [0, 1, 1], [0.5, 0.5, 0.5], [0.5, 0.5, 0.5], [0, 1, 1], [0.5, 1.5, 0.5], [0, 0, 0]),
    [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 -sqrt(2/3) 0 sqrt(1/3) 0; 0 0 0 1 0 0; 0 0 sqrt(1/3) 0 sqrt(2/3) 0; 0 0 0 0 0 1]
)
```
"""
struct Model <: Models
    name::String
    size::Int
    lsterms::Vector{String}
    core::Vector{Bool}
    rydbergritz::Vector{Bool}
    thresholds::Vector{Float64}
    defects::Matrix{Float64}
    mixing::Vector{String}
    angles::Matrix{Float64}
    SLJ::Matrix{Int}
    channels::Channels
    unitary::Matrix{Float64}
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

function get_lr(T::Model)
    l = zeros(Int, T.size)
    l[findall(T.core)] = T.channels.lr
    return l
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

    BasisState(energy::Float64, parity::Int, f::Float64, nu::Vector{Float64}, lr::Vector{Int}, S::Vector{Float64}, coeff::Vector{Float64}, channels::Channels)

Type to store all relevant information of multi-channel bound states for a given Rydberg series. 
`BasisState` is generated by the `basisarray` function using a `Model` and `EigenStates`.
"""
struct BasisState
    energy::Float64
    parity::Int
    f::Float64
    nu::Vector{Float64}
    lr::Vector{Int}
    S::Vector{Float64}
    coeff::Vector{Float64}
    channels::Channels
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

    DataBaseState(nu::Float64, parity::Int, f::Float64, nui::Vector{Float64}, nui_all::Vector{Float64}, ai::Vector{Float64}, ai_all::Vector{Float64}, 
        S::Vector{Float64}, L::Vector{Float64}, J::Vector{Float64}, lr::Vector{Int}, Jr::Vector{Float64}, neg::Float64)

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
    L::Vector{Float64}
    J::Vector{Float64}
    lr::Vector{Int}
    Jr::Vector{Float64}
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

function Base.size(T::DataBaseArray)
    return length(T.states)
end

function Base.union(T1::DataBaseArray, T2::DataBaseArray)
    return DataBaseArray(vcat(T1.states, T2.states))
end
