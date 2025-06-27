module Yb174

using ..MQDT:
    Parameters,
    fModel,
    lsChannels,
    jjChannels,
    lsQuantumNumbers,
    jjQuantumNumbers,
    test_model

export PARA,
    FMODEL_HIGHN_S0,
    FMODEL_HIGHN_S1,
    FMODEL_LOWN_P0,
    FMODEL_HIGHN_P0,
    FMODEL_HIGHN_P1,
    FMODEL_HIGHN_P2,
    FMODEL_HIGHN_D1,
    FMODEL_HIGHN_D2,
    FMODEL_HIGHN_D3

# Isotope data
PARA = Parameters(
    :Yb174,
    1822.88848628*173.9388621, # nuclear mass
    0, # nuclear spin
    109736.9695858, # Rydberg constant in 1/cm
    50443.070393, # lowest ionization threshold in 1/cm
    0, # hyperfine constant in 1/cm
    2.1, # nuclear dipole
)

# MQDT Models
FMODEL_HIGHN_S0 = fModel(
    "S J=0, ν > 2", # fit for states 6s7s upward [Phys. Rev. X 15, 011009 (2025)]
    6,
    [
        "6sns 1S0",
        "4f13 5d 6snl a",
        "6pnp 1S0",
        "4f13 5d 6snl b",
        "6pnp 3P0",
        "4f13 5d 6snl c",
    ],
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
        jjQuantumNumbers(0.5, 1, 1.5, 1, 1.5, 0),
        jjQuantumNumbers(0.5, 1, 0.5, 1, 0.5, 0),
    ]),
    [
        1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 -sqrt(2/3) 0 sqrt(1/3) 0;
        0 0 0 1 0 0;
        0 0 sqrt(1/3) 0 sqrt(2/3) 0;
        0 0 0 0 0 1
    ], # according to paper & rydcalc code
    #[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 sqrt(2/3) 0 sqrt(1/3) 0; 0 0 0 1 0 0; 0 0 -sqrt(1/3) 0 sqrt(2/3) 0; 0 0 0 0 0 1] # correct frame transformation
)

FMODEL_HIGHN_S1 = fModel(
    "S J=1, ν > 26", # fit only valid from 28s upward [Phys. Rev. Lett. 128, 033201 (2022)]
    1,
    ["6sns 3S1"],
    Bool[1],
    [50443.070393],
    #[4.4382 4 -1e4 8e6 -3e9], # rydberg ritz
    [0.4382 4 -1e4 8e6 -3e9], # adjusted for non-rydberg ritz
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 1)]),
    [1;;],
)

FMODEL_LOWN_P0 = fModel(
    "P J=0, 1.5 < ν < 5.5", # fit to NIST data between 6p and 9p
    1,
    ["6snp 3P0"],
    Bool[1],
    [50443.070393],
    [0.969279 0.288219 1.36228],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 0)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 0)]),
    [1;;],
)

FMODEL_HIGHN_P0 = fModel(
    "P J=0, ν > 6", # [Phys. Rev. X 15, 011009 (2025)]
    2,
    ["6snp 3P0", "4f13 5d 6snl"],
    Bool[1, 0],
    [50443.070393, 83967.7],
    [0.95356884 -0.28602498; 0.19845903 0],
    ["12"],
    [0.16328854 0],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 0)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 0)]),
    [1 0; 0 1],
)

FMODEL_HIGHN_P1 = fModel(
    "P J=1, ν > 6", # [Phys. Rev. X 15, 011009 (2025)]
    6,
    [
        "6snp 1P1",
        "6snp 3P1",
        "4f13 5d 6snl a",
        "4f13 5d 6snl b",
        "4f13 5d 6snl c",
        "4f13 5d 6snl d",
    ],
    Bool[1, 1, 0, 0, 0, 0],
    [50443.070393, 50443.070393, 83967.7, 83967.7, 83967.7, 83967.7],
    [
        0.92271098 2.6036257;
        0.98208719 -5.4562725;
        0.22851720 0;
        0.20607759 0;
        0.19352751 0;
        0.18153094 0
    ],
    ["12", "13", "14", "15", "16", "23", "24", "25", "26"],
    [
        -0.08410871 120.37555 -9314.23;
        -0.07318156 0 0;
        -0.06651977 0 0;
        -0.02210989 0 0;
        -0.10451698 0 0;
        0.02477048 0 0;
        0.05765807 0 0;
        0.08606276 0 0;
        0.04994363 0 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 1),
        jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 1),
    ]),
    [
        sqrt(2/3) sqrt(1/3) 0 0 0 0;
        -sqrt(1/3) sqrt(2/3) 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1
    ],
)

FMODEL_HIGHN_P2 = fModel(
    "P J=2, ν > 5", # [Phys. Rev. X 15, 011009 (2025)]
    4,
    ["6snp 3P2", "4f13 5d 6snl a", "4f13 5d 6snl b", "4f13 5d 6snl c"],
    Bool[1, 0, 0, 0],
    [50443.070393, 83967.7, 83967.7, 83967.7],
    [0.925121305 -2.73247165 74.664989; 0.230133261 0 0; 0.209638118 0 0; 0.186228192 0 0],
    ["12", "13", "14"],
    [0.0706666127 0; 0.0232711158 0; -0.0292153659 0],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 2)]),
    [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],
)

FMODEL_HIGHN_D1 = fModel(
    "D J=1, ν > 26", # fit only valid from 30d upward [Phys. Rev. X 15, 011009 (2025)]
    1,
    ["6snd 3D1"],
    Bool[1],
    [50443.070393],
    [0.75258093 0.3826 -483.1],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 1)]),
    [1;;],
)

FMODEL_HIGHN_D2 = fModel(
    "D J=2, ν > 5", # [Phys. Rev. X 15, 011009 (2025)]
    5,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2"],
    Bool[1, 1, 0, 0, 1],
    [50443.070393, 50443.070393, 83967.7, 83967.7, 79725.35],
    [
        0.729500971 -0.0284447537;
        0.75229161 0.0967044398;
        0.196120406 0;
        0.233821165 0;
        0.152890218 0
    ],
    ["12", "13", "14", "24", "15", "25"],
    [
        0.21157531 -15.3844;
        0.00522534111 0;
        0.0398754262 0;
        -0.00720265975 0;
        0.104784389 0;
        0.0721775002 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 2),
        jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 2),
        jjQuantumNumbers(0.5, 1, 0.5, 1, 1.5, 2), # Jc and Jr could also be switched or both be 3/2
    ]),
    [
        sqrt(3/5) sqrt(2/5) 0 0 0;
        -sqrt(2/5) sqrt(3/5) 0 0 0;
        0 0 1 0 0;
        0 0 0 1 0;
        0 0 0 0 1
    ],
)

FMODEL_HIGHN_D3 = fModel(
    "D J=3, ν > 14", # fit only valid from 30d upward [Phys. Rev. X 15, 011009 (2025)], provides good match around 18d
    1,
    ["6snd 3D3"],
    Bool[1],
    [50443.070393],
    [0.72895315 -0.2065 220.5],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 3)]),
    [1;;],
)

MODELS = [
    FMODEL_HIGHN_S0,
    FMODEL_HIGHN_S1,
    FMODEL_LOWN_P0,
    FMODEL_HIGHN_P0,
    FMODEL_HIGHN_P1,
    FMODEL_HIGHN_P2,
    FMODEL_HIGHN_D1,
    FMODEL_HIGHN_D2,
    FMODEL_HIGHN_D3,
]
test_model(MODELS)

end
