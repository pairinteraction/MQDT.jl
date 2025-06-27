module Yb171

using ..MQDT:
    Parameters,
    fModel,
    lsChannels,
    fjChannels,
    lsQuantumNumbers,
    fjQuantumNumbers,
    test_model

export PARA,
    FMODEL_HIGHN_S05,
    FMODEL_HIGHN_S15,
    FMODEL_LOWN_P05,
    FMODEL_HIGHN_P05,
    FMODEL_HIGHN_P15,
    FMODEL_HIGHN_P25,
    FMODEL_HIGHN_D05,
    FMODEL_HIGHN_D15,
    FMODEL_HIGHN_D25,
    FMODEL_HIGHN_D35

# Isotope data
PARA = Parameters(
    :Yb171,
    1822.88848628*170.9363258, # nuclear mass
    0.5, # nuclear spin
    109736.9635066, # Rydberg constant in 1/cm
    # 50442.795744, # lowest ionization threshold in 1/cm (Fc=0)
    50443.217463, # higher hyperfine ionization threshold in 1/cm (Fc=1)
    0.421719, # hyperfine constant in 1/cm
    # 12.64 # hyperfine constant in GHz
    0.49367, # nuclear dipole
)

# MQDT Models
FMODEL_HIGHN_S05 = fModel(
    "S F=1/2, ν > 26", # [Phys. Rev. X 15, 011009 (2025)]
    7,
    [
        "6sns 1S0",
        "4f13 5d 6snl a",
        "6pnp 1S0",
        "4f13 5d 6snl b",
        "6pnp 3P0",
        "4f13 5d 6snl c",
        "6sns 3S1",
    ],
    Bool[1, 0, 1, 0, 1, 0, 1],
    [50442.795744, 83967.7, 80835.39, 83967.7, 77504.98, 83967.7, 50443.217463],
    [
        0.357519763 0.298712849 0 0 0;
        0.203907536 0 0 0 0;
        0.116803536 0 0 0 0;
        0.286731074 0 0 0 0;
        0.248113946 0 0 0 0;
        0.148678953 0 0 0 0;
        0.438426851 3.91762642 -10612.6828 8017432.38 -2582622910.0
    ],
    ["12", "13", "14", "34", "35", "16"],
    [
        0.131810463 0;
        0.297612147 0;
        0.055508821 0;
        0.101030515 0;
        0.102911159 0;
        0.137723736 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 0, 0, 0),
        lsQuantumNumbers(0.5, 0, 1, 1, 0, 0),
        lsQuantumNumbers(0.5, 1, 1, 1, 1, 0),
        lsQuantumNumbers(0.5, 1, 0, 0, 0, 1),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 0, 0, 0.5, 0.5),
        fjQuantumNumbers(0.5, 1, 1.5, 1, 1, 1.5, 0.5),
        fjQuantumNumbers(0.5, 1, 0.5, 0, 1, 0.5, 0.5), # Fc of the 3P0 state could be both 0 or 1
        fjQuantumNumbers(0.5, 0, 0.5, 1, 0, 0.5, 0.5),
    ]),
    [
        1/2 0 0 0 0 0 sqrt(3)/2;
        0 1 0 0 0 0 0;
        0 0 -sqrt(2/3) 0 sqrt(1/3) 0 0;
        0 0 0 1 0 0 0;
        0 0 sqrt(1/3) 0 sqrt(2/3) 0 0;
        0 0 0 0 0 1 0;
        sqrt(3)/2 0 0 0 0 0 -1/2
    ], # according to paper & rydcalc code
    #[1/2 0 0 0 0 0 sqrt(3)/2; 0 1 0 0 0 0 0; 0 0 sqrt(2/3) 0 sqrt(1/3) 0 0; 0 0 0 1 0 0 0; 0 0 -sqrt(1/3) 0 sqrt(2/3) 0 0; 0 0 0 0 0 1 0; sqrt(3)/2 0 0 0 0 0 -1/2] # correct frame transformation
)

FMODEL_HIGHN_S15 = fModel(
    "S F=3/2, ν > 26", # [Phys. Rev. X 15, 011009 (2025)]
    1,
    ["6sns 3S1"],
    Bool[1],
    [50443.217463],
    [0.438426851 3.91762642 -10612.6828 8017432.38 -2582622910.0],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 0, 0.5, 1.5)]),
    [1;;],
)

FMODEL_LOWN_P05 = fModel(
    "P F=1/2, 1.5 < ν < 2.5", # fit to Yb174 NIST data
    3,
    ["6snp 1P1", "6snp 3P1", "6snp 3P0"],
    Bool[1, 1, 1],
    [50443.217463, 50443.217463, 50442.795744],
    [0.161083 0; 0.920424 0; 0.180701 0],
    ["12"],
    [-0.426128 6.272986],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 0.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 0.5, 0.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 1, 0.5, 0.5),
    ]),
    [-sqrt(2/3) -sqrt(1/3) 0; 1/(2sqrt(3)) -sqrt(1/6) sqrt(3)/2; -1/2 1/sqrt(2) 1/2],
)

FMODEL_HIGHN_P05 = fModel(
    "P F=1/2, ν > 5.7", # [Phys. Rev. X 15, 011009 (2025)] # fit for ν > 28, but extrapolates nicely down for ν > 5.7
    8,
    [
        "6snp 1P1",
        "6snp 3P1",
        "4f13 5d 6snl a",
        "4f13 5d 6snl b",
        "4f13 5d 6snl c",
        "4f13 5d 6snl d",
        "6snp 3P0",
        "4f13 5d 6snl e",
    ],
    Bool[1, 1, 0, 0, 0, 0, 1, 0],
    [50443.217463, 50443.217463, 83967.7, 83967.7, 83967.7, 83967.7, 50442.795744, 83967.7],
    [
        0.921706585 2.6036257;
        0.979638580 -5.4562725;
        0.228828720 0;
        0.205484818 0;
        0.193528629 0;
        0.181385000 0;
        0.953071282 0.131025247;
        0.198445928 0
    ],
    ["12", "27", "13", "14", "15", "16", "23", "24", "25", "26", "78"],
    [
        -0.087127227 135.400009 -12985.0162;
        -0.001430175 0 0;
        -0.073904060 0 0;
        -0.063632668 0 0;
        -0.021924569 0 0;
        -0.106678810 0 0;
        0.032556999 0 0;
        0.054105142 0 0;
        0.086127672 0 0;
        0.053804487 0 0;
        0.163043619 0 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 0.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 0.5, 0.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 1, 0.5, 0.5),
    ]),
    [
        -sqrt(2/3) -sqrt(1/3) 0 0 0 0 0 0;
        1/(2sqrt(3)) -sqrt(1/6) 0 0 0 0 sqrt(3)/2 0;
        0 0 1 0 0 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        -1/2 1/sqrt(2) 0 0 0 0 1/2 0;
        0 0 0 0 0 0 0 1
    ],
)

FMODEL_HIGHN_P15 = fModel(
    "P F=3/2, ν > 18", # [Phys. Rev. X 15, 011009 (2025)] # fit for ν > 28, but extrapolates nicely down for ν > 18
    10,
    [
        "6snp 1P1",
        "6snp 3P1",
        "4f13 5d 6snl a",
        "4f13 5d 6snl b",
        "4f13 5d 6snl c",
        "4f13 5d 6snl d",
        "6snp 3P2",
        "4f13 5d 6snl e",
        "4f13 5d 6snl f",
        "4f13 5d 6snl g",
    ],
    Bool[1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
    [
        50443.217463,
        50443.217463,
        83967.7,
        83967.7,
        83967.7,
        83967.7,
        50442.795744,
        83967.7,
        83967.7,
        83967.7,
    ],
    [
        0.921706585 2.56569459 0;
        0.979638580 -5.239904224 0;
        0.228828720 0 0;
        0.205484818 0 0;
        0.193528629 0 0;
        0.181385000 0 0;
        0.924825736 -3.542481644 81.5334687;
        0.236866903 0 0;
        0.221055883 0 0;
        0.185599376 0 0
    ],
    ["12", "13", "14", "15", "16", "23", "24", "25", "26", "78", "79", "710"],
    [
        -0.087127227 135.400009 -12985.0162;
        -0.073904060 0 0;
        -0.063632668 0 0;
        -0.021924569 0 0;
        -0.106678810 0 0;
        0.032556999 0 0;
        0.054105142 0 0;
        0.086127672 0 0;
        0.053804487 0 0;
        0.071426685 0 0;
        0.027464110 0 0;
        -0.029741862 0 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 0.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 1, 1.5, 1.5),
    ]),
    [
        sqrt(5/3)/2 sqrt(5/6)/2 0 0 0 0 -sqrt(3/2)/2 0 0 0;
        -sqrt(1/3) sqrt(2/3) 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0;
        1/2 1/(2sqrt(2)) 0 0 0 0 sqrt(5/2)/2 0 0 0;
        0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 1
    ],
)

FMODEL_HIGHN_P25 = fModel(
    "P F=5/2, ν > 28", # [Phys. Rev. X 15, 011009 (2025)]
    4,
    ["6snp 3P2", "4f13 5d 6snl a", "4f13 5d 6snl b", "4f13 5d 6snl c"],
    Bool[1, 0, 0, 0],
    [50443.217463, 83967.7, 83967.7, 83967.7],
    [
        0.924825736 -3.542481644 81.5334687;
        0.236866903 0 0;
        0.221055883 0 0;
        0.185599376 0 0
    ],
    ["12", "13", "14"],
    [0.071426685 0; 0.027464110 0; -0.029741862 0],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 2.5)]),
    [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],
)

FMODEL_HIGHN_D05 = fModel(
    "D F=1/2, ν > 30", # [Phys. Rev. X 15, 011009 (2025)]
    1,
    ["6snd 3D1"],
    Bool[1],
    [50443.21746],
    [0.75258093 0.382628525 -483.120633],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 1)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 0.5)]),
    [1;;],
)

FMODEL_HIGHN_D15 = fModel(
    "D F=3/2, ν > 30", # [Phys. Rev. X 15, 011009 (2025)]
    6,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2", "6snd 3D1"],
    Bool[1, 1, 0, 0, 1, 1],
    [50443.21746, 50443.21746, 83967.7, 83967.7, 79725.35, 50442.795744],
    [
        0.730537124 -0.000186828866 0;
        0.751591782 -0.00114049637 0;
        0.196120394 0 0;
        0.233742396 0 0;
        0.152905343 0 0;
        0.75258093 0.382628525 -483.120633
    ],
    ["12", "13", "14", "24", "15", "25"],
    [
        0.205496654 0;
        0.00522401624 0;
        0.0409502343 0;
        -0.00378075773 0;
        0.108563952 0;
        0.0665700438 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 1.5),
        fjQuantumNumbers(0.5, 1, 1.5, 1, 1, 1.5, 1.5), # Fc of the 1D2 state could be both 0 or 1
        fjQuantumNumbers(0.5, 0, 0.5, 0, 2, 1.5, 1.5),
    ]),
    [
        -sqrt(3/5) -sqrt(2/5) 0 0 0 0;
        sqrt(3/5)/2 -3/(2sqrt(10)) 0 0 0 sqrt(5/2)/2;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        -1/2 sqrt(3/2)/2 0 0 0 sqrt(3/2)/2
    ],
)

FMODEL_HIGHN_D25 = fModel(
    "D F=5/2, ν > 30", # [Phys. Rev. X 15, 011009 (2025)]
    6,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2", "6snd 3D3"],
    Bool[1, 1, 0, 0, 1, 1],
    [50443.21746, 50443.21746, 83967.7, 83967.7, 79725.35, 50442.795744],
    [
        0.730537124 -0.000186828866 0;
        0.751591782 -0.00114049637 0;
        0.196120394 0 0;
        0.233742396 0 0;
        0.152905343 0 0;
        0.72895315 -0.20653489 220.484722
    ],
    ["12", "13", "14", "24", "15", "25"],
    [
        0.205496654 0;
        0.00522401624 0;
        0.0409502343 0;
        -0.00378075773 0;
        0.108563952 0;
        0.0665700438 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 2.5),
        fjQuantumNumbers(0.5, 1, 1.5, 1, 1, 1.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 2, 2.5, 2.5),
    ]),
    [
        sqrt(7/5)/2 sqrt(7/30) 0 0 0 -sqrt(5/3)/2;
        -sqrt(2/5) sqrt(3/5) 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        1/2 sqrt(1/6) 0 0 0 sqrt(7/3)/2
    ],
)

FMODEL_HIGHN_D35 = fModel(
    "D F=7/2, ν > 14", # fit only valid from 30d upward [Phys. Rev. X 15, 011009 (2025)], provides good match around 18d
    1,
    ["6snd 3D3"],
    Bool[1],
    [50443.21746],
    [0.72895315 -0.20653489 220.484722],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 3.5)]),
    [1;;],
)

MODELS = [
    FMODEL_HIGHN_S05,
    FMODEL_HIGHN_S15,
    FMODEL_LOWN_P05,
    FMODEL_HIGHN_P05,
    FMODEL_HIGHN_P15,
    FMODEL_HIGHN_P25,
    FMODEL_HIGHN_D05,
    FMODEL_HIGHN_D15,
    FMODEL_HIGHN_D25,
    FMODEL_HIGHN_D35,
]
test_model(MODELS)

end
