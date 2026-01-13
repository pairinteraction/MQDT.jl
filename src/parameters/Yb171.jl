module Yb171

using ..MQDT:
    Parameters,
    fModel,
    lsChannels,
    jjChannels,
    fjChannels,
    lsQuantumNumbers,
    jjQuantumNumbers,
    fjQuantumNumbers,
    coreQuantumNumbers

export PARA,
    FMODEL_HIGHN_S05,
    FMODEL_HIGHN_S15,
    FMODEL_HIGHN_P05,
    FMODEL_HIGHN_P15,
    FMODEL_HIGHN_D05,
    FMODEL_HIGHN_D15,
    FMODEL_HIGHN_D25,
    FMODEL_HIGHN_D35,
    FMODEL_HIGHN_F25,
    FMODEL_HIGHN_F35,
    FMODEL_HIGHN_F45,
    FMODEL_HIGHN_G25,
    FMODEL_HIGHN_G35,
    FMODEL_HIGHN_G45,
    FMODEL_HIGHN_G55,
    FMODEL_LOWN_S05,
    FMODEL_LOWN_S15,
    FMODEL_LOWEST_P05,
    FMODEL_LOWN_P05,
    FMODEL_LOWEST_P15,
    FMODEL_LOWN_P15,
    FMODEL_LOWEST_P25,
    FMODEL_LOWN_P25,
    FMODEL_LOWN_D05,
    FMODEL_LOWN_D15,
    FMODEL_LOWN_D25,
    FMODEL_LOWN_D35

# Isotope data
THRESHOLDS = Dict(
    coreQuantumNumbers(0, 0.5, 0) => 50442.795744,
    coreQuantumNumbers(0, 0.5, 1) => 50443.217463,
    coreQuantumNumbers(1, 0.5, NaN) => 77504.98,
    coreQuantumNumbers(1, NaN, NaN) => 79725.35,
    coreQuantumNumbers(1, 1.5, NaN) => 80835.39,
    "4f13 5d 6s" => 83967.7,
)

PARA = Parameters(
    :Yb171,
    1822.88848628*170.9363258, # nuclear mass
    0.5, # nuclear spin
    109736.9635066, # Rydberg constant in 1/cm
    THRESHOLDS[coreQuantumNumbers(0, 0.5, 1)], # higher hyperfine ionization threshold in 1/cm (Fc=1)
    0.421719, # hyperfine constant in 1/cm
    # 12.64 # hyperfine constant in GHz
    0.49367, # nuclear dipole
    THRESHOLDS,
)

# --------------------------------------------------------
# MQDT models valid at large n
# --------------------------------------------------------

FMODEL_HIGHN_S05 = fModel(
    :Yb171,
    "S F=1/2, ν > 26", # [Phys. Rev. A 112, 042817 (2025), Phys. Rev. X 15, 011009 (2025)]
    7,
    ["6sns 1S0", "4f13 5d 6snl a", "6pnp 1S0", "4f13 5d 6snl b", "6pnp 3P0", "4f13 5d 6snl c", "6sns 3S1"],
    Bool[1, 0, 1, 0, 1, 0, 1],
    [
        0.357488757 0.165981371 0 0 0;
        0.203918644 0 0 0 0;
        0.116819032 0 0 0 0;
        0.287350241 0 0 0 0;
        0.247621114 0 0 0 0;
        0.148681324 0 0 0 0;
        0.438542187 3.78366407 -10709.7378 8054542.58 -2523011670
    ],
    ["1.2", "1.3", "1.4", "3.4", "3.5", "1.6"],
    [
        0.131755467 0;
        0.297504211 0;
        0.055421439 0;
        0.100871756 0;
        0.103123032 0;
        0.137753117 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 0, 0, 0, 0.5),
        lsQuantumNumbers(0.5, 0, 1, 1, 0, 0, 0.5),
        lsQuantumNumbers(0.5, 1, 1, 1, 1, 0, 0.5),
        lsQuantumNumbers(0.5, 1, 0, 0, 0, 1, 0.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 0, 0, 0.5, 0.5),
        fjQuantumNumbers(0.5, 1, 1.5, 1, 1, 1.5, 0.5),
        fjQuantumNumbers(0.5, 1, 0.5, NaN, 1, 0.5, 0.5), # Fc of the 3P0 state could be both 0 or 1
        fjQuantumNumbers(0.5, 0, 0.5, 1, 0, 0.5, 0.5),
    ]),
    [
        1/2 0 0 0 0 0 sqrt(3)/2;
        0 1 0 0 0 0 0;
        0 0 sqrt(2/3) 0 -sqrt(1/3) 0 0;
        0 0 0 1 0 0 0;
        0 0 sqrt(1/3) 0 sqrt(2/3) 0 0;
        0 0 0 0 0 1 0;
        sqrt(3)/2 0 0 0 0 0 -1/2
    ], # corrected in [arXiv:2507.11487v1]
)

FMODEL_HIGHN_S15 = fModel(
    :Yb171,
    "S F=3/2, ν > 26", # [Phys. Rev. X 15, 011009 (2025)]
    1,
    ["6sns 3S1"],
    Bool[1],
    [0.438426851 3.91762642 -10612.6828 8017432.38 -2582622910.0;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1, 1.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 0, 0.5, 1.5)]),
    [1;;],
)

FMODEL_HIGHN_P05 = fModel(
    :Yb171,
    "P F=1/2, ν > 5.7", # [Phys. Rev. A 112, 042817 (2025), Phys. Rev. X 15, 011009 (2025)] # fit for ν > 28, but extrapolates nicely down for ν > 5.7
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
    [
        0.922094502 2.12370136;
        0.981191543 -4.54209175;
        0.229094016 0;
        0.206073107 0;
        0.193527627 0;
        0.181165673 0;
        0.953185132 0.0277444042;
        0.198448494 0
    ],
    ["1.2", "2.7", "1.3", "1.4", "1.5", "1.6", "2.3", "2.4", "2.5", "2.6", "7.8"],
    [
        -0.102285383 153.521338 -15393.2283;
        -0.00168607392 0 0;
        -0.0719467433 0 0;
        -0.0673315968 0 0;
        -0.0221077377 0 0;
        -0.107638329 0 0;
        0.0416653549 0 0;
        0.0590660991 0 0;
        0.0861585559 0 0;
        0.0566417469 0 0;
        0.163113423 0 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1, 0.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1, 0.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0, 0.5),
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
    :Yb171,
    "P F=3/2, ν > 10", # [Phys. Rev. A 112, 042817 (2025)] # fit for ν > 10
    11,
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
        "6snf 3F2",
    ],
    Bool[1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1],
    [
        0.922094502 2.12370136 0;
        0.981191543 -4.54209175 0;
        0.229094016 0 0;
        0.206073107 0 0;
        0.193527627 0 0;
        0.181165673 0 0;
        0.925345494 -3.23594086 80.2535181;
        0.232649227 0 0;
        0.210070444 0 0;
        0.185699031 0 0;
        0.0718955585 -1.0913707 -38.4618954
    ],
    ["1.2", "1.3", "1.4", "1.5", "1.6", "2.3", "2.4", "2.5", "2.6", "7.8", "7.9", "7.10", "7.11"],
    [
        -0.102285383 153.251338 -15393.2283;
        -0.0719467433 0 0;
        -0.0673315968 0 0;
        -0.0221077377 0 0;
        -0.107638329 0 0;
        0.0416653549 0 0;
        0.0590660991 0 0;
        0.0861585559 0 0;
        0.0566417469 0 0;
        0.0703574701 0 0;
        0.0235308506 0 0;
        -0.0295876723 0 0;
        0.018377516 0 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 2, 1.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 0.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 1, 1.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 3, 2.5, 1.5),
    ]),
    [
        sqrt(5/3)/2 sqrt(5/6)/2 0 0 0 0 -sqrt(3/2)/2 0 0 0 0;
        -sqrt(1/3) sqrt(2/3) 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0;
        1/2 1/(2sqrt(2)) 0 0 0 0 sqrt(5/2)/2 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 0 -1
    ],
)

FMODEL_HIGHN_D05 = fModel(
    :Yb171,
    "D F=1/2, ν > 30", # [Phys. Rev. X 15, 011009 (2025)]
    1,
    ["6snd 3D1"],
    Bool[1],
    [0.75258093 0.382628525 -483.120633;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 1, 0.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 0.5)]),
    [-1;;],
)

FMODEL_HIGHN_D15 = fModel(
    :Yb171,
    "D F=3/2, ν > 30", # [Phys. Rev. A 112, 042817 (2025), Phys. Rev. X 15, 011009 (2025)]
    6,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2", "6snd 3D1"],
    Bool[1, 1, 0, 0, 1, 1],
    [
        0.73056016 -0.108286264 0;
        0.75155852 0.000367204397 0;
        0.195831577 0 0;
        0.236133225 0 0;
        0.147506921 0 0;
        0.75336354 -1.84349555 994.210321
    ],
    ["1.2", "1.3", "1.4", "2.4", "1.5", "2.5"],
    [
        0.22146327 -16.2798928;
        0.00431695191 0;
        0.0381576181 0;
        -0.00708200703 0;
        0.109346659 0;
        0.0636016813 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2, 1.5),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1, 1.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 1.5),
        fjQuantumNumbers(0.5, 1, NaN, NaN, 1, NaN, 1.5), # Jc of the 1D2 state could be both 1/2 or 3/2
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
    :Yb171,
    "D F=5/2, ν > 30", # [Phys. Rev. A 112, 042817 (2025), Phys. Rev. X 15, 011009 (2025)]
    6,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2", "6snd 3D3"],
    Bool[1, 1, 0, 0, 1, 1],
    [
        0.73056016 -0.108286264 0;
        0.75155852 0.000367204397 0;
        0.195831577 0 0;
        0.236133225 0 0;
        0.147506921 0 0;
        0.72861481 0.79979111 -484.236631
    ],
    ["1.2", "1.3", "1.4", "2.4", "1.5", "2.5"],
    [
        0.22146327 -16.2798928;
        0.00431695191 0;
        0.0381576181 0;
        -0.00708200703 0;
        0.109346659 0;
        0.0636016813 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2, 2.5),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2, 2.5),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2, 2.5),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3, 2.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 2.5),
        fjQuantumNumbers(0.5, 1, NaN, NaN, 1, NaN, 2.5), # Jc of the 1D2 state could be both 1/2 or 3/2
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
    :Yb171,
    "D F=7/2, ν > 14", # fit only valid from 30d upward [Phys. Rev. X 15, 011009 (2025)], provides good match around 18d
    1,
    ["6snd 3D3"],
    Bool[1],
    [0.72895315 -0.20653489 220.484722;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3, 3.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 3.5)]),
    [1;;],
)

FMODEL_HIGHN_F25 = fModel(
    :Yb171,
    "F F=5/2, ν > 20", # [Phys. Rev. A 112, 042817 (2025)]
    12,
    [
        "6snf 1F3",
        "6snf 3F3",
        "4f13 5d 6snl a",
        "4f13 5d 6snl b",
        "4f13 5d 6snl c",
        "4f13 5d 6snl d",
        "4f13 5d 6snl e",
        "6snf 3P2",
        "4f13 5d 6snl f",
        "4f13 5d 6snl g",
        "4f13 5d 6snl h",
        "6snf 3F2",
    ],
    Bool[1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
    [
        0.277086649 -13.2829133 0;
        0.0719837014 -0.741492064 0;
        0.251457795 0 0;
        0.227434828 0 0;
        0.175780645 0 0;
        0.196547521 0 0;
        0.21440857 0 0;
        0.925345494 -3.23594086 80.2535181;
        0.232649227 0 0;
        0.210070444 0 0;
        0.185699031 0 0;
        0.0718955585 -1.0913707 -38.4618954
    ],
    ["1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "2.3", "2.4", "2.5", "2.6", "2.7", "8.9", "8.10", "8.11", "8.12"],
    [
        -0.0209955122 0.251041249;
        -0.0585753224 0;
        -0.0750574327 0;
        0.122671919 0;
        -0.0401036164 0;
        0.0654271994 0;
        -0.0683007974 0;
        0.035415976 0;
        -0.0327625807 0;
        -0.050225071 0;
        0.0455759316 0;
        0.0703574701 0;
        0.0235308506 0;
        -0.0295876723 0;
        0.018377516 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 3, 3, 3, 2.5),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 3, 2.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2, 2.5),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 2, 2.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 3, 3.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 3, 2.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 3, 2.5, 2.5),
    ]),
    [
        -sqrt(4/7) -sqrt(3/7) 0 0 0 0 0 0 0 0 0 0;
        sqrt(5/7)/2 -sqrt(5/21) 0 0 0 0 0 0 0 0 0 sqrt(7/3)/2;
        0 0 1 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0;
        -1/2 1/sqrt(3) 0 0 0 0 0 0 0 0 0 sqrt(5/3)/2
    ],
)

FMODEL_HIGHN_F35 = fModel(
    :Yb171,
    "F F=7/2, ν > 20", # [Phys. Rev. A 112, 042817 (2025)]
    8,
    [
        "6snf 1F3",
        "6snf 3F3",
        "4f13 5d 6snl a",
        "4f13 5d 6snl b",
        "4f13 5d 6snl c",
        "4f13 5d 6snl d",
        "4f13 5d 6snl e",
        "6snf 3F4",
    ],
    Bool[1, 1, 0, 0, 0, 0, 0, 1],
    [
        0.277086649 -13.290196 0;
        0.0719837014 -0.754736076 0;
        0.251457795 0 0;
        0.227434828 0 0;
        0.175780645 0 0;
        0.196547521 0 0;
        0.21440857 0 0;
        0.0834193873 -1.11453386 -1545.71844
    ],
    ["1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "2.3", "2.4", "2.5", "2.6", "2.7"],
    [
        0.0209955122 0.251041249;
        -0.0585753224 0;
        -0.0750574327 0;
        0.122671919 0;
        -0.0401036164 0;
        0.0654271994 0;
        -0.0683007974 0;
        0.035415976 0;
        -0.0327625807 0;
        -0.050225071 0;
        0.0455759316 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 3, 3, 3, 3.5),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 3, 3.5),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 4, 3.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 3, 3.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 3, 2.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 3, 3.5, 3.5),
    ]),
    [
        3/(2sqrt(7)) 3sqrt(3)/(4sqrt(7)) 0 0 0 0 0 -sqrt(7)/4;
        -sqrt(3/7) sqrt(4/7) 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;
        1/2 sqrt(3)/4 0 0 0 0 0 3/4
    ],
)

FMODEL_HIGHN_F45 = fModel(
    :Yb171,
    "F F=9/2, ν > 20", # extracted from FMODEL_HIGHN_F35
    1,
    ["6snf 3F4"],
    Bool[1],
    [0.0834193873 -1.11453386 -1545.71844;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 3, 3, 4, 4.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 3, 3.5, 4.5)]),
    [1;;],
)

FMODEL_HIGHN_G25 = fModel(
    :Yb171,
    "G F=5/2, ν > 25", # [Phys. Rev. A 112, 042817 (2025)]
    1,
    ["6sng 3G3"],
    Bool[1],
    [0.02613255 -0.14203905;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 4, 4, 3, 2.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 4, 3.5, 2.5)]),
    [-1;;],
)

FMODEL_HIGHN_G35 = fModel(
    :Yb171,
    "G F=7/2, ν > 25", # [Phys. Rev. A 112, 042817 (2025)]
    3,
    ["6sng +G4", "6sng -G4", "6sng 3G3"],
    Bool[1, 1, 1],
    [
        0.02628545 -0.13182564;
        0.02548145 -0.12028462;
        0.02613255 -0.14203905
    ],
    ["1.2"],
    [-0.089123698 0;],
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 4, 4.5, 4, 3.5),
        jjQuantumNumbers(0.5, 0, 0.5, 4, 3.5, 4, 3.5),
        jjQuantumNumbers(0.5, 0, 0.5, 4, 3.5, 3, 3.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 4, 4.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 4, 3.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 4, 3.5, 3.5),
    ]),
    [
        -1 0 0;
        0 -sqrt(7)/4 3/4;
        0 3/4 sqrt(7)/4
    ],
)

FMODEL_HIGHN_G45 = fModel(
    :Yb171,
    "G F=9/2, ν > 25", # [Phys. Rev. A 112, 042817 (2025)]
    3,
    ["6sng +G4", "6sng -G4", "6sng 3G5"],
    Bool[1, 1, 1],
    [
        0.02628545 -0.13182564;
        0.02548145 -0.12028462;
        0.02536571 -0.18507079
    ],
    ["1.2"],
    [-0.089123698 0;],
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 4, 4.5, 4, 4.5),
        jjQuantumNumbers(0.5, 0, 0.5, 4, 3.5, 4, 4.5),
        jjQuantumNumbers(0.5, 0, 0.5, 4, 4.5, 5, 4.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 4, 4.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 4, 3.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 4, 4.5, 4.5),
    ]),
    [
        sqrt(11/20) 0 -sqrt(9/20);
        0 1 0;
        sqrt(9/20) 0 sqrt(11/20)
    ],
)

FMODEL_HIGHN_G55 = fModel(
    :Yb171,
    "G F=11/2, ν > 25", # extracted from FMODEL_HIGHN_G45
    1,
    ["6sng 3G5"],
    Bool[1],
    [0.02536571 -0.18507079;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 4, 4, 5, 5.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 4, 4.5, 5.5)]),
    [1;;],
)

# --------------------------------------------------------
# MQDT models valid at small n
# --------------------------------------------------------

FMODEL_LOWN_S05 = fModel(
    :Yb171,
    "S F=1/2, 2 < ν < 26", # model from [Phys. Rev. X 15, 011009 (2025)] with modified 3S1 data taken from Yb174 fit to NIST
    7,
    ["6sns 1S0", "4f13 5d 6snl a", "6pnp 1S0", "4f13 5d 6snl b", "6pnp 3P0", "4f13 5d 6snl c", "6sns 3S1"],
    Bool[1, 0, 1, 0, 1, 0, 1],
    [
        0.357488757 0.163255076 0;
        0.203917828 0 0;
        0.116813499 0 0;
        0.287210377 0 0;
        0.247550262 0 0;
        0.148686263 0 0;
        0.432841 0.724559 -1.95424
    ],
    ["1.2", "1.3", "1.4", "3.4", "3.5", "1.6"],
    [
        0.13179534 0;
        0.29748039 0;
        0.0553920359 0;
        0.100843905 0;
        0.10317753 0;
        0.137709223 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 0, 0, 0, 0.5),
        lsQuantumNumbers(0.5, 0, 1, 1, 0, 0, 0.5),
        lsQuantumNumbers(0.5, 1, 1, 1, 1, 0, 0.5),
        lsQuantumNumbers(0.5, 1, 0, 0, 0, 1, 0.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 0, 0, 0.5, 0.5),
        fjQuantumNumbers(0.5, 1, 1.5, 1, 1, 1.5, 0.5),
        fjQuantumNumbers(0.5, 1, 0.5, NaN, 1, 0.5, 0.5), # Fc of the 3P0 state could be both 0 or 1
        fjQuantumNumbers(0.5, 0, 0.5, 1, 0, 0.5, 0.5),
    ]),
    [
        1/2 0 0 0 0 0 sqrt(3)/2;
        0 1 0 0 0 0 0;
        0 0 sqrt(2/3) 0 -sqrt(1/3) 0 0;
        0 0 0 1 0 0 0;
        0 0 sqrt(1/3) 0 sqrt(2/3) 0 0;
        0 0 0 0 0 1 0;
        sqrt(3)/2 0 0 0 0 0 -1/2
    ], # corrected in [arXiv:2507.11487v1]
)

FMODEL_LOWN_S15 = fModel(
    :Yb171,
    "S F=3/2, 2 < ν < 26", # taken from Yb174 fit to NIST
    1,
    ["6sns 3S1"],
    Bool[1],
    [0.432841 0.724559 -1.95424;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1, 1.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 0, 0.5, 1.5)]),
    [1;;],
)

FMODEL_LOWEST_P05 = fModel(
    :Yb171,
    "P F=1/2, 1.5 < ν < 2.5", # fit to Yb174 NIST data
    3,
    ["6snp 1P1", "6snp 3P1", "6snp 3P0"],
    Bool[1, 1, 1],
    [0.161083 0; 0.920424 0; 0.180701 0],
    ["1.2"],
    [-0.426128 6.272986;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1, 0.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1, 0.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0, 0.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 0.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 0.5, 0.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 1, 0.5, 0.5),
    ]),
    [
        -sqrt(2/3) -sqrt(1/3) 0;
        1/(2sqrt(3)) -sqrt(1/6) sqrt(3)/2;
        -1/2 1/sqrt(2) 1/2
    ],
)

FMODEL_LOWN_P05 = fModel(
    :Yb171,
    "P F=1/2, 2.9 < ν < 5.9", # data taken from Yb174 NIST fits
    3,
    ["6snp 1P1", "6snp 3P1", "6snp 3P0"],
    Bool[1, 1, 1],
    [
        0.967223 -3.03997 0.569205;
        0.967918 0.25116 0.868505;
        0.969279 0.288219 1.36228
    ],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1, 0.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1, 0.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0, 0.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 0.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 0.5, 0.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 1, 0.5, 0.5),
    ]),
    [
        -sqrt(2/3) -sqrt(1/3) 0;
        1/(2sqrt(3)) -sqrt(1/6) sqrt(3)/2;
        -1/2 1/sqrt(2) 1/2
    ],
)

FMODEL_LOWEST_P15 = fModel(
    :Yb171,
    "P F=3/2, 1.5 < ν < 2.5", # fit to Yb174 NIST data
    3,
    ["6snp 1P1", "6snp 3P1", "6snp 3P2"],
    Bool[1, 1, 1],
    [0.161083 0; 0.920424 0; 0.110501 0],
    ["1.2"],
    [-0.426128 6.272986;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2, 1.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 0.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 0, 1, 1.5, 1.5),
    ]),
    [
        sqrt(5/3)/2 sqrt(5/6)/2 -sqrt(3/2)/2;
        -sqrt(1/3) sqrt(2/3) 0;
        1/2 1/(2sqrt(2)) sqrt(5/2)/2
    ],
)

FMODEL_LOWN_P15 = fModel(
    :Yb171,
    "P F=3/2, 3 < ν < 10", # data taken from Yb174 NIST fits
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
        0.967223 -3.03997 0.569205;
        0.967918 0.25116 0.868505;
        0.228828720 0 0;
        0.205484818 0 0;
        0.193528629 0 0;
        0.181385000 0 0;
        0.906105 0.383471 1.23512;
        0.236866903 0 0;
        0.221055883 0 0;
        0.185599376 0 0
    ],
    ["1.2", "1.3", "1.4", "1.5", "1.6", "2.3", "2.4", "2.5", "2.6", "7.8", "7.9", "7.10"],
    [
        -0.08410871 120.37555 -9314.23;
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
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2, 1.5),
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

FMODEL_LOWEST_P25 = fModel(
    :Yb171,
    "P F=5/2, 1.5 < ν < 4.5", # data taken from Yb174 NIST fit
    1,
    ["6snp 3P2"],
    Bool[1],
    [0.906105 0.383471 1.23512;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2, 2.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 2.5)]),
    [1;;],
)

FMODEL_LOWN_P25 = fModel(
    :Yb171,
    "P F=5/2, 5 < ν < 20", # data taken from Yb174 [Phys. Rev. X 15, 011009 (2025)] fits
    4,
    ["6snp 3P2", "4f13 5d 6snl a", "4f13 5d 6snl b", "4f13 5d 6snl c"],
    Bool[1, 0, 0, 0],
    [
        0.925121305 -2.73247165 74.664989;
        0.230133261 0 0;
        0.209638118 0 0;
        0.186228192 0 0
    ],
    ["1.2", "1.3", "1.4"],
    [
        0.0706666127 0;
        0.0232711158 0;
        -0.0292153659 0
    ],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2, 2.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 1, 1.5, 2.5)]),
    [
        1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1
    ],
)

FMODEL_LOWN_D05 = fModel(
    :Yb171,
    "D F=1/2, 2 < ν < 30", # data taken from Yb174 NIST fit for 2 < ν < 5
    1,
    ["6snd 3D1"],
    Bool[1],
    [0.758222 -0.017906 3.392161;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 1, 0.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 0.5)]),
    [-1;;],
)

FMODEL_LOWN_D15 = fModel(
    :Yb171,
    "D F=3/2, 2 < ν < 30", # model from [arXiv:2507.11487v1] with modified 3D1 data taken from Yb174 fit to NIST
    6,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2", "6snd 3D1"],
    Bool[1, 1, 0, 0, 1, 1],
    [
        0.730541589 -0.0967938662 0;
        0.751542685 0.00038836127 0;
        0.195864083 0 0;
        0.235944408 0 0;
        0.147483609 0 0;
        0.758222 -0.017906 3.392161
    ],
    ["1.2", "1.3", "1.4", "2.4", "1.5", "2.5"],
    [
        0.220048245 0;
        0.00427599 0;
        0.0381563093 0;
        -0.00700797918 0;
        0.109380331 0;
        0.0635544456 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2, 1.5),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2, 1.5),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1, 1.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 1.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 1.5),
        fjQuantumNumbers(0.5, 1, NaN, NaN, 1, NaN, 1.5), # Jc of the 1D2 state could be both 1/2 or 3/2
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

FMODEL_LOWN_D25 = fModel(
    :Yb171,
    "D F=5/2, 2 < ν < 30", # model from [arXiv:2507.11487v1] with modified 3D3 data taken from Yb174 fit to NIST
    6,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2", "6snd 3D3"],
    Bool[1, 1, 0, 0, 1, 1],
    [
        0.730541589 -0.0967938662 0;
        0.751542685 0.00038836127 0;
        0.195864083 0 0;
        0.235944408 0 0;
        0.147483609 0 0;
        0.734512 -0.019501 3.459114
    ],
    ["1.2", "1.3", "1.4", "2.4", "1.5", "2.5"],
    [
        0.220048245 -14.9486;
        0.00427599 0;
        0.0381563093 0;
        -0.00700797918 0;
        0.109380331 0;
        0.0635544456 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2, 2.5),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2, 2.5),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2, 2.5),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3, 2.5),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 2.5),
        fjQuantumNumbers(0.5, 1, NaN, NaN, 1, NaN, 2.5), # Jc of the 1D2 state could be both 1/2 or 3/2
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

FMODEL_LOWN_D35 = fModel(
    :Yb171,
    "D F=7/2, 2 < ν < 14", # data taken from Yb174 NIST fit for 2 < ν < 5, provides good match around 18d
    1,
    ["6snd 3D3"],
    Bool[1],
    [0.734512 -0.019501 3.459114;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3, 3.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 3.5)]),
    [1;;],
)

end
