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
    FMODEL_HIGHN_G55

# Isotope data
THRESHOLDS = Dict(
    coreQuantumNumbers(0, 0.5, 0) => 50442.795744,
    coreQuantumNumbers(0, 0.5, 1) => 50443.217463,
    coreQuantumNumbers(1, 0.5, NaN) => 77504.98,
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

# MQDT Models
FMODEL_HIGHN_S05 = fModel(
    :Yb171,
    "S F=1/2, ν > 26", # [Phys. Rev. X 15, 011009 (2025)]
    7,
    ["6sns 1S0", "4f13 5d 6snl a", "6pnp 1S0", "4f13 5d 6snl b", "6pnp 3P0", "4f13 5d 6snl c", "6sns 3S1"],
    Bool[1, 0, 1, 0, 1, 0, 1],
    [
        0.357488757 0.163255076 0 0 0;
        0.203917828 0 0 0 0;
        0.116813499 0 0 0 0;
        0.287210377 0 0 0 0;
        0.247550262 0 0 0 0;
        0.148686263 0 0 0 0;
        0.438803844 3.54561559 -10673.2496 7702824.55 -2255973430
    ],
    ["12", "13", "14", "34", "35", "16"],
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
        fjQuantumNumbers(0.5, 1, 0.5, 0, 1, 0.5, 0.5), # Fc of the 3P0 state could be both 0 or 1 TODO wouldnt this change the ionization threshold?
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
    [0.438426851 3.91762642 -10612.6828 8017432.38 -2582622910.0],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1, 1.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 0, 0.5, 1.5)]),
    [1;;],
)

FMODEL_HIGHN_P05 = fModel(
    :Yb171,
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
    [
        0.922046003 2.17932092;
        0.981076381 -4.48995034;
        0.229019212 0;
        0.205999445 0;
        0.19352756 0;
        0.181271472 0;
        0.953288122 -0.0298396847;
        0.198451139 0
    ],
    ["12", "27", "13", "14", "15", "16", "23", "24", "25", "26", "78"],
    [
        -0.100292816 149.140925 -13487.7692;
        -0.0016592076 0 0;
        -0.0727917308 0 0;
        -0.0669120237 0 0;
        -0.0221321759 0 0;
        -0.107302569 0 0;
        0.0396527798 0 0;
        0.0596597186 0 0;
        0.0861416497 0 0;
        0.0565415641 0 0;
        0.163175562 0 0
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
    "P F=3/2, ν > 10", # [arXiv:2507.11487v1] # fit for ν > 10
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
        "62nf 3F2",
    ],
    Bool[1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1],
    [
        0.922046003 2.17932092 0;
        0.981076381 -4.48995034 0;
        0.229019212 0 0;
        0.205999445 0 0;
        0.19352756 0 0;
        0.181271472 0 0;
        0.926908753 -4.6228148 453.565911;
        0.223859222 0 0;
        0.240057861 0 0;
        0.184824687 0 0;
        0.0718810048 -1.08216233 -38.2507093
    ],
    ["12", "13", "14", "15", "16", "23", "24", "25", "26", "78", "79", "710", "117"],
    [
        -0.100292816 149.140925 -13487.7692;
        -0.0727917308 0 0;
        -0.0669120237 0 0;
        -0.0221321759 0 0;
        -0.107302569 0 0;
        0.0396527798 0 0;
        0.0596597186 0 0;
        0.0861416497 0 0;
        0.0565415641 0 0;
        0.0657620037 0 0;
        0.0215995148 0 0;
        -0.028252844 0 0;
        -0.0184267375 0 0; # the minus sign is intentional, to account for "backward" coupling, because "117" couples 11 with 7, coupling to channels 1 and 2 is neglected so far
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
    [0.75258093 0.382628525 -483.120633],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 1, 0.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 1.5, 0.5)]),
    [-1;;],
)

FMODEL_HIGHN_D15 = fModel(
    :Yb171,
    "D F=3/2, ν > 30", # [arXiv:2507.11487v1]
    6,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2", "6snd 3D1"],
    Bool[1, 1, 0, 0, 1, 1],
    [
        0.730541589 -0.0967938662 0;
        0.751542685 0.00038836127 0;
        0.195864083 0 0;
        0.235944408 0 0;
        0.147483609 0 0;
        0.75334159 -1.80187083 986.918851
    ],
    ["12", "13", "14", "24", "15", "25"],
    [
        0.220048245 -14.9486;
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
    Dict(coreQuantumNumbers(1, 1.5) => 79725.35),
)

FMODEL_HIGHN_D25 = fModel(
    :Yb171,
    "D F=5/2, ν > 30", # [arXiv:2507.11487v1]
    6,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2", "6snd 3D3"],
    Bool[1, 1, 0, 0, 1, 1],
    [
        0.730541589 -0.0967938662 0;
        0.751542685 0.00038836127 0;
        0.195864083 0 0;
        0.235944408 0 0;
        0.147483609 0 0;
        0.72865616 0.793185994 -523.858959
    ],
    ["12", "13", "14", "24", "15", "25"],
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
        fjQuantumNumbers(0.5, 1, 1.5, 1, 1, 1.5, 2.5), # Fc of the 1D2 state could be both 0 or 1
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
    Dict(coreQuantumNumbers(1, 1.5) => 79725.35),
)

FMODEL_HIGHN_D35 = fModel(
    :Yb171,
    "D F=7/2, ν > 14", # fit only valid from 30d upward [Phys. Rev. X 15, 011009 (2025)], provides good match around 18d
    1,
    ["6snd 3D3"],
    Bool[1],
    [0.72895315 -0.20653489 220.484722],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3, 3.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 2, 2.5, 3.5)]),
    [1;;],
)

FMODEL_HIGHN_F25 = fModel(
    :Yb171,
    "F F=5/2, ν > 20", # [arXiv:2507.11487v1]
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
        0.27707535 -13.2829133 0;
        0.0719827367 -0.741492064 0;
        0.250829254 0 0;
        0.227042016 0 0;
        0.17576256 0 0;
        0.19654694 0 0;
        0.214351625 0 0;
        0.926908753 -4.6228148 453.565911;
        0.223859222 0 0;
        0.240057861 0 0;
        0.184824687 0 0;
        0.0718810048 -1.08216233 -38.2507093
    ],
    ["12", "13", "14", "15", "16", "17", "23", "24", "25", "26", "27", "89", "810", "118", "128"],
    [
        -0.0208481417 0.239045493;
        -0.0593782568 0;
        -0.0755947274 0;
        0.122678758 0;
        -0.0400418576 0;
        0.0646373252 0;
        -0.0684538786 0;
        0.0352130279 0;
        -0.0326572035 0;
        -0.050215299 0;
        0.0453892695 0;
        0.0657620037 0;
        0.0215995148 0;
        0.028252844 0; # note the intentional minus sign to account for inverse coupling
        -0.0184267375 0 # note the intentional minus sign to account for inverse coupling
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
    "F F=7/2, ν > 20", # [arXiv:2507.11487v1]
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
        0.27707535 -13.2829133 0;
        0.0719827367 -0.741492064 0;
        0.250829254 0 0;
        0.227042016 0 0;
        0.17576256 0 0;
        0.19654694 0 0;
        0.214351625 0 0;
        0.0834066138 -1.11349283 -1539.63739
    ],
    ["12", "13", "14", "15", "16", "17", "23", "24", "25", "26", "27"],
    [
        -0.0208481417 0.239045493;
        -0.0593782568 0;
        -0.0755947274 0;
        0.122678758 0;
        -0.0400418576 0;
        0.0646373252 0;
        -0.0684538786 0;
        0.0352130279 0;
        -0.0326572035 0;
        -0.050215299 0;
        0.0453892695 0
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
    [0.0834066138 -1.11349283 -1539.63739;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 3, 3, 4, 4.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 3, 3.5, 4.5)]),
    [1;;],
)

FMODEL_HIGHN_G25 = fModel(
    :Yb171,
    "G F=5/2, ν > 25", # [arXiv:2507.11487v1]
    1,
    ["6sng 3G3"],
    Bool[1],
    [0.02611627 -0.13079392],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 4, 4, 3, 2.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 4, 3.5, 2.5)]),
    [-1;;],
)

FMODEL_HIGHN_G35 = fModel(
    :Yb171,
    "G F=7/2, ν > 25", # [arXiv:2507.11487v1]
    3,
    ["6sng +G4", "6sng -G4", "6sng 3G3"],
    Bool[1, 1, 1],
    [
        0.02627822 -0.13288481;
        0.02536987 -0.11668088;
        0.02611627 -0.13079392
    ],
    ["12"],
    [-0.08222676 0;],
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
    "G F=9/2, ν > 25", # [arXiv:2507.11487v1]
    3,
    ["6sng +G4", "6sng -G4", "6sng 3G5"],
    Bool[1, 1, 1],
    [
        0.02627822 -0.13288481;
        0.02536987 -0.11668088;
        0.02547095 -0.2007948
    ],
    ["12"],
    [-0.08222676 0;],
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
    [0.02547095 -0.2007948],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 4, 4, 5, 5.5)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 1, 4, 4.5, 5.5)]),
    [1;;],
)

end
