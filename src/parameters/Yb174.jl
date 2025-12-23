module Yb174

using ..MQDT: Parameters, fModel, lsChannels, jjChannels, lsQuantumNumbers, jjQuantumNumbers, coreQuantumNumbers

export PARA,
    FMODEL_HIGHN_S0,
    FMODEL_HIGHN_S1,
    FMODEL_HIGHN_P0,
    FMODEL_HIGHN_P1,
    FMODEL_HIGHN_P2,
    FMODEL_HIGHN_D1,
    FMODEL_HIGHN_D2,
    FMODEL_HIGHN_D3,
    FMODEL_HIGHN_F2,
    FMODEL_HIGHN_F3,
    FMODEL_HIGHN_F4,
    FMODEL_HIGHN_G3,
    FMODEL_HIGHN_G4,
    FMODEL_HIGHN_G5,
    FMODEL_LOWN_S0,
    FMODEL_LOWN_S1,
    FMODEL_LOWN_P0,
    FMODEL_LOWEST_P1,
    FMODEL_LOWN_P1,
    FMODEL_LOWN_P2,
    FMODEL_LOWN_D1,
    FMODEL_LOWN_D2,
    FMODEL_LOWN_D3

# Isotope data
THRESHOLDS = Dict(
    coreQuantumNumbers(0, 0.5) => 50443.070393,
    coreQuantumNumbers(1, 0.5) => 77504.98,
    coreQuantumNumbers(1, 1.5) => 80835.39,
    "4f13 5d 6s" => 83967.7,
)

PARA = Parameters(
    :Yb174,
    1822.88848628*173.9388621, # nuclear mass
    0, # nuclear spin
    109736.9695858, # Rydberg constant in 1/cm
    THRESHOLDS[coreQuantumNumbers(0, 0.5)], # lowest ionization threshold in 1/cm
    0, # hyperfine constant in 1/cm
    2.1, # nuclear dipole
    THRESHOLDS,
)

# --------------------------------------------------------
# MQDT models valid at large n
# --------------------------------------------------------

FMODEL_HIGHN_S0 = fModel(
    :Yb174,
    "S J=0, ν > 2", # fit for states 6s7s upward [Phys. Rev. X 15, 011009 (2025)]
    6,
    ["6sns 1S0", "4f13 5d 6snl a", "6pnp 1S0", "4f13 5d 6snl b", "6pnp 3P0", "4f13 5d 6snl c"],
    Bool[1, 0, 1, 0, 1, 0],
    [
        0.355101645 0.277673956;
        0.204537535 0;
        0.116393648 0;
        0.295439966 0;
        0.257664798 0;
        0.155797119 0
    ],
    ["1.2", "1.3", "1.4", "3.4", "3.5", "1.6"],
    [
        0.126557575 0;
        0.300103593 0;
        0.056987912 0;
        0.114312578 0;
        0.0986363362 0;
        0.142498543 0
    ],
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
        0 0 sqrt(2/3) 0 -sqrt(1/3) 0;
        0 0 0 1 0 0;
        0 0 sqrt(1/3) 0 sqrt(2/3) 0;
        0 0 0 0 0 1
    ], # updated in [arXiv:2507.11487v1]
)

FMODEL_HIGHN_S1 = fModel(
    :Yb174,
    "S J=1, ν > 26", # fit only valid from 28s upward [Phys. Rev. Lett. 128, 033201 (2022)]
    1,
    ["6sns 3S1"],
    Bool[1],
    [0.4382 4 -1e4 8e6 -3e9;], # adjusted for non-rydberg ritz
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 1)]),
    [1;;],
)

FMODEL_HIGHN_P0 = fModel(
    :Yb174,
    "P J=0, ν > 6", # [Phys. Rev. X 15, 011009 (2025)]
    2,
    ["6snp 3P0", "4f13 5d 6snd"],
    Bool[1, 0],
    [0.953661478 -0.287531374; 0.198460766 0],
    ["1.2"],
    [0.163343232 0],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 0)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 0)]),
    [1 0; 0 1],
)

FMODEL_HIGHN_P1 = fModel(
    :Yb174,
    "P J=1, ν > 6", # [Phys. Rev. X 15, 011009 (2025)]
    6,
    ["6snp 1P1", "6snp 3P1", "4f13 5d 6snl a", "4f13 5d 6snl b", "4f13 5d 6snl c", "4f13 5d 6snl d"],
    Bool[1, 1, 0, 0, 0, 0],
    [
        0.92271098 2.6036257;
        0.98208719 -5.4562725;
        0.22851720 0;
        0.20607759 0;
        0.19352751 0;
        0.18153094 0
    ],
    ["1.2", "1.3", "1.4", "1.5", "1.6", "2.3", "2.4", "2.5", "2.6"],
    [
        -0.08410871 120.37555 -9314.23;
        -0.07317986 0 0;
        -0.06651879 0 0;
        -0.02212194 0 0;
        -0.10452109 0 0;
        0.02477464 0 0;
        0.05763934 0 0;
        0.0860644 0 0;
        0.04993818 0 0
    ],
    lsChannels([lsQuantumNumbers(0.5, 0, 0, 1, 1, 1), lsQuantumNumbers(0.5, 1, 0, 1, 1, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 1), jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 1)]),
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
    :Yb174,
    "P J=2, ν > 5", # [Phys. Rev. X 15, 011009 (2025)]
    4,
    ["6snp 3P2", "4f13 5d 6snl a", "4f13 5d 6snl b", "4f13 5d 6snl c"],
    Bool[1, 0, 0, 0],
    [
        0.925150932 -2.69197178 66.7159709;
        0.230028034 0 0;
        0.209224174 0 0;
        0.186236574 0 0
    ],
    ["1.2", "1.3", "1.4"],
    [
        0.0706189664 0;
        0.0231221428 0;
        -0.0291730345 0
    ],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 2)]),
    [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],
)

FMODEL_HIGHN_D1 = fModel(
    :Yb174,
    "D J=1, ν > 26", # fit only valid from 30d upward [Phys. Rev. X 15, 011009 (2025)]
    1,
    ["6snd 3D1"],
    Bool[1],
    [0.75258093 0.3826 -483.1;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 1)]),
    [1;;],
)

FMODEL_HIGHN_D2 = fModel(
    :Yb174,
    "D J=2, ν > 5", # [Phys. Rev. X 15, 011009 (2025)]
    5,
    ["6snd 1D2", "6snd 3D2", "4f13 5d 6snl a", "4f13 5d 6snl b", "6pnp 1D2"],
    Bool[1, 1, 0, 0, 1],
    [
        0.729513646 -0.0377841183;
        0.752292223 0.104072325;
        0.19612036 0;
        0.233752026 0;
        0.152911249 0
    ],
    ["1.2", "1.3", "1.4", "2.4", "1.5", "2.5"],
    [
        0.21157531 -15.3844;
        0.00521559431 0;
        0.0398131577 0;
        -0.0071658109 0;
        0.10481227 0;
        0.0721660042 0
    ],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 2),
        jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 2),
        jjQuantumNumbers(0.5, 1, NaN, 1, NaN, 2), # Jc and Jr could also be switched or both be 3/2
    ]),
    [
        sqrt(3/5) sqrt(2/5) 0 0 0;
        -sqrt(2/5) sqrt(3/5) 0 0 0;
        0 0 1 0 0;
        0 0 0 1 0;
        0 0 0 0 1
    ],
    Dict(coreQuantumNumbers(1, NaN) => 79725.35),
)

FMODEL_HIGHN_D3 = fModel(
    :Yb174,
    "D J=3, ν > 18", # fit only valid from 30d upward [Phys. Rev. X 15, 011009 (2025)], provides good match around 21d
    1,
    ["6snd 3D3"],
    Bool[1],
    [0.72902016 -0.705328923 829.238844;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 3)]),
    [1;;],
)

FMODEL_HIGHN_F2 = fModel(
    :Yb174,
    "F J=2, ν > 25", # [arXiv:2507.11487v1]
    1,
    ["6snf 3F2"],
    Bool[1],
    [0.0718252326 -1.00091963 -106.291066;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 3, 3, 2)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 3, 2.5, 2)]),
    [1;;],
)

FMODEL_HIGHN_F3 = fModel(
    :Yb174,
    "F J=3, ν > 7", # [Phys. Rev. A 112, 042817 (2025)]
    7,
    ["6snf 1F3", "6snf 3F3", "4f13 5d 6snl a", "4f13 5d 6snl b", "4f13 5d 6snl c", "4f13 5d 6snl d", "4f13 5d 6snl e"],
    Bool[1, 1, 0, 0, 0, 0, 0],
    [
        0.276158949 -12.7258012;
        0.0715123712 -0.768462937;
        0.239015576 0;
        0.226770354 0;
        0.175354845 0;
        0.196660618 0;
        0.21069642 0
    ],
    ["1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "2.3", "2.4", "2.5", "2.6", "2.7"],
    [
        -0.0209955122 0.251041249; # updated since [arXiv:2507.11487v1]
        -0.00411835457 0;
        -0.0962784945 0;
        0.132826901 0;
        -0.0439244317 0;
        0.0508460294 0;
        -0.0376574252 0;
        0.026944623 0;
        -0.0148474857 0;
        -0.0521244126 0;
        0.0349516329 0
    ],
    lsChannels([lsQuantumNumbers(0.5, 0, 0, 3, 3, 3), lsQuantumNumbers(0.5, 1, 0, 3, 3, 3)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 3, 3.5, 3), jjQuantumNumbers(0.5, 0, 0.5, 3, 2.5, 3)]),
    [
        sqrt(4/7) sqrt(3/7) 0 0 0 0 0;
        -sqrt(3/7) sqrt(4/7) 0 0 0 0 0;
        0 0 1 0 0 0 0;
        0 0 0 1 0 0 0;
        0 0 0 0 1 0 0;
        0 0 0 0 0 1 0;
        0 0 0 0 0 0 1
    ],
)

FMODEL_HIGHN_F4 = fModel(
    :Yb174,
    "F J=4, ν > 25", # [arXiv:2507.11487v1]
    1,
    ["6snf 3F4"],
    Bool[1],
    [0.0839027969 -2.91009023;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 3, 3, 4)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 3, 3.5, 4)]),
    [1;;],
)

FMODEL_HIGHN_G3 = fModel(
    :Yb174,
    "G J=3, ν > 25", # [arXiv:2507.11487v1]
    1,
    ["6sng 3G3"],
    Bool[1],
    [0.0260964574 -0.14139526;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 4, 4, 3)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 4, 3.5, 3)]),
    [1;;],
)

FMODEL_HIGHN_G4 = fModel(
    :Yb174,
    "G J=4, ν > 25", # [Phys. Rev. A 112, 042817 (2025)]
    2,
    ["6sng +G4", "6sng -G4"],
    Bool[1, 1],
    [0.0262659964 -0.148808463; 0.0254568575 -0.134219071;],
    ["1.2"],
    [-0.089123698 0;], # updated since [arXiv:2507.11487v1]
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 4, 4.5, 4), jjQuantumNumbers(0.5, 0, 0.5, 4, 3.5, 4)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 4, 4.5, 4), jjQuantumNumbers(0.5, 0, 0.5, 4, 3.5, 4)]),
    [1 0; 0 1], # this series is well described in jj coupling. singlet-triplet mixing is approximately atan(sqrt(4/5))
)

FMODEL_HIGHN_G5 = fModel(
    :Yb174,
    "G J=5, ν > 25", # [Phys. Rev. A 112, 042817 (2025)]
    1,
    ["6sng 3G5"],
    Bool[1],
    [0.02536571 -0.18507079;], # updated since [arXiv:2507.11487v1]
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 4, 4, 5)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 4, 4.5, 5)]),
    [1;;],
)

# --------------------------------------------------------
# MQDT models valid at small n
# --------------------------------------------------------

FMODEL_LOWN_S0 = fModel(
    :Yb174,
    "S J=0, 1 < ν < 2", # fit to the 6s^2 ground state
    1,
    ["6sns 1S0"],
    Bool[1],
    [0.525055 0;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 0, 0, 0, 0, 0)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 0)]),
    [1;;],
)

FMODEL_LOWN_S1 = fModel(
    :Yb174,
    "S J=1, 2 < ν < 26", # fit to NIST data between 7s and 13s, extrapolation seems good up to 30s
    1,
    ["6sns 3S1"],
    Bool[1],
    [0.432841 0.724559 -1.95424;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 1)]),
    [1;;],
)

FMODEL_LOWN_P0 = fModel(
    :Yb174,
    "P J=0, 1.5 < ν < 5.5", # fit to NIST data between 6p and 9p
    1,
    ["6snp 3P0"],
    Bool[1],
    [0.969279 0.288219 1.36228;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 0)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 0)]),
    [1;;],
)

FMODEL_LOWEST_P1 = fModel(
    :Yb174,
    "P J=1, 1.7 < ν < 2.7", # fit to NIST data for the 6p states
    2,
    ["6snp 1P1", "6snp 3P1"],
    Bool[1, 1],
    [0.161083 0; 0.920424 0],
    ["1.2"],
    [-0.426128 6.272986;],
    lsChannels([lsQuantumNumbers(0.5, 0, 0, 1, 1, 1), lsQuantumNumbers(0.5, 1, 0, 1, 1, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 1), jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 1)]),
    [sqrt(2/3) sqrt(1/3); -sqrt(1/3) sqrt(2/3)],
)

FMODEL_LOWN_P1 = fModel(
    :Yb174,
    "P J=1, 2.7 < ν < 5.7", # fit to NIST data between 7p and 9p
    2,
    ["6snp 1P1", "6snp 3P1"],
    Bool[1, 1],
    [0.967223 -3.03997 0.569205; 0.967918 0.25116 0.868505;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 0, 0, 1, 1, 1), lsQuantumNumbers(0.5, 1, 0, 1, 1, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 1), jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 1)]),
    [sqrt(2/3) sqrt(1/3); -sqrt(1/3) sqrt(2/3)],
)

FMODEL_LOWN_P2 = fModel(
    :Yb174,
    "P J=2, 1.5 < ν < 4.5", # fit to NIST data between 6p and 8p
    1,
    ["6snp 3P2"],
    Bool[1],
    [0.906105 0.383471 1.23512;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 2)]),
    [1;;],
)

FMODEL_LOWN_D1 = fModel(
    :Yb174,
    "D J=1, 2 < ν < 5", # fit to NIST data between 5d and 8d, causes a µ=0.005 difference at the 30d state
    1,
    ["6snd 3D1"],
    Bool[1],
    [0.758222 -0.017906 3.392161;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 1)]),
    [1;;],
)

FMODEL_LOWN_D2 = fModel(
    :Yb174,
    "D J=2, 2 < ν < 5", # fit to NIST data between 5d and 7d
    2,
    ["6snd 1D2", "6snd 3D2"],
    Bool[1, 1],
    [0.703156 0.973192; 0.724546 0.372621],
    ["1.2"],
    [0.948409 2.121270],
    lsChannels([lsQuantumNumbers(0.5, 0, 0, 2, 2, 2), lsQuantumNumbers(0.5, 1, 0, 2, 2, 2)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 2), jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 2)]),
    [sqrt(3/5) sqrt(2/5); -sqrt(2/5) sqrt(3/5)],
)

FMODEL_LOWN_D3 = fModel(
    :Yb174,
    "D J=3, 2 < ν < 18", # fit to NIST data between 5d and 8d, provides good match around 21d
    1,
    ["6snd 3D3"],
    Bool[1],
    [0.734512 -0.019501 3.459114;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 3)]),
    [1;;],
)

end
