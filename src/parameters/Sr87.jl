module Sr87

using ..MQDT: Parameters, fModel, lsChannels, fjChannels, lsQuantumNumbers, fjQuantumNumbers, coreQuantumNumbers

export PARA,
    FMODEL_HIGHN_S35,
    FMODEL_HIGHN_S45,
    FMODEL_HIGHN_S55,
    FMODEL_LOWN_P45,
    FMODEL_HIGHN_P25,
    FMODEL_HIGHN_P35,
    FMODEL_HIGHN_P45,
    FMODEL_HIGHN_P55,
    FMODEL_HIGHN_P65,
    FMODEL_HIGHN_D15,
    FMODEL_HIGHN_D25,
    FMODEL_HIGHN_D35,
    FMODEL_HIGHN_D45,
    FMODEL_HIGHN_D55,
    FMODEL_HIGHN_D65,
    FMODEL_HIGHN_D75,
    FMODEL_HIGHN_F45

# Isotope data
THRESHOLDS = Dict(coreQuantumNumbers(0, 0.5, 4) => 45932.287373577, coreQuantumNumbers(0, 0.5, 5) => 45932.120512528)

PARA = Parameters(
    :Sr87,
    1822.888486192*86.9088774970, # nuclear mass
    4.5, # nuclear spin
    109736.62301005243, # Rydberg constant in 1/cm
    45932.1956, # weighted ionization threshold in 1/cm
    -0.03337220975052096, # hyperfine constant in 1/cm
    -1.0936030, # nuclear dipole
    THRESHOLDS,
)

# MQDT Models
# Defect data as compiled in [F Robicheaux 2019 J. Phys. B: At. Mol. Opt. Phys. 52 244001]
FMODEL_HIGHN_S35 = fModel(
    :Sr87,
    "S F=7/2, ν > 11",
    1,
    ["5sns 3S1"],
    # ["(6s1/2)(ns1/2)"],
    Bool[1],
    [3.370778 0.418 -0.3;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 4, 0, 0.5, 3.5)]),
    [-1;;],
)

FMODEL_HIGHN_S45 = fModel(
    :Sr87,
    "S F=9/2, ν > 11",
    2,
    ["5sns 1S0", "5sns 3S1"],
    # ["(6s1/2)(ns1/2)", "(6s1/2)(ns1/2)"],
    Bool[1, 1],
    [3.26896 -0.138 0.9; 3.370778 0.418 -0.3],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 0, 0, 0, 0, 0), lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 4, 0, 0.5, 4.5), fjQuantumNumbers(0.5, 0, 0.5, 5, 0, 0.5, 4.5)]),
    [3/2/sqrt(5) sqrt(11)/2/sqrt(5); sqrt(11)/2/sqrt(5) -3/2/sqrt(5)],
)

FMODEL_HIGHN_S55 = fModel(
    :Sr87,
    "S F=11/2, ν > 11",
    1,
    ["5sns 3S1"],
    # ["(6s1/2)(ns1/2)"],
    Bool[1],
    [3.370778 0.418 -0.3;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 5, 0, 0.5, 5.5)]),
    [1;;],
)

FMODEL_LOWN_P45 = fModel(
    :Sr87,
    "clock, 1.8 < ν < 2.2",
    4,
    ["5snp 1P1", "5snp 3P0", "5snp 3P1", "5snp 3P2"],
    # ["(6s1/2)(np1/2)", "(6s1/2)(np1/2)", "(6s1/2)(np3/2)", "(6s1/2)(np3/2)"],
    Bool[1, 1, 1, 1],
    [0.8720737 0; 0.13689075 0; 0.13143188 0; 0.11955235 0],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 4.5),
    ]),
    [
        -0.428174 -0.516398 0.387298 -0.632456;
        0.67082 0.0 0.74162 0.0;
        0.60553 -0.365148 -0.547723 -0.447214;
        0.0 -0.774597 0.0 0.632456
    ]',
)

FMODEL_HIGHN_P25 = fModel(
    :Sr87,
    "P F=5/2, ν > 5",
    1,
    ["5snp 3P2"],
    # ["(6s1/2)(np3/2)"],
    Bool[1],
    [2.882 0.446 -1.9;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 2.5)]),
    [-1;;],
)

FMODEL_HIGHN_P35 = fModel(
    :Sr87,
    "P F=7/2, ν > 5",
    3,
    ["5snp 1P1", "5snp 3P1", "5snp 3P2"],
    # ["(6s1/2)(np1/2)", "(6s1/2)(np3/2)", "(6s1/2)(np3/2)"],
    Bool[1, 1, 1],
    [2.724 -4.67 -157; 2.8826 0.407 -1.3; 2.882 0.446 -1.9],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 0.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 3.5),
    ]),
    [0.57735 0.341565 0.74162; -0.816497 0.241523 0.524404; 0.0 0.908295 -0.41833]',
)

FMODEL_HIGHN_P45 = fModel(
    :Sr87,
    "P F=9/2, ν > 7",
    4,
    ["5snp 1P1", "5snp 3P0", "5snp 3P1", "5snp 3P2"],
    # ["(6s1/2)(np1/2)", "(6s1/2)(np1/2)", "(6s1/2)(np3/2)", "(6s1/2)(np3/2)"],
    Bool[1, 1, 1, 1],
    [2.724 -4.67 -157; 2.8867 0.44 -1.9; 2.8826 0.407 -1.3; 2.882 0.446 -1.9],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 4.5),
    ]),
    [
        -0.428174 -0.516398 0.387298 -0.632456;
        0.67082 0.0 0.74162 0.0;
        0.60553 -0.365148 -0.547723 -0.447214;
        0.0 -0.774597 0.0 0.632456
    ]',
)

FMODEL_HIGHN_P55 = fModel(
    :Sr87,
    "P F=11/2, ν > 5",
    3,
    ["5snp 1P1", "5snp 3P1", "5snp 3P2"],
    # ["(6s1/2)(np3/2)", "(6s1/2)(np1/2)", "(6s1/2)(np3/2)"],
    Bool[1, 1, 1],
    [2.724 -4.67 -157; 2.8826 0.407 -1.3; 2.882 0.446 -1.9],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 0.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 5.5),
    ]),
    [0.67082 -0.57735 0.465475; 0.474342 0.816497 0.32914; 0.570088 0.0 -0.821584]',
)

FMODEL_HIGHN_P65 = fModel(
    :Sr87,
    "P F=13/2, ν > 5",
    1,
    ["5snp 3P2"],
    # ["(6s1/2)(np3/2)"],
    Bool[1],
    [2.882 0.446 -1.9;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 6.5)]),
    [1;;],
)

FMODEL_HIGHN_D15 = fModel(
    :Sr87,
    "D F=3/2, ν > 47",
    1,
    ["5snd 3D3"],
    # ["(6s1/2)(nd5/2)"],
    Bool[1],
    [2.655 -41.4 -15363;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 1.5)]),
    [-1;;],
)

FMODEL_HIGHN_D25 = fModel(
    :Sr87,
    "D F=5/2, ν > 47",
    3,
    ["5snd 1D2", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1],
    [2.3847 -39.41 -1090; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["12"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 1.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 2.5),
    ]),
    [0.632456 0.223607 0.74162; -0.774597 0.182574 0.60553; 0.0 0.957427 -0.288675]',
)

FMODEL_HIGHN_D35 = fModel(
    :Sr87,
    "D F=7/2, ν > 47",
    4,
    ["5snd 1D2", "5snd 3D1", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1, 1],
    [2.3847 -39.41 -1090; 2.67524 -13.15 -4444; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["13"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 1.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 1.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 3.5),
    ]),
    [
        -0.574456 -0.34641 0.264575 -0.69282;
        0.41833 0.0 0.908295 0.0;
        0.703562 -0.282843 -0.324037 -0.565685;
        0.0 -0.894427 0.0 0.447214
    ]',
)

FMODEL_HIGHN_D45 = fModel(
    :Sr87,
    "D F=9/2, ν > 47",
    4,
    ["5snd 1D2", "5snd 3D1", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1, 1],
    [2.3847 -39.41 -1090; 2.67524 -13.15 -4444; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["13"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 1.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 1.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 4.5),
    ]),
    [
        0.489898 0.458258 -0.4 0.6245;
        -0.632456 0.0 -0.774597 0.0;
        -0.6 0.374166 0.489898 0.509902;
        0.0 0.806226 0.0 -0.591608
    ]',
)

FMODEL_HIGHN_D55 = fModel(
    :Sr87,
    "D F=11/2, ν > 47",
    4,
    ["5snd 1D2", "5snd 3D1", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1, 1],
    [2.3847 -39.41 -1090; 2.67524 -13.15 -4444; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["13"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 1.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 1.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 5.5),
    ]),
    [
        -0.360555 -0.565685 0.519615 -0.52915;
        0.821584 0.0 0.570088 0.0;
        0.441588 -0.46188 -0.636396 -0.432049;
        0.0 -0.68313 0.0 0.730297
    ]',
)

FMODEL_HIGHN_D65 = fModel(
    :Sr87,
    "D F=13/2, ν > 47",
    3,
    ["5snd 1D2", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1],
    [2.3847 -39.41 -1090; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["12"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 6.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 1.5, 6.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 6.5),
    ]),
    [0.67082 -0.632456 0.387298; 0.547723 0.774597 0.316228; 0.5 0.0 -0.866025]',
)

FMODEL_HIGHN_D75 = fModel(
    :Sr87,
    "D F=15/2, ν > 47",
    1,
    ["5snd 3D3"],
    # ["(6s1/2)(nd5/2)"],
    Bool[1],
    [2.655 -41.4 -15363;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)]),
    fjChannels([fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 7.5)]),
    [1;;],
)

FMODEL_HIGHN_F45 = fModel(
    :Sr87,
    "F F=9/2, ν > 9",
    4,
    ["5snf 1F3", "5snf 3F2", "5snf 3F3", "5snf 3F4"],
    # ["(6s1/2)(nf5/2)", "(6s1/2)(nf5/2)", "(6s1/2)(nf7/2)", "(6s1/2)(nf7/2)"],
    Bool[1, 1, 1, 1],
    [0.12 -2.2 120; 0.12 -2.2 120; 0.12 -2.2 120; 0.12 -2.2 120],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 3, 3, 3),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 2),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 3),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 4),
    ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 3, 2.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 3, 3.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 3, 2.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 3, 3.5, 4.5),
    ]),
    [
        -0.527799 -0.414039 0.387298 -0.632456;
        0.591608 0.0 0.806226 0.0;
        0.609449 -0.358569 -0.447214 -0.547723;
        0.0 -0.83666 0.0 0.547723
    ]',
)

# TODO: other F state models

end
