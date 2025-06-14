module Sr87

using ..MQDT: Parameters, fModel, lsChannels, fjChannels, lsQuantumNumbers, fjQuantumNumbers, test_model

export PARA, MODEL_S35, MODEL_S45, MODEL_S55, MODEL_CL, MODEL_P25, MODEL_P35, MODEL_P45, MODEL_P55, MODEL_P65, MODEL_D15, MODEL_D25, MODEL_D35, MODEL_D45, MODEL_D55, MODEL_D65, MODEL_D75, MODEL_F45

# Isotope data 
PARA = Parameters(
    :Sr87,
    1822.888486192*86.9088774970, # nuclear mass
    4.5, # nuclear spin
    109736.62301005243, # Rydberg constant in 1/cm
    45932.1956, # weighted ionization threshold in 1/cm
    -0.03337220975052096, # hyperfine constant in 1/cm
    -1.0936030 # nuclear dipole
)

# MQDT Models
MODEL_S35 = fModel(
    "S F=7/2",
    1,
    ["5sns 3S1"],
    # ["(6s1/2)(ns1/2)"],
    Bool[1], 
    [45932.287373577],
    [3.370778 0.418 -0.3;],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 0, 0.5, 3.5)
        ]),
    [1;;]
)

MODEL_S45 = fModel(
    "S F=9/2",
    2,
    ["5sns 1S0", "5sns 3S1"],
    # ["(6s1/2)(ns1/2)", "(6s1/2)(ns1/2)"],
    Bool[1, 1], 
    [45932.287373577, 45932.120512528],
    [3.26896 -0.138 0.9; 3.370778 0.418 -0.3],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 0, 0, 1),
        lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 0, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 0, 0.5, 4.5)
        ]),
    [3/2/sqrt(5) sqrt(11)/2/sqrt(5); sqrt(11)/2/sqrt(5) -3/2/sqrt(5)]
)

MODEL_S55 = fModel(
    "S F=11/2",
    1,
    ["5sns 3S1"],
    # ["(6s1/2)(ns1/2)"],
    Bool[1], 
    [45932.120512528],
    [3.370778 0.418 -0.3;],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 5, 0, 0.5, 5.5)
        ]),
    [1;;]
)

MODEL_CL = fModel(
    "clock",
    4,
    ["5snp 1P1", "5snp 3P0", "5snp 3P1", "5snp 3P2"],
    # ["(6s1/2)(np1/2)", "(6s1/2)(np1/2)", "(6s1/2)(np3/2)", "(6s1/2)(np3/2)"],
    Bool[1, 1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528, 45932.120512528],
    [0.8720737 0; 0.13689075 0; 0.13143188 0; 0.11955235 0],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 4.5)
        ]),
    [-0.428174 -0.516398 0.387298 -0.632456; 0.67082 0.0 0.74162 0.0; 0.60553 -0.365148 -0.547723 -0.447214; 0.0 -0.774597 0.0 0.632456]'
)

MODEL_P25 = fModel(
    "P F=5/2",
    1,
    ["5snp 3P2"],
    # ["(6s1/2)(np3/2)"],
    Bool[1], 
    [45932.287373577],
    [2.882 0.446 -1.9;],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 2.5)
        ]),
    [1;;]
)

MODEL_P35 = fModel(
    "P F=7/2",
    3,
    ["5snp 1P1", "5snp 3P1", "5snp 3P2"],
    # ["(6s1/2)(np1/2)", "(6s1/2)(np3/2)", "(6s1/2)(np3/2)"],
    Bool[1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528],
    [2.724 -4.67 -157; 2.8826 0.407 -1.3; 2.882 0.446 -1.9],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 0.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 3.5)
        ]),
    [0.57735 0.341565 0.74162; -0.816497 0.241523 0.524404; 0.0 0.908295 -0.41833]'
)

MODEL_P45 = fModel(
    "P F=9/2",
    4,
    ["5snp 1P1", "5snp 3P0", "5snp 3P1", "5snp 3P2"],
    # ["(6s1/2)(np1/2)", "(6s1/2)(np1/2)", "(6s1/2)(np3/2)", "(6s1/2)(np3/2)"],
    Bool[1, 1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528, 45932.120512528],
    [2.724 -4.67 -157; 2.8867 0.44 -1.9; 2.8826 0.407 -1.3; 2.882 0.446 -1.9],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 0.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 4.5)
        ]),
    [-0.428174 -0.516398 0.387298 -0.632456; 0.67082 0.0 0.74162 0.0; 0.60553 -0.365148 -0.547723 -0.447214; 0.0 -0.774597 0.0 0.632456]'
)

MODEL_P55 = fModel(
    "P F=11/2",
    3,
    ["5snp 1P1", "5snp 3P1", "5snp 3P2"],
    # ["(6s1/2)(np3/2)", "(6s1/2)(np1/2)", "(6s1/2)(np3/2)"],
    Bool[1, 1, 1], 
    [45932.287373577, 45932.120512528, 45932.120512528],
    [2.724 -4.67 -157; 2.8826 0.407 -1.3; 2.882 0.446 -1.9],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 1, 1.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 0.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 5.5)
        ]),
    [0.67082 -0.57735 0.465475; 0.474342 0.816497 0.32914; 0.570088 0.0 -0.821584]'
)

MODEL_P65 = fModel(
    "P F=13/2",
    1,
    ["5snp 3P2"],
    # ["(6s1/2)(np3/2)"],
    Bool[1], 
    [45932.120512528],
    [2.882 0.446 -1.9;],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 5, 1, 1.5, 6.5)
        ]),
    [1;;]
)

MODEL_D15 = fModel(
    "D F=3/2",
    1,
    ["5snd 3D3"],
    # ["(6s1/2)(nd5/2)"],
    Bool[1], 
    [45932.120512528],
    [2.655 -41.4 -15363;],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 1.5)
        ]),
    [1;;]
)

MODEL_D25 = fModel(
    "D F=5/2",
    3,
    ["5snd 1D2", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528],
    [2.3847 -39.41 -1090; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["12"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 1.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 2.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 2.5)
        ]),
    [0.632456 0.223607 0.74162; -0.774597 0.182574 0.60553; 0.0 0.957427 -0.288675]'
)

MODEL_D35 = fModel(
    "D F=7/2",
    4,
    ["5snd 1D2", "5snd 3D1", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528, 45932.120512528],
    [2.3847 -39.41 -1090; 2.67524 -13.15 -4444; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["13"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 1.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 1.5, 3.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 3.5)
        ]),
    [-0.574456 -0.34641 0.264575 -0.69282; 0.41833 0.0 0.908295 0.0; 0.703562 -0.282843 -0.324037 -0.565685; 0.0 -0.894427 0.0 0.447214]'
)

MODEL_D45 = fModel(
    "D F=9/2",
    4,
    ["5snd 1D2", "5snd 3D1", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528, 45932.120512528],
    [2.3847 -39.41 -1090; 2.67524 -13.15 -4444; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["13"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 1.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 1.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 4.5)
        ]),
    [0.489898 0.458258 -0.4 0.6245; -0.632456 0.0 -0.774597 0.0; -0.6 0.374166 0.489898 0.509902; 0.0 0.806226 0.0 -0.591608]'
)

MODEL_D55 = fModel(
    "D F=11/2",
    4,
    ["5snd 1D2", "5snd 3D1", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528, 45932.120512528],
    [2.3847 -39.41 -1090; 2.67524 -13.15 -4444; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["13"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 1.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 1.5, 5.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 5.5)
        ]),
    [-0.360555 -0.565685 0.519615 -0.52915; 0.821584 0.0 0.570088 0.0; 0.441588 -0.46188 -0.636396 -0.432049; 0.0 -0.68313 0.0 0.730297]'
)

MODEL_D65 = fModel(
    "D F=13/2",
    3,
    ["5snd 1D2", "5snd 3D2", "5snd 3D3"],
    # ["(6s1/2)(nd3/2)", "(6s1/2)(nd5/2)", "(6s1/2)(nd5/2)"],
    Bool[1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528],
    [2.3847 -39.41 -1090; 2.66149 -16.77 -6656; 2.655 -41.4 -15363],
    ["12"],
    [-0.14 0;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 2, 2.5, 6.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 1.5, 6.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 6.5)
        ]),
    [0.67082 -0.632456 0.387298; 0.547723 0.774597 0.316228; 0.5 0.0 -0.866025]'
)

MODEL_D75 = fModel(
    "D F=15/2",
    1,
    ["5snd 3D3"],
    # ["(6s1/2)(nd5/2)"],
    Bool[1], 
    [45932.120512528],
    [2.655 -41.4 -15363;],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 5, 2, 2.5, 7.5)
        ]),
    [1;;]
)

MODEL_F45 = fModel(
    "F F=9/2",
    4,
    ["5snf 1F3", "5snf 3F2", "5snf 3F3", "5snf 3F4"],
    # ["(6s1/2)(nf5/2)", "(6s1/2)(nf5/2)", "(6s1/2)(nf7/2)", "(6s1/2)(nf7/2)"],
    Bool[1, 1, 1, 1], 
    [45932.287373577, 45932.287373577, 45932.120512528, 45932.120512528],
    [0.12 -2.2 120; 0.12 -2.2 120; 0.12 -2.2 120; 0.12 -2.2 120],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 3, 3, 3),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 2),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 3),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 4)
        ]),
    fjChannels([
        fjQuantumNumbers(0.5, 0, 0.5, 4, 3, 2.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 4, 3, 3.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 3, 2.5, 4.5),
        fjQuantumNumbers(0.5, 0, 0.5, 5, 3, 3.5, 4.5)
        ]),
    [-0.527799 -0.414039 0.387298 -0.632456; 0.591608 0.0 0.806226 0.0; 0.609449 -0.358569 -0.447214 -0.547723; 0.0 -0.83666 0.0 0.547723]'
)

# TODO: other F state models

MODELS = [
    MODEL_CL,
    MODEL_S35,
    MODEL_S45,
    MODEL_S55,
    MODEL_P25,
    MODEL_P35,
    MODEL_P45,
    MODEL_P55,
    MODEL_P65,
    MODEL_D15,
    MODEL_D25,
    MODEL_D35,
    MODEL_D45,
    MODEL_D55,
    MODEL_D65,
    MODEL_D75
]
test_model(MODELS)

end