module Sr88

using ..MQDT:
    Parameters,
    fModel,
    kModel,
    lsChannels,
    jjChannels,
    lsQuantumNumbers,
    jjQuantumNumbers,
    test_model

export PARA,
    FMODEL_LOWN_P1,
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
    FMODEL_HIGHN_F4

# Isotope data
PARA = Parameters(
    :Sr88,
    1822.888486192*87.9056122571, # nuclear mass
    0, # nuclear spin
    109736.63086399352, # Rydberg constant in 1/cm
    # 45932.2002, # threshold used in [JPB 47 155001 (2014)]
    45932.1956, # lowest ionization threshold in 1/cm
    0, # hyperfine constant in 1/cm
    2.3, # nuclear dipole
)

# MQDT Models
FMODEL_HIGHN_S0 = fModel(
    :Sr88,
    "S J=0",
    1,
    ["5sns 1S0"],
    Bool[1],
    [45932.1956],
    [3.26896 -0.138 0.9;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 0, 0, 0, 0, 0)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 0)]),
    [1;;],
)

FMODEL_HIGHN_S1 = fModel(
    :Sr88,
    "S J=1",
    1,
    ["5sns 3S1"],
    Bool[1],
    [45932.1956],
    [3.370778 0.418 -0.3;],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 0, 0, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 1)]),
    [1;;],
)

FMODEL_LOWN_P1 = fModel(
    :Sr88,
    "recombination",
    2,
    ["5snp 1P1", "5snp 3P1"],
    Bool[1, 1],
    [45932.1956, 45932.1956],
    [0.87199081 0; 0.13140955 0],
    [""],
    [1.31169947 -4.48280597;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 1),
        jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 1),
    ]),
    [-sqrt(1/3) sqrt(2/3); sqrt(2/3) sqrt(1/3)],
)

FMODEL_HIGHN_P0 = fModel(
    :Sr88,
    "P J=0",
    1,
    ["5snp 3P0"],
    Bool[1],
    [45932.1956],
    [2.8867 0.44 -1.9],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 0)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 0)]),
    [1;;],
)

FMODEL_HIGHN_P1 = fModel(
    :Sr88,
    "P J=1",
    2,
    ["5snp 1P1", "5snp 3P1"],
    Bool[1, 1],
    [45932.1956, 45932.1956],
    [2.724 -4.67 -157; 2.8826 0.407 -1.3],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 1),
        jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 1),
    ]),
    [-sqrt(1/3) sqrt(2/3); sqrt(2/3) sqrt(1/3)],
)

FMODEL_HIGHN_P2 = fModel(
    :Sr88,
    "P J=2",
    1,
    ["5snp 3P2"],
    Bool[1],
    [45932.1956],
    [2.882 0.446 -1.9],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 1, 1, 2)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 2)]),
    [1;;],
)

FMODEL_HIGHN_D1 = fModel(
    :Sr88,
    "D J=1",
    1,
    ["5snd 3D1"],
    Bool[1],
    [45932.1956],
    [2.67524 -13.15 -4444],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 1)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 1)]),
    [1;;],
)

FMODEL_HIGHN_D2 = fModel(
    :Sr88,
    "D J=2",
    2,
    ["5snd 1D2", "5snd 3D2"],
    Bool[1, 1],
    [45932.1956, 45932.1956],
    [2.3847 -39.41 -1090; 2.66149 -16.77 -6656],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 2),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 2),
        jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 2),
    ]),
    [-sqrt(2/5) sqrt(3/5); sqrt(3/5) sqrt(2/5)],
)

FMODEL_HIGHN_D3 = fModel(
    :Sr88,
    "D J=3",
    1,
    ["5snd 3D3"],
    Bool[1],
    [45932.1956],
    [2.655 -41.4 -15363],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 2, 2, 3)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 3)]),
    [1;;],
)

FMODEL_HIGHN_F2 = fModel(
    :Sr88,
    "F J=2",
    1,
    ["5snd 3F2"],
    Bool[1],
    [45932.1956],
    [0.12 -2.2 120],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 3, 3, 2)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 3, 2.5, 2)]),
    [1;;],
)

FMODEL_HIGHN_F3 = fModel(
    :Sr88,
    "F J=3",
    2,
    ["5snd 1F3", "5snd 3F3"],
    Bool[1, 1],
    [45932.1956, 45932.1956],
    [0.089 -2 30; 0.12 -2.2 120],
    [""],
    [0;;],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 3, 3, 3),
        lsQuantumNumbers(0.5, 1, 0, 3, 3, 3),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 3, 2.5, 3),
        jjQuantumNumbers(0.5, 0, 0.5, 3, 3.5, 3),
    ]),
    [-sqrt(3/7) sqrt(4/7); sqrt(4/7) sqrt(3/7)],
)

FMODEL_HIGHN_F4 = fModel(
    :Sr88,
    "F J=4",
    1,
    ["5snd 3F4"],
    Bool[1],
    [45932.1956],
    [0.12 -2.2 120],
    [""],
    [0;;],
    lsChannels([lsQuantumNumbers(0.5, 1, 0, 3, 3, 4)]),
    jjChannels([jjQuantumNumbers(0.5, 0, 0.5, 3, 3.5, 4)]),
    [1;;],
)

# FMODELS = [
#     FMODEL_LOWN_P1,
#     FMODEL_HIGHN_S0,
#     FMODEL_HIGHN_S1,
#     FMODEL_HIGHN_P0,
#     FMODEL_HIGHN_P1,
#     FMODEL_HIGHN_P2,
#     FMODEL_HIGHN_D1,
#     FMODEL_HIGHN_D2,
#     FMODEL_HIGHN_D3,
#     FMODEL_HIGHN_F2,
#     FMODEL_HIGHN_F3,
#     FMODEL_HIGHN_F4,
# ]
# test_model(FMODELS)

KMODEL_S0 = kModel(
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
)

KMODEL_S1 = kModel(
    :Sr88,
    "3S1",
    2,
    ["5sns 3S1", "5pnp 3P1"],
    Bool[0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 0, 0, 1),
        lsQuantumNumbers(0.5, 1, 1, 1, 1, 1),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 0, 0.5, 1),
        jjQuantumNumbers(0.5, 1, NaN, 1, NaN, 1),
    ]),
    [45932.2002, 70048.11],
    [-103.9244 -133.4517; -133.4517 -168.0452],
    [-27.66912, 55.17184],
)

KMODEL_P0 = kModel(
    :Sr88,
    "3P0",
    2,
    ["5snp 3P0", "4dnp 3P0"],
    Bool[0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 0),
        lsQuantumNumbers(0.5, 1, 2, 1, 1, 0),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 1, 0.5, 0),
        jjQuantumNumbers(0.5, 2, 1.5, 1, 1.5, 0),
    ]),
    [45932.2002, 60628.26],
    [-0.4009565 -0.2220569; -0.2220569 0.4025180],
    [1.039923, -1.021696],
)

KMODEL_1P1 = kModel(
    :Sr88,
    "1P1",
    2,
    ["5snp 1P1", "4dnp 1P1"],
    Bool[0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 0, 2, 1, 1, 1),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 1, NaN, 1),
        jjQuantumNumbers(0.5, 2, NaN, 1, NaN, 1),
    ]),
    [45932.2002, 60628.26],
    [11.16809 16.16933; 16.16933 22.39617],
    [-0.9097862, 4.272626],
)

KMODEL_3P1 = kModel(
    :Sr88,
    "3P1",
    2,
    ["5snp 3P1", "4dnp 3P1"],
    Bool[0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 1),
        lsQuantumNumbers(0.5, 1, 2, 1, 1, 1),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 1, NaN, 1),
        jjQuantumNumbers(0.5, 2, NaN, 1, NaN, 1),
    ]),
    [45932.2002, 60628.26],
    [-0.4199067 -0.2292304; -0.2292304 -0.3526179],
    [1.082615, -1.304779],
)

KMODEL_P2 = kModel(
    :Sr88,
    "3P2",
    2,
    ["5snp 3P2", "4dnp 3P2"],
    Bool[0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 1, 1, 2),
        lsQuantumNumbers(0.5, 1, 2, 1, 1, 2),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 1, 1.5, 2),
        jjQuantumNumbers(0.5, 2, NaN, 1, NaN, 2),
    ]),
    [45932.2002, 60628.26],
    [-0.4531133 -0.2179619; -0.2179619 -0.5285102],
    [1.050866, -0.4051199],
)

KMODEL_D1 = kModel(
    :Sr88,
    "3D1",
    2,
    ["5snd 3D1", "4dns 3D1"],
    Bool[0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 1),
        lsQuantumNumbers(0.5, 1, 2, 0, 2, 1),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 1),
        jjQuantumNumbers(0.5, 2, 1.5, 0, 0.5, 1),
    ]),
    [45932.2002, 60628.26],
    [-0.7403359 0.5504572; 0.5504572 1.461400],
    [0.9684681, 0.2777353],
)

KMODEL_D2 = kModel(
    :Sr88,
    "1D2 / 3D2",
    6,
    [
        "(5s1/2)(nd5/2)",
        "(5s1/2)(nd3/2)",
        "(4d5/2)(ns1/2)",
        "(4d3/2)(ns1/2)",
        "5pnp 1D2",
        "4dnd 3P2",
    ],
    Bool[1, 1, 1, 1, 0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, NaN, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, NaN, 0, 2, 2, 2),
        lsQuantumNumbers(0.5, NaN, 2, 0, 2, 2),
        lsQuantumNumbers(0.5, NaN, 2, 0, 2, 2),
        lsQuantumNumbers(0.5, 0, 1, 1, 2, 2),
        lsQuantumNumbers(0.5, 1, 2, 2, 1, 2),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 2),
        jjQuantumNumbers(0.5, 0, 0.5, 2, 1.5, 2),
        jjQuantumNumbers(0.5, 2, 2.5, 0, 0.5, 2),
        jjQuantumNumbers(0.5, 2, 1.5, 0, 0.5, 2),
        jjQuantumNumbers(0.5, 1, NaN, 1, NaN, 2),
        jjQuantumNumbers(0.5, 2, NaN, 2, NaN, 2),
    ]),
    [45932.2002, 45932.2002, 60628.26, 60628.26, 70048.11, 60628.26],
    [
        -0.3853883 0.2308103 -0.2996898 0.6248391 -0.2381621 -0.08944624;
        0.2308103 -0.4881877 -0.6411698 0.000008101262 -0.4849582 0.002427350;
        -0.2996898 -0.6411698 1.136225 0.2078805 0.0 0.0;
        0.6248391 0.000008101262 0.2078805 1.123831 0.0 0.0;
        -0.2381621 -0.4849582 0.0 0.0 0.6117878 0.0;
        -0.08944624 0.002427350 0.0 0.0 0.0 2.205400
    ],
    [-1.775326, 2.052554, 4.733804, 3.989162, 5.292869, 6.079562],
)

KMODEL_D3 = kModel(
    :Sr88,
    "3D3",
    3,
    ["5snd 3D3", "4dns 3D3", "4dnd 3D3"],
    Bool[0, 0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, 1, 0, 2, 2, 3),
        lsQuantumNumbers(0.5, 1, 2, 0, 2, 3),
        lsQuantumNumbers(0.5, 1, 2, 2, 2, 3),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 2, 2.5, 3),
        jjQuantumNumbers(0.5, 2, 2.5, 0, 0.5, 3),
        jjQuantumNumbers(0.5, 2, NaN, 2, NaN, 3),
    ]),
    [45932.2002, 60628.26, 60628.26],
    [
        -0.7793857 0.4360198 0.2229788;
        0.4360198 1.212314 -0.0001683225;
        0.2229788 -0.0001683225 -0.2238265
    ],
    [1.071997, 8.514161, 5.544426],
)

KMODEL_1F3 = kModel(
    :Sr88,
    "1F3",
    2,
    ["5snf 1F3", "4dnp 1F3"],
    Bool[0, 0],
    lsChannels([
        lsQuantumNumbers(0.5, 0, 0, 3, 3, 3),
        lsQuantumNumbers(0.5, 0, 2, 1, 3, 3),
    ]),
    jjChannels([
        jjQuantumNumbers(0.5, 0, 0.5, 3, NaN, 3),
        jjQuantumNumbers(0.5, 2, NaN, 1, NaN, 3),
    ]),
    [45932.2002, 60628.26],
    [0.1711631 0.4505951; 0.4505951 -0.6978294],
    [-0.3530368, -1.318505],
)

# KMODELS = [
#     KMODEL_S0,
#     KMODEL_S1,
#     KMODEL_P0,
#     KMODEL_1P1,
#     KMODEL_3P1,
#     KMODEL_P2,
#     KMODEL_D1,
#     KMODEL_D2,
#     KMODEL_D3,
#     KMODEL_1F3,
# ]
# test_model(KMODELS)

end
