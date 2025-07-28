# MQDT

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pairinteraction.github.io/MQDT.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pairinteraction.github.io/MQDT.jl/dev/)
[![Build Status](https://github.com/pairinteraction/MQDT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pairinteraction/MQDT.jl/actions/workflows/CI.yml?query=branch%3Amain)

MQDT.jl is a Julia package designed to calculate Rydberg states and matrix elements using multi-channel quantum defect theory (MQDT).

## Installing Julia

```bash
curl -fsSL https://install.julialang.org | sh
```

## Update the package dependencies

```bash
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Building docs

From within the `docs` directory, you can build the documentation with:
```bash
julia --project=. -e 'import Pkg; Pkg.develop(path=".."); Pkg.instantiate(); include("make.jl")'
```

## Running examples

From within the `examples` directory, first initialize a new julia project and install the MQDT.jl package in development mode:
```bash
julia --project=. -e 'import Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
```

Then you can run a script (from within the `examples` directory) like this:
```bash
julia --project=. generate_lu_fano.jl
```

## Citation

If you use MQDT.jl in your research, please consider citing our paper:

> Frederic Hummel, Sebastian Weber, Johannes Moegerle, Henri Menke, Jonathan King, Benjamin Bloom, Sebastian Hofferberth, Ming Li, *Engineering Rydberg-pair interactions in divalent atoms with hyperfine-split ionization thresholds*, [Phys. Rev. A 110, 042821 (2024)][journal-link]

[journal-link]: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.110.042821
