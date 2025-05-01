# MQDT

[![Build Status](https://github.com/atom-pairinteraction/MQDT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/atom-pairinteraction/MQDT.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Running tests from within the repository

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Running examples from within the repository

```bash
julia --project=examples -e 'import Pkg; Pkg.instantiate(); include("examples/generate_lu_fano.jl")'
```
