# MQDT

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pairinteraction.github.io/MQDT.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pairinteraction.github.io/MQDT.jl/dev/)
[![Build Status](https://github.com/atom-pairinteraction/MQDT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/atom-pairinteraction/MQDT.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Installing Julia

```bash
curl -fsSL https://install.julialang.org | sh
```

## Running tests from within the repository

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Building docs from within the repository

```bash
julia --project=docs -e 'import Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
```

## Running examples from within the repository

```bash
julia --project=examples -e 'import Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("examples/generate_lu_fano.jl")'
```
