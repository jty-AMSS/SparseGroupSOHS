# SparseGroupSOHS

This repository contains Julia code for the numerical experiments in the paper: Sparse Sum of Hermitian Squares in Group Algebras of Finite Groups.

The current implementation focuses on cyclic groups $ \mathbb{Z}_N $. It generates random self-adjoint elements in the group algebra $ \mathbb{C}[\mathbb{Z}_N] $, rescales them to satisfy $ 1 \succeq f \succeq  \alpha $, and solves semidefinite feasibility problems to find sparse Fourier sum-of-squares certificates.

The experiments were reproduced in Linux. To reproduce the experiment, start Julia and execute `include("Ex0.jl")`.
## Files
- `ZN.jl`: Basic implementation of cyclic group elements and elements of   $ \mathbb{C}[\mathbb{Z}_N] $, including algebra operations and Fourier transforms.
- `Ex0.jl`: Numerical experiment for the sparse FSOS hierarchy over  $ Z_N$.

## Requirements

The code is written in Julia and uses the following packages:

```julia
using JuMP
using COSMO
using FFTW
using LinearAlgebra
using Random
