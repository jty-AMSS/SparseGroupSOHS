# SparseGroupSOHS

This repository contains Julia code for the numerical experiments in the paper: Sparse Sum of Hermitian Squares in Group Algebras of Finite Groups.

The current implementation focuses on cyclic groups $ Z_N $. It generates random self-adjoint elements in the group algebra $ \mathbb{C}[Z_N] $, rescales them to satisfy $ \alpha \succeq f \succeq1$, and solves semidefinite feasibility problems to find sparse Fourier sum-of-squares certificates.

## Files
- `ZN.jl`: Basic implementation of cyclic group elements and elements of   $ \mathbb{C}[Z_N] $, including algebra operations and Fourier transforms.
- `Ex0.jl`: Numerical experiment for the sparse FSOShierarchy over  $ Z_N$.

## Requirements

The code is written in Julia and uses the following packages:

```julia
using JuMP
using COSMO
using FFTW
using LinearAlgebra
using Random
