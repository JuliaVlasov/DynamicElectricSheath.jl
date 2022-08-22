# DynamicElectricSheath.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaVlasov.github.io/DynamicElectricSheath.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaVlasov.github.io/DynamicElectricSheath.jl/dev/)
[![Build Status](https://github.com/JuliaVlasov/DynamicElectricSheath.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaVlasov/DynamicElectricSheath.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaVlasov/DynamicElectricSheath.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaVlasov/DynamicElectricSheath.jl)

Two-species Vlasov-Poisson solver on $[-1,1] \times \RR$ using
 
 with :
 - non-emitting boundary conditions on the distribution functions
 - Source term for the ions
 - Dynamic electric field at the boundary
 
 
 Simulation of the symmetric sheath problem with ionization term.
 
 Discretization :
 - Upwind for the transport
 - Integration for the Poisson equation -\\lambda^2 dx E = rho with trapeze formula for rho.  The integration  preserves the symmetry of rho if it is even.
 
 
 Initial data :
 
 - Initial data are in the file initial_data.hpp
 
  * Author : Mehdi BASDI. 
  * Date : 27/07/2022.
  * Translation in Julia : Pierre Navaro (@pnavaro) & Averil Prost (@averil-prost).
  * Date : 05/08/2022.

```bash
git clone https://github.com/juliavlasov/DynamicElectricSheath.jl
cd DynamicElectricSheath.jl
julia --project
julia> using Pkg
julia> Pkg.instantiate()
julia> include("example.jl")
```
