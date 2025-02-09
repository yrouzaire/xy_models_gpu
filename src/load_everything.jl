# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

using Parameters, BenchmarkTools
using JLD2
using CUDA
using Distributions
using FFTW # for the Fourier Transform of the spatial correlation function

Tf = Float32


include("init.jl")
include("measurements.jl")
include("auxiliary.jl")

include("update_XY.jl")
include("update_NRXY.jl")
include("update_Kuramoto.jl")
include("update_NRKuramoto.jl")

include("load_CairoMakie.jl")
include("visualisation.jl")
# using GLMakie 
# GLMakie is only useful for movies, so only loaded in the corresponding file to avoid conflicts with CairoMakie


"All Good !"