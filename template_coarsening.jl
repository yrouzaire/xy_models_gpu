# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include("src/load_everything.jl")

#= 
This file is a template for a simple simulation of the coarsening of the Kuramoto model.
It shows how to 
    - define the parameters 
    - run the simulation (update the system) 
    - perform measurements on the system 
    - save the results in a .jld2 file (and, separately, the configurations themselves) 
    - load the data back 
=#

## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##

lattice_type = "triangular"

Lx = 256
Ly = 256
R = 4 # number of independent realisations

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)

tmax = Tf(1E2)
dt = Tf(1E-1)
time_spacing = "linear" # choose among "linear", "log", "quadratic"
if time_spacing == "linear"
    every = tmax / 50
    times = Tf.(collect(every:every:tmax)) # linear
elseif time_spacing == "log"
    times = logspace(10dt, tmax, nb_frames) # log
elseif time_spacing == "quadratic"
    times = Tf.(range(sqrt(dt), sqrt(tmax), nb_frames)) .^ 2  # quadratic
else
    error("time_spacing should be 'linear', 'log' or 'quadratic'")
end

inits = ["hightemp"]
Ts = [0.05]
sigmas = [0.1]
distribution_types = ["uniform", "gaussian"]

thetas_saved_cpu = zeros(Float16, Lx, Ly, R, length(Ts), length(sigmas), length(inits), length(distribution_types), length(times)) # in Float16 for storage reasons
magnetisations = zeros(R, length(Ts), length(sigmas), length(inits), length(distribution_types), length(times))
number_of_defects = zeros(R, length(Ts), length(sigmas), length(inits), length(distribution_types), length(times))

m = 0
M = length(Ts) * length(inits) * length(sigmas) * length(distribution_types)

z = @elapsed for i in each(Ts), j in each(sigmas), k in each(inits), mm in each(distribution_types)

    T = Tf(Ts[i])
    sigma = Tf(sigmas[j])
    init = inits[k]
    distribution_type = distribution_types[mm]
    m += 1

    println("Simulation $m/$M : T = $T, σ = $(sigma), Init = $init, ω ~ $distribution_type")

    thetas = create_thetas(Lx, Ly, R, init)
    thetas_new = similar(thetas) # allocates a similar array to thetas 
    omegas = create_omegas(Lx, Ly, R, sigma, distribution_type)

    t = Float64(0)
    for tt in 1:length(times)
        thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, times[tt], lattice_type)
        
        thetas_saved_cpu[:, :, :, i, j, k, mm, tt] = Array(thetas)
        # magnetisations[:, i, j, k, mm, tt] = OP(thetas)[1]
        # number_of_defects[:, i, j, k, mm, tt] = number_defects(thetas)
    end
end
println("Simulation over.")
prinz(z) # prints the runtime in seconds, minutes and hours. Defined in `auxiliary.jl` file.

# plot_thetas(thetas_saved_cpu[:, :, 1, 1, 1, 1, 1, end]) # for a quick visual check if needed 

## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##

comments = "$(time_spacing)ly spaced times. $(lattice_type) lattice."

pwd()
filepath = pwd() * "/models/kuramoto/data/"
filename = "coarsening_Lx$(Lx)_Ly$(Ly)_R$(R)_Ts$(Ts)_sigmas$(sigmas)_inits_$(join(inits, "_"))_distributions_$(join(distribution_types, "_"))_tmax$(tmax)"

@save filepath * filename * ".jld2" Lx Ly R Ts sigmas times tmax dt inits distribution_types comments runtime = z
@save filepath * filename * "_with_thetas.jld2" thetas_saved_cpu Lx Ly R Ts sigmas tmax times dt inits distribution_types comments runtime = z

println("Data saved in $(filepath * filename).jld2")


## ------------------ If you want to load back the data ------------------ ##
filepath = pwd() * "/models/kuramoto/data/" # same as above
filename = "has-it-relaxed_coarsening_Lx1024_Ly1024_R1_Ts[0.05]_sigmas[0.1]_inits_hightemp_distributions_gaussian_tmax100000.0_with_thetas"
# @load filepath * filename * ".jld2" Lx Ly R Ts sigmas times tmax dt inits distribution_types comments runtime # without the "=z"
@load filepath * filename * ".jld2" thetas_saved_cpu Lx Ly R Ts sigmas tmax times dt inits distribution_types comments runtime

thetas_saved_cpu

plot_thetas(thetas_saved_cpu[:, :, :, 1, 1, 1, 1, end], defects=false) # for a quick visual check if needed 

