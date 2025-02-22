# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include("src/load_everything.jl")

## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##

lattice_type = "square"

Lx = Ly = 512
R = 8 # number of independent realisations
R_simus = 2
R_tot = R * R_simus

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)

tmax = Tf(1E4)
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
params_init = ("dummy") # only useful for the "pair" init (a pair of defects)

#= Parameters (α, σ) (intrinsic frequency, non-reciprocity)
(0,0) is the XY model 
(0,0.15) is the NRXY model 
(alpha_L2,0) is the Kuramoto model 
(alpha_L2,0.15) is the non-reciprocal Kuramoto model 
=#
alpha_L2 = 2 * 2.3 / Lx # in my Frontiers in Physics, I've shown that αξ = 2.3 for t → ∞. If we want ξ ≈ L/2, then one obtains α = 2 * 2.3 / L
alphas_sigmas = [(0, 0), (0, 0.15), (alpha_L2, 0), (alpha_L2, 0.15)]
ratios_TKT = [0.2, 0.8]
distribution_types = ["gaussian"]

thetas_saved_cpu_final_time = zeros(Float16, Lx, Ly, R, length(ratios_TKT), length(alphas_sigmas), length(inits), length(distribution_types)) # in Float16 for storage reasons
magnetisations = zeros(R_tot, length(ratios_TKT), length(alphas_sigmas), length(inits), length(distribution_types), length(times))
number_of_defects = zeros(R_tot, length(ratios_TKT), length(alphas_sigmas), length(inits), length(distribution_types), length(times))
xs = Array{Vector{Float64}}(undef, R_tot, length(ratios_TKT), length(alphas_sigmas), length(inits), length(distribution_types), length(times)) # positions of the defects
ys = Array{Vector{Float64}}(undef, R_tot, length(ratios_TKT), length(alphas_sigmas), length(inits), length(distribution_types), length(times)) # positions of the defects
qs = Array{Vector{Float64}}(undef, R_tot, length(ratios_TKT), length(alphas_sigmas), length(inits), length(distribution_types), length(times)) # positions of the defects


m = 0
M = length(ratios_TKT) * length(inits) * length(alphas_sigmas) * length(distribution_types) * R_simus

z = @elapsed for i in each(ratios_TKT), j in each(alphas_sigmas), k in each(inits), mm in each(distribution_types), r in 1:R_simus

    alpha, sigma = Tf.(alphas_sigmas[j])
    T = Tf(round(ratios_TKT[i] *  0.22 * (1 - 0.43 * alpha * log(Lx)), digits=3))
        # in the XY model, TKT_xy = 0.89/4 = 0.22 is the critical temperature of the phase transition to disorder.
        # in the Kuramoto model, T_c = TKT_xy * (1 - 0.43 α log(L)) is the critical temperature below which the defects are bounded
    init = inits[k]
    distribution_type = distribution_types[mm]
    
    global m += 1
    println("Simulation $m/$M : T = $T, α = $(alpha), σ = $(sigma), init = $init, ω ~ $distribution_type")

    thetas = create_thetas(Lx, Ly, R, init, params_init)
    thetas_new = similar(thetas) 
    omegas = create_omegas(Lx, Ly, R, alpha, distribution_type)

    t = Float64(0)
    for tt in eachindex(times)
        thetas, t = evolve_NRKuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, t, dt, times[tt], lattice_type)
        
        # thetas_saved_cpu[:, :, :, i, j, k, mm, tt] = Array(thetas)
        # magnetisations[:, i, j, k, mm, tt] = OP(thetas)[1]
        
        dico = spot_defects(thetas) # for a quick visual check if needed
        qs[1+(r-1)*R:(r*R), i, j, k, mm, tt] = dico["qs"]
        xs[1+(r-1)*R:(r*R), i, j, k, mm, tt] = dico["xs"]
        ys[1+(r-1)*R:(r*R), i, j, k, mm, tt] = dico["ys"]
        number_of_defects[1+(r-1)*R:(r*R), i, j, k, mm, tt] = dico["number_defects"]
    end
    
    thetas_saved_cpu_final_time[:, :, :, i, j, k, mm] = Array(thetas)
end
println("Simulation over.")
prinz(z) # prints the runtime in seconds, minutes and hours. Defined in `auxiliary.jl` file

# plot_thetas(thetas_saved_cpu_final_time[:, :, 1, 1, 1, 1, 1], defects=true) # for a quick visual check if needed 

## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##

comments = "$(time_spacing)ly spaced times. $(lattice_type) lattice."

pwd()
filepath = pwd() * "/models/nrkuramoto/data/"
filename = "test_hyperuniformity_Lx$(Lx)_Ly$(Ly)_Rtot$(R_tot)_ratios_TKT$(ratios_TKT)_inits_$(join(inits, "_"))_distributions_$(join(distribution_types, "_"))_tmax$(tmax)"

@save filepath * filename * ".jld2" number_of_defects qs xs ys Lx Ly R R_simus R_tot ratios_TKT alphas_sigmas times tmax dt inits distribution_types comments runtime = z
# @save filepath * filename * "_with_thetas_final_time.jld2" thetas_saved_cpu_final_time Lx Ly R R_simus R_tot ratios_TKT alphas_sigmas tmax times dt inits distribution_types comments runtime = z

println("Data saved in $(filepath * filename).jld2")
