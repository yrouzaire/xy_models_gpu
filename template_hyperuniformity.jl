# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include("src/load_everything.jl")

## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##

lattice_type = "square"

Lx = 256
Ly = 256
@assert Lx == Ly "for this script, Lx must be equal to Ly"
R = 2 # number of independent realisations

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)

tmax = Tf(1E3)
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
params_init = (x_plus=Lx / 4, y_plus=Ly / 2, r0=Lx / 2, mu_plus=0.0) # only useful for the "pair" init (a pair of defects)

ratios_TKT = [0.2]
alpha_max = 2 * 2.3 / Lx # in my Frontiers in Physics, I've shown that αξ = 2.3 for t → ∞. If we want ξ ≈ L/2, then one obtains α = 2 * 2.3 / L
alphas = [0, alpha_max]
distribution_types = ["gaussian"]

thetas_saved_cpu = zeros(Float16, Lx, Ly, R, length(ratios_TKT), length(alphas), length(inits), length(distribution_types), length(times)) # in Float16 for storage reasons
magnetisations = zeros(R, length(ratios_TKT), length(alphas), length(inits), length(distribution_types), length(times))
number_of_defects = zeros(R, length(ratios_TKT), length(alphas), length(inits), length(distribution_types), length(times))
xs = Array{Vector{Float64}}(undef, R, length(ratios_TKT), length(alphas), length(inits), length(distribution_types), length(times)) # positions of the defects
ys = Array{Vector{Float64}}(undef, R, length(ratios_TKT), length(alphas), length(inits), length(distribution_types), length(times)) # positions of the defects
qs = Array{Vector{Float64}}(undef, R, length(ratios_TKT), length(alphas), length(inits), length(distribution_types), length(times)) # positions of the defects

m = 0
M = length(ratios_TKT) * length(inits) * length(alphas) * length(distribution_types)

z = @elapsed for i in each(ratios_TKT), j in each(alphas), k in each(inits), mm in each(distribution_types)

    alpha = Tf(alphas[j])
    T = Tf(round(ratios_TKT[i] *  0.22 * (1 - 0.43 * alpha * log(Lx)), digits=3))
        # in the XY model, TKT_xy = 0.89/4 = 0.22 is the critical temperature of the phase transition to disorder.
        # in the Kuramoto model, T_c = TKT_xy * (1 - 0.43 α log(L)) is the critical temperature below which the defects are bounded
    init = inits[k]
    distribution_type = distribution_types[mm]
    m += 1

    println("Simulation $m/$M : T = $T, α = $(alpha), Init = $init, ω ~ $distribution_type")

    thetas = create_thetas(Lx, Ly, R, init, params_init)
    thetas_new = similar(thetas) # allocates a similar array to thetas 
    omegas = create_omegas(Lx, Ly, R, alpha, distribution_type)

    t = Float64(0)
    for tt in 1:length(times)
        thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, times[tt], lattice_type)
        
        thetas_saved_cpu[:, :, :, i, j, k, mm, tt] = Array(thetas)
        # magnetisations[:, i, j, k, mm, tt] = OP(thetas)[1]
        # number_of_defects[:, i, j, k, mm, tt] = number_defects(thetas)

        # dico = spot_defects(thetas, lattice_type) # for a quick visual check if needed
        # qs[:, i, j, k, mm, tt] = dico["qs"]
        # xs[:, i, j, k, mm, tt] = dico["xs"]
        # ys[:, i, j, k, mm, tt] = dico["ys"]
    end
end
println("Simulation over.")
prinz(z) # prints the runtime in seconds, minutes and hours. Defined in `auxiliary.jl` file.

fig = plot_thetas(thetas_saved_cpu[:, :, :, 1, 1, 1, 1, end], defects=false) # for a quick visual check if needed 
fig
dico["xs"]
dico["number_defects"]

# @btime spot_defects(thetas_saved_cpu[:, :, :, 1, 1, 1, 1, end], lattice_type)

## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##

comments = "$(time_spacing)ly spaced times. $(lattice_type) lattice."

pwd()
filepath = pwd() * "/models/kuramoto/data/"
filename = "coarsening_Lx$(Lx)_Ly$(Ly)_R$(R)_ratios_TKT$(ratios_TKT)_alphas$(alphas)_inits_$(join(inits, "_"))_distributions_$(join(distribution_types, "_"))_tmax$(tmax)"

@save filepath * filename * ".jld2" Lx Ly R ratios_TKT alphas times tmax dt inits distribution_types comments runtime = z
# @save filepath * filename * "_with_thetas.jld2" thetas_saved_cpu Lx Ly R ratios_TKT alphas tmax times dt inits distribution_types comments runtime = z

println("Data saved in $(filepath * filename).jld2")


## ------------------ If you want to load back the data ------------------ ##
# filepath = pwd() * "/models/kuramoto/data/" # same as above
# filename = "coarsening_Lx256_Ly256_R4_Ts[0.05]_alphas[0.1]_inits_hightemp_distributions_uniform_gaussian_tmax10000.0"
# @load filepath * filename * ".jld2" Lx Ly R Ts alphas times tmax dt inits distribution_types comments runtime # without the "=z"