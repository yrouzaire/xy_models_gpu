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

Lx = 512
Ly = 512
R = 1 # number of independent realisations
Lmin = min(Lx, Ly)

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

inits = ["hightemp"] # choose among "lowtemp", "hightemp", "pair"
params_init = (x_plus=Lx / 4, y_plus=Ly / 2, r0=Lx / 2, mu_plus=0.0) # only useful for the "pair" init (a pair of defects)
Ts = [0.8] * 0.22 # in the XY model, T = 0.22 is the critical temperature of the phase transition to disorder. 
alphas = [0.] # try first between 0 and 0.2 
distribution_types = ["gaussian"] # choose among "uniform", "gaussian", "exponential", "cauchy" and "truncated_cauchy"

thetas_saved_cpu = zeros(Float16, Lx, Ly, R, length(Ts), length(alphas), length(inits), length(distribution_types), length(times)) # in Float16 for storage reasons
magnetisations = zeros(R, length(Ts), length(alphas), length(inits), length(distribution_types), length(times))
number_of_defects = zeros(R, length(Ts), length(alphas), length(inits), length(distribution_types), length(times))
Cs = zeros(round(Int,Lmin/2-1) , R, length(Ts), length(alphas), length(inits), length(distribution_types), length(times)) # correlation function C(r,t)
xis = zeros(R, length(Ts), length(alphas), length(inits), length(distribution_types), length(times)) # correlation length ξ(t)

mm = 0
M = length(Ts) * length(inits) * length(alphas) * length(distribution_types)

z = @elapsed for i in each(Ts), j in each(alphas), k in each(inits), m in each(distribution_types)

    T = Tf(Ts[i])
    alpha = Tf(alphas[j])
    init = inits[k]
    distribution_type = distribution_types[m]
    
    global mm += 1
    println("Simulation $mm/$M : T = $T, α = $(alpha), Init = $init, ω ~ $distribution_type")

    thetas = create_thetas(Lx, Ly, R, init, params_init)
    thetas_new = similar(thetas) # allocates a similar array to thetas 
    omegas = create_omegas(Lx, Ly, R, alpha, distribution_type)

    t = Float64(0)
    for tt in 1:length(times)
        thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, times[tt], lattice_type)
        
        thetas_saved_cpu[:, :, :, i, j, k, m, tt] = Array(thetas)

        magnetisations[:, i, j, k, m, tt] = Array(OP(thetas)[2])
        number_of_defects[:, i, j, k, m, tt] = Array(number_defects(thetas))
        tmp = Array(corr_fft(thetas))
        Cs[:, :, i, j, k, m, tt] = tmp
        xis[:, i, j, k, m, tt] = correlation_length(tmp)
    
    end
end
println("Simulation over.")
prinz(z) # prints the runtime in seconds, minutes and hours. Defined in `auxiliary.jl` file.


plot_thetas(thetas_saved_cpu[:, :, :, 1, 1, 1, 1, end], defects=false)
##

fig = Figure() 
for i in 1:R 
    ax_lowtemp = Axis(fig[1, i], aspect=DataAspect())
    # if i == 1 
    #     ax_lowtemp.ylabel = "Ordered init. config"
    # end
    hidedecorations!(ax_lowtemp)
    data = mod.(thetas_saved_cpu[:, :, i, 1, 1, 1, 1, end], 2pi)
    h = heatmap!(ax_lowtemp, data, 
        colormap=cols_thetas, 
        colorrange=(0, 2pi), )

    ax_hightemp = Axis(fig[2, i], aspect=DataAspect())
    # if i == 1 
    #     ax_hightemp.ylabel = "Random init. config"
    # end
    hidedecorations!(ax_hightemp)
    data = mod.(thetas_saved_cpu[:, :, i, 1, 1, 2, 1, end], 2pi)
    h = heatmap!(ax_hightemp, data, 
    colormap=cols_thetas, 
        colorrange=(0, 2pi), )
end
sideinfo_low = Label(fig[1, 0], "Initially Ordered", rotation = pi/2)
sideinfo_high = Label(fig[2, 0], "Initially Random", rotation = pi/2)
supertitle = Label(fig[0, :], L"L = 256, J = 1, T = 0.2, \Omega = 0, dt = 0.1, t_{max} = 1\cdot 10^5 ", fontsize = 30)
Colorbar(fig[1:2, R+1],label="θ",colormap=cols_thetas,
    ticklabelsize=20, labelsize=50, tickalign=1,colorrange=(0, 2pi), 
    ticks=([0, pi/2, pi, 3pi/2, 2pi], [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),)
for i in 1:R 
    colsize!(fig.layout, i, Aspect(1, Lx / Ly))
    colgap!(fig.layout, i, 6)
end
colgap!(fig.layout, R+1, 6)
rowgap!(fig.layout, 1, 6)
rowgap!(fig.layout, 2, 6)

resize_to_layout!(fig)
save(pwd()*"/snapshots_Li_tmax$(tmax).pdf", fig)
fig




## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##

comments = "$(time_spacing)ly spaced times. $(lattice_type) lattice."

pwd()
filepath = pwd() * "/models/kuramoto/data/"
filename = "coarsening_Lx$(Lx)_Ly$(Ly)_R$(R)_Ts$(Ts)_alphas$(alphas)_inits_$(join(inits, "_"))_distributions_$(join(distribution_types, "_"))_tmax$(tmax)"

@save filepath * filename * ".jld2" Lx Ly R Ts alphas times tmax dt inits distribution_types comments runtime = z
@save filepath * filename * "_with_thetas.jld2" thetas_saved_cpu Lx Ly R Ts alphas tmax times dt inits distribution_types comments runtime = z

println("Data saved in $(filepath * filename).jld2")


## ------------------ If you want to load back the data ------------------ ##
# filepath = pwd() * "/models/kuramoto/data/" # same as above
# filename = "coarsening_Lx256_Ly256_R4_Ts[0.05]_alphas[0.1]_inits_hightemp_distributions_uniform_gaussian_tmax10000.0"
# @load filepath * filename * ".jld2" Lx Ly R Ts alphas times tmax dt inits distribution_types comments runtime # without the "=z"
