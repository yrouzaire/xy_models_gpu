# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include("src/load_everything.jl")


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
duration_movies_in_seconds = 12
frame_per_seconds = 30
nb_frames = frame_per_seconds * duration_movies_in_seconds
every = tmax / nb_frames
time_spacing = "linear" # "linear", "log", "quadratic"
if time_spacing == "linear"
    times = Tf.(collect(every:every:tmax)) # linear
elseif time_spacing == "log"
    times = logspace(10dt, tmax, nb_frames) # log
elseif time_spacing == "quadratic"
    times = Tf.(range(sqrt(dt), sqrt(tmax), nb_frames)) .^ 2  # quadratic
else
    error("time_spacing should be 'linear', 'log' or 'quadratic'")
end


inits = ["lowtemp", "hightemp"]
params_init = (x_plus=Lx / 4, y_plus=Ly / 2, r0=Lx / 2, mu_plus=0.0) # only useful for the "pair" init (a pair of defects) 
Ts = [0.05, 0.1]
sigmas = [0, 0.1]
distribution_types = ["uniform", "gaussian"]

thetas_saved_cpu = zeros(Float16, Lx, Ly, R, length(Ts), length(sigmas), length(inits), length(distribution_types), length(times))
thetas_saved_cpu
m = 0
M = length(Ts) * length(inits) * length(sigmas) * length(distribution_types)

z = @elapsed for i in each(Ts), j in each(sigmas), k in each(inits), mm in each(distribution_types)

    T = Tf(Ts[i])
    sigma = Tf(sigmas[j])
    init = inits[k]
    distribution_type = distribution_types[mm]
    m += 1

    println("Simulation $m/$M : T = $T, σ = $(sigma), Init = $init, ω ~ $distribution_type")

    thetas = create_thetas(Lx, Ly, R, init, params_init)
    thetas_new = similar(thetas)
    omegas = create_omegas(Lx, Ly, R, sigma, distribution_type)

    t = Float64(0)
    for tt in each(times)
        # println("t = $(round(t, digits=2)), $(round(100t/tmax, digits=2)) %")
        thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, times[tt], lattice_type)
        thetas_saved_cpu[:, :, :, i, j, k, mm, tt] = Array(thetas)
    end
end
prinz(z)


## ------------------------ Save data ------------------------ ##
## ------------------------ Save data ------------------------ ##
## ------------------------ Save data ------------------------ ##
## ------------------------ Save data ------------------------ ##
pwd()

comments = "$(time_spacing)ly spaced times. $(lattice_type) lattice. "
filepath = pwd() * "/models/kuramoto/movies/data_for_movies/"
filename = "Lx$(Lx)_Ly$(Ly)_R$(R)_Ts$(Ts)_sigmas$(sigmas)_inits_$(join(inits, "_"))_distributions_$(join(distribution_types, "_"))_tmax$(tmax)_.jld2"
@save filepath * filename thetas_saved_cpu times Lx Ly R Ts sigmas tmax times dt inits distribution_types comments runtime = z

## ------------------------ Load data ------------------------ ##
## ------------------------ Load data ------------------------ ##
## ------------------------ Load data ------------------------ ##
## ------------------------ Load data ------------------------ ##
pwd()

filepath = pwd() * "/models/kuramoto/movies/data_for_movies/"
filename = ".jld2"
@load filepath * filename thetas_saved_cpu times Lx Ly R Ts sigmas tmax times dt inits distribution_types comments runtime


## ------------------------ Make the movie ------------------------ ##
## ------------------------ Make the movie ------------------------ ##
## ------------------------ Make the movie ------------------------ ##
## ------------------------ Make the movie ------------------------ ##
using GLMakie
GLMakie.activate!()

zm = @elapsed for ind_T in each(Ts), ind_init in each(inits), ind_sig in each(sigmas), ind_distrib in each(distribution_types), rr in 1:R

    T = Ts[ind_T]
    sigma = sigmas[ind_sig]
    init = inits[ind_init]
    distribution_type = distribution_types[ind_distrib]

    fig = Figure()
    ax1 = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y",
        ylabelrotation=0, xtickalign=0, ytickalign=0,
        title=L"t = 0, σ = %$(sigma), T = %$(T)")
    hidedecorations!(ax1)

    data_to_plot = Observable(mod.(thetas_saved_cpu[:, :, rr, ind_T, ind_sig, ind_init, ind_distrib, 1], Float32(2pi)))
    h = CairoMakie.heatmap!(ax1, data_to_plot,
        colormap=cols_thetas,
        colorrange=(0, 2pi))


    Colorbar(fig[1, 2], h,
        label=L"θ",
        ticklabelsize=20, labelrotation=0, labelsize=40,
        tickalign=1,
        labelpadding=-10, # distance between the colorbar ticks and the label
        ticks=([0, pi / 2, pi, 3pi / 2, 2pi], [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),
    )
    colsize!(fig.layout, 1, Aspect(1, Lx / Ly))
    resize_to_layout!(fig)
    fig


    # nframes = 10 # to test on small number of frames if needed
    nframes = length(times)
    filename = "/$(lowercase(distribution_type))_Lx$(Lx)_Ly$(Ly)_T$(T)_sigma$(sigma)_tmax$(tmax)_r$(rr)"

    GLMakie.record(fig, filepath * "/$(lowercase(init))/" * filename * ".mp4", 1:nframes, fps=frame_per_seconds) do tt
        println("Frame $(round(100tt / nframes,digits=2)) %")
        data_to_plot = Observable(mod.(thetas_saved_cpu[:, :, rr, ind_T, ind_sig, ind_init, ind_distrib, tt], Float32(2pi)))
        ax1.title = L"t = %$(round(times[tt], digits=1)), σ = %$(sigma), T = %$(T)"
    end

end
prinz(zm)
