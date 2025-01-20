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

tmax = Tf(1E3)
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
    times = Tf.(range(sqrt(dt), sqrt(tmax), nb_frames)) .^2  # quadratic
else 
    error("time_spacing should be 'linear', 'log' or 'quadratic'")
end


inits = ["lowtemp", "hightemp"]
Ts = [0.05, 0.1]

thetas_saved_cpu = zeros(Float16, Lx, Ly, R, length(Ts), length(inits), length(times))

m = 0
M = length(Ts) * length(inits) 

z = @elapsed for i in each(Ts), j in each(inits)
    
    T = Tf(Ts[i])
    init = inits[j]
    m += 1

    println("Simulation $m/$M : T = $T, Init = $init")

    thetas = create_thetas(Lx, Ly, R, init)
    thetas_new = similar(thetas)
    
    t = Float64(0)
    for tt in each(times)
        # println("t = $(round(t, digits=2)), $(round(100t/tmax, digits=2)) %")
        thetas, t = evolve_XY!(thetas, thetas_new, Lx, Ly, R, T, t, dt, times[tt], lattice_type)
        thetas_saved_cpu[:, :, :, i, j, tt] = Array(thetas)
    end
end
prinz(z)


## ------------------------ Save data ------------------------ ##
## ------------------------ Save data ------------------------ ##
## ------------------------ Save data ------------------------ ##
## ------------------------ Save data ------------------------ ##
pwd()

comments = "$(time_spacing)ly spaced times. $(lattice_type) lattice. "
filepath = pwd()*"/models/xy/movies/data_for_movies/"
filename = "Lx$(Lx)_Ly$(Ly)_R$(R)_Ts$(Ts)_inits_$(join(inits, "_"))_tmax$(tmax)_.jld2"
@save filepath*filename thetas_saved_cpu times Lx Ly R Ts tmax times dt inits comments runtime=z 

## ------------------------ Load data ------------------------ ##
## ------------------------ Load data ------------------------ ##
## ------------------------ Load data ------------------------ ##
## ------------------------ Load data ------------------------ ##
pwd() 

filepath = pwd() * "/models/xy/movies/data_for_movies/"
filename = "Lx64_Ly64_R4_Ts[0.05]_inits_lowtemp_hightemp_tmax40.0.jld2"
@load filepath * filename thetas_saved_cpu times Lx Ly R Ts tmax times dt inits comments runtime


## ------------------------ Make the movie ------------------------ ##
## ------------------------ Make the movie ------------------------ ##
## ------------------------ Make the movie ------------------------ ##
## ------------------------ Make the movie ------------------------ ##
using GLMakie
GLMakie.activate!()

zm = @elapsed for ind_T in 1:length(Ts), ind_init in 1:length(inits), rr in 1:R

    T = Ts[ind_T]
    init = inits[ind_init]

    fig = Figure()
    ax1 = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y",
        ylabelrotation=0, xtickalign=0, ytickalign=0,
        title=L"t = 0, T = %$(T)")
    hidedecorations!(ax1)

    data_to_plot = Observable(mod.(thetas_saved_cpu[:, :, rr, ind_T, ind_init, 1], Float32(2pi)))
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
    filename = "Lx$(Lx)_Ly$(Ly)_T$(T)_tmax$(tmax)_r$(rr)"

    GLMakie.record(fig, filepath * "/$(lowercase(init))/" * filename * ".mp4", 1:nframes, fps=frame_per_seconds) do tt
        println("Frame $(round(100tt / nframes,digits=2)) %")
        data_to_plot[] = mod.(thetas_saved_cpu[:, :, rr, ind_T, ind_init, tt], Float32(2pi))
        L"t = %$(round(times[tt], digits=1)), T = %$(T)"
    end

end
prinz(zm)
