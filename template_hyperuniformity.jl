# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include("src/load_everything.jl")

## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##

lattice_type = "square"
lattice_type = "triangular"

Lx = 512
Ly = 512
@assert Lx == Ly "For this script, Lx must be equal to Ly"
R = 16 # number of independent realisations

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
params_init = (x_plus=Lx / 4, y_plus=Ly / 2, r0=Lx / 2, mu_plus=0.0) # only useful for the "pair" init (a pair of defects)

ratios_TKT = [0.2]
alpha_L2 = 2 * 2.3 / Lx # in my Frontiers in Physics, I've shown that αξ = 2.3 for t → ∞. If we want ξ ≈ L/2, then one obtains α = 2 * 2.3 / L
alphas = [alpha_L2]
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
    for tt in eachindex(times)
        thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, times[tt], lattice_type)
        
        thetas_saved_cpu[:, :, :, i, j, k, mm, tt] = Array(thetas)
        # magnetisations[:, i, j, k, mm, tt] = OP(thetas)[1]
        
        dico = spot_defects(thetas) # for a quick visual check if needed
        qs[:, i, j, k, mm, tt] = dico["qs"]
        xs[:, i, j, k, mm, tt] = dico["xs"]
        ys[:, i, j, k, mm, tt] = dico["ys"]
        number_of_defects[:, i, j, k, mm, tt] = dico["number_defects"]
    end
end
println("Simulation over.")
prinz(z) # prints the runtime in seconds, minutes and hours. Defined in `auxiliary.jl` file
titre = "Total n = "*string(round(sum(Int, number_of_defects[:, 1, 1, 1, 1, end])))
plot_thetas(thetas_saved_cpu[:, :, 1, 1, 1, 1, 1, end], defects=true, title=titre) # for a quick visual check if needed 

## Structure Factor S(q)

xs_plus, ys_plus, xs_minus, ys_minus = split_plus_minus(qs[:, 1, 1, 1, 1, end], xs[:, 1, 1, 1, 1, end], ys[:, 1, 1, 1, 1, end])
NN = 50
# NN = round(Int, sum(Int, number_of_defects[:, 1, 1, 1, 1, end])/R) # to have approx the same number of defects in total compared to the real defect number 
# xs_plus, ys_plus, xs_minus, ys_minus = [Lx * rand(NN) for r in 1:R], [Ly * rand(NN) for r in 1:R], [Lx * rand(NN) for r in 1:R], [Ly * rand(NN) for r in 1:R]
using Sobol
# s = SobolSeq(2)
# p = Lx * reduce(hcat, next!(s) for i = 1:2NN)'
# xs_plus, ys_plus, xs_minus, ys_minus = [p[1:NN, 1] for r in 1:R], [p[1:NN, 2] for r in 1:R], [p[NN+1:end, 1] for r in 1:R], [p[NN+1:end, 2] for r in 1:R]


fig = Figure(size=(500, 500))
ax = Axis(fig[1, 1], title="Defects positions", limits=(0, Lx, 0, Ly))
scatter!(ax, xs_plus[1], ys_plus[1], color=:red)
scatter!(ax, xs_minus[1], ys_minus[1], color=:green)
scatter!(ax, xs_plus[1], ys_minus[1], color=:blue)
scatter!(ax, xs_minus[1], ys_plus[1], color=:gold)
resize_to_layout!(fig)
fig



# #

dr = 3
ks = reverse(2pi ./ collect(1:dr:Lx))



dico_SF = structure_factor(xs_plus, ys_plus, xs_minus, ys_minus, Lx, R, dr = dr)
SF_all_avg = dico_SF["all_avg"]
SF_plus_minus_avg = dico_SF["plus_minus_avg"]
SF_plus_plus_avg = dico_SF["plus_plus_avg"]
SF_minus_minus_avg = dico_SF["minus_minus_avg"]

fig = Figure(size=(500, 1200))
ax = Axis(fig[1, 1], xlabel=L"k", ylabel=L"S(k)", title="Structure factor",
    xscale=log10) #, yscale = log10)
lines!(ax, ks, SF_all_avg, label="All")
lines!(ax, ks, SF_plus_minus_avg, label="+/-")
lines!(ax, ks, SF_plus_plus_avg, label="+/+")
lines!(ax, ks, SF_minus_minus_avg, label="-/-")
axislegend(ax, position=:rt)


ax = Axis(fig[2,1], xlabel=L"k", ylabel=L"S(k)", title="Structure factor via g(r)",
    xscale=log10) #, yscale = log10)
dico_SF_gr = structure_factor_via_gr(xs_plus, ys_plus, xs_minus, ys_minus, Lx, R, dr=dr)
SF_all_avg = dico_SF_gr["all_avg"]
SF_plus_minus_avg = dico_SF_gr["plus_minus_avg"]
SF_plus_plus_avg = dico_SF_gr["plus_plus_avg"]
SF_minus_minus_avg = dico_SF_gr["minus_minus_avg"]
lines!(ax, ks[round(Int,end/2 + 1):end], SF_all_avg, label="All")
lines!(ax, ks[round(Int,end/2 + 1):end], SF_plus_minus_avg, label="+/-")
lines!(ax, ks[round(Int,end/2 + 1):end], SF_plus_plus_avg, label="+/+")
lines!(ax, ks[round(Int,end/2 + 1):end], SF_minus_minus_avg, label="-/-")
axislegend(ax, position=:lt)


# # Pair correlation function g(r)
rs = collect(1:dr:Lx/2)

dico_GR = gr(xs_plus, ys_plus, xs_minus, ys_minus, Lx, R, dr = dr)
GR_all_avg = dico_GR["all_avg"]
GR_plus_minus_avg = dico_GR["plus_minus_avg"]
GR_plus_plus_avg = dico_GR["plus_plus_avg"]
GR_minus_minus_avg = dico_GR["minus_minus_avg"]

# fig = Figure()
ax = Axis(fig[3, 1], xlabel=L"r", ylabel=L"g(r)", title="Pair correlation function", 
    # limits=(0, Lx/2, 0, 2.3),  yticks=0:0.5:2
    )
lines!(ax, rs, GR_all_avg, label="All")
lines!(ax, rs, GR_plus_minus_avg, label="+/-")
lines!(ax, rs, GR_plus_plus_avg, label="+/+")
lines!(ax, rs, GR_minus_minus_avg, label="-/-")
hlines!(ax, [1], color=:grey80)
axislegend(ax)
fig














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
