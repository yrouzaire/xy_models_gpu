# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include("../src/load_everything.jl")

## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##

pwd()
filepath = pwd() * "/models/nrkuramoto/data/"
filename = "test_hyperuniformity_Lx512_Ly512_Rtot1024_ratios_TKT[0.2]_inits_hightemp_distributions_gaussian_tmax10000.0"

@load filepath * filename * ".jld2" Lx Ly R R_simus R_tot ratios_TKT alphas_sigmas times tmax dt inits distribution_types comments runtime





## ------------------ Analysis ------------------ ##
## ------------------ Analysis ------------------ ##
## ------------------ Analysis ------------------ ##
## ------------------ Analysis ------------------ ##



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



dico_SF = structure_factor(xs_plus, ys_plus, xs_minus, ys_minus, Lx, R, dr=dr)
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


ax = Axis(fig[2, 1], xlabel=L"k", ylabel=L"S(k)", title="Structure factor via g(r)",
    xscale=log10) #, yscale = log10)
dico_SF_gr = structure_factor_via_gr(xs_plus, ys_plus, xs_minus, ys_minus, Lx, R, dr=dr)
SF_all_avg = dico_SF_gr["all_avg"]
SF_plus_minus_avg = dico_SF_gr["plus_minus_avg"]
SF_plus_plus_avg = dico_SF_gr["plus_plus_avg"]
SF_minus_minus_avg = dico_SF_gr["minus_minus_avg"]
lines!(ax, ks[round(Int, end / 2 + 1):end], SF_all_avg, label="All")
lines!(ax, ks[round(Int, end / 2 + 1):end], SF_plus_minus_avg, label="+/-")
lines!(ax, ks[round(Int, end / 2 + 1):end], SF_plus_plus_avg, label="+/+")
lines!(ax, ks[round(Int, end / 2 + 1):end], SF_minus_minus_avg, label="-/-")
axislegend(ax, position=:lt)


# # Pair correlation function g(r)
rs = collect(1:dr:Lx/2)

dico_GR = gr(xs_plus, ys_plus, xs_minus, ys_minus, Lx, R, dr=dr)
GR_all_avg = dico_GR["all_avg"]
GR_plus_minus_avg = dico_GR["plus_minus_avg"]
GR_plus_plus_avg = dico_GR["plus_plus_avg"]
GR_minus_minus_avg = dico_GR["minus_minus_avg"]

# fig = Figure()
ax = Axis(
    fig[3, 1],
    xlabel=L"r",
    ylabel=L"g(r)",
    title="Pair correlation function",
    # limits=(0, Lx/2, 0, 2.3),  yticks=0:0.5:2
)
lines!(ax, rs, GR_all_avg, label="All")
lines!(ax, rs, GR_plus_minus_avg, label="+/-")
lines!(ax, rs, GR_plus_plus_avg, label="+/+")
lines!(ax, rs, GR_minus_minus_avg, label="-/-")
hlines!(ax, [1], color=:grey80)
axislegend(ax)
fig













