# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include("../src/load_everything.jl")


## ------------------ Structure Factor ------------------ ##
## ------------------ Structure Factor ------------------ ##
## ------------------ Structure Factor ------------------ ##
## ------------------ Structure Factor ------------------ ##
# from ChatGPT
hist(fftfreq(128,1))
function compute_structure_factor(xs::Vector{Tf}, ys::Vector{Tf}, L, grid_size) where {Tf}
    N = length(xs)  # Number of particles
    dx = L / grid_size  # Grid spacing

    # Step 1: Discretize positions into a 2D density field
    rho = zeros(grid_size, grid_size)  # Density grid

    for i in 1:N
        ix = mod1(floor(Int, xs[i] / dx) + 1, grid_size)
        iy = mod1(floor(Int, ys[i] / dx) + 1, grid_size)
        rho[ix, iy] += 1  # Count particles in each grid cell
    end

    rho .-= mean(rho)  # Remove mean to focus on fluctuations

    # Step 2: Compute FFT of the density field
    rho_k = fft(rho)

    # Step 3: Compute power spectrum |rho_k|²
    S_k = abs.(rho_k) .^ 2 / N  # Normalize by particle number

    # Step 4: Average over circular bins to get S(q)
    qx = fftfreq(grid_size, 2 / L)  # Frequency grid (qx)
    qy = fftfreq(grid_size, 2 / L)  # Frequency grid (qy)
    q_magnitudes = sqrt.(qx' .^ 2 .+ qy .^ 2)  # Compute |q|

    q_bins = range(0, maximum(q_magnitudes), length=grid_size ÷ 2)
    S_q = zeros(length(q_bins) - 1)
    counts = zeros(length(q_bins) - 1)

    for i in eachindex(q_magnitudes)
        q_val = q_magnitudes[i]
        bin_idx = searchsortedfirst(q_bins, q_val) - 1
        if bin_idx > 0
            S_q[bin_idx] += S_k[i]
            counts[bin_idx] += 1
        end
    end

    S_q ./= counts  # Normalize
    q_centers = (q_bins[1:end-1] .+ q_bins[2:end]) / 2  # Midpoints of bins

    return q_centers, S_q
end
function compute_structure_factor(xs::Vector{Vector{Tf}}, ys::Vector{Vector{Tf}}, L, grid_size) where {Tf}
    R = length(xs)
    S_q_vals_all = zeros(round(Int, grid_size / 2 - 1), R)
    for r in 1:R
        q_vals_all, S_q_vals_all[:, r] = compute_structure_factor(xs[r], ys[r], L, grid_size)
    end
    q_vals = compute_structure_factor(xs[1], ys[1], L, grid_size)[1]
    return q_vals, S_q_vals_all
end


## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##

pwd()
filepath = pwd() * "/models/nrkuramoto/data/"
filename = "test_hyperuniformity_Lx512_Ly512_Rtot1024_ratios_TKT[0.2]_inits_hightemp_distributions_gaussian_tmax10000.0"
filename = "test_hyperuniformity_Lx1024_Ly1024_Rtot4_ratios_TKT[0.3]_inits_lowtemp_distributions_gaussian_tmax30000.0_with_thetas_final_time"
filename = "test_hyperuniformity_Lx1024_Ly1024_Rtot64_ratios_TKT[0.3]_inits_lowtemp_hightemp_distributions_gaussian_tmax50000.0_with_thetas_final_time"
filename = "test_hyperuniformity_Lx512_Ly512_Rtot256_ratios_TKT[0.5, 0.8]_inits_hightemp_distributions_gaussian_tmax10000.0"
@load filepath * filename * ".jld2"  number_of_defects qs xs ys Lx Ly R R_simus R_tot ratios_TKT alphas_sigmas times tmax dt inits distribution_types comments runtime
# thetas_saved_cpu_final_time

#= 
Size of qs, xs, ys 
R_tot, ratios_TKT, alphas_sigmas, inits, distribution_types, times 
=#



## ------------------ Number of defects ------------------ ##
## ------------------ Number of defects ------------------ ##
## ------------------ Number of defects ------------------ ##
## ------------------ Number of defects ------------------ ##

fig = Figure(size=(1000, 500))
ax = Axis(fig[1, 1], xlabel = L"t", ylabel =L"n", 
    xscale = log10, yscale = log10) 
data = mean(number_of_defects, dims=1)
for (ind_alfsig, value) in enumerate(alphas_sigmas)
    labb = "α = $(round(value[1], digits=2)), σ = $(value[2])"
    scatterlines!(ax, times, data[1, ind_T, ind_alfsig, 1, 1, :], label=labb)
    # scatterlines!(ax, times, data[1, ind_T, ind_alfsig, 2, 1, :], label=labb)
end
vlines!(ax, 1e4, color=:grey80)
lines!(ax, times[5:end], x -> 3E4 * log(x) / x, color=:black)
# axislegend()


ax2 = Axis(fig[1, 2], xlabel=L"n", ylabel="PDF(n)", 
    #xscale = log10, yscale = log10, 
    # xticks= 0:10:80
    )
for i in each(alphas_sigmas)
    data = number_of_defects[:, 1, i, 1, 1, end]
    hist!(ax2, data, bins=15, normalization=:pdf)
    vlines!(ax2, [mean(data)], color=:black)
end
# text!(ax2, 1, 1, text="⟨n⟩ = "*string(round(mean(data), digits=2)), space=:relative, offset=(-8, -4), align=(:right, :top))
# text!(ax, 1, 1, text="(a)",
#     fontsize=24, font=:bold_italic,
#     align=(:right, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)

fig


## ------------------ Visualise field theta ------------------ ##
## ------------------ Visualise field theta ------------------ ##
## ------------------ Visualise field theta ------------------ ##
## ------------------ Visualise field theta ------------------ ##
data = thetas_saved_cpu_final_time[:, :, 1, 1, 1, 1, 1]
plot_thetas(data, defects=false)

## ------------------ Analysis ------------------ ##
## ------------------ Analysis ------------------ ##
## ------------------ Analysis ------------------ ##
## ------------------ Analysis ------------------ ##

ind_T = 2
ind_alfsig = 2
dr = 1
ind_time = length(times) 
indices_R = 1:1
indices_R = 1:R_tot
last_n = 1 # one can merge the last_n time defect config to obtain better stats 
grid_size = round(Int, Lx/2)  # Resolution of density grid


qss = qs[indices_R, ind_T, ind_alfsig, 1, 1, ind_time]
xss = xs[indices_R, ind_T, ind_alfsig, 1, 1, ind_time]
yss = ys[indices_R, ind_T, ind_alfsig, 1, 1, ind_time]

# qss = vcat([qs[indices_R, ind_T, ind_alfsig, 1, 1, ttt] for ttt in length(times)-last_n+1:length(times)]...)
# xss = vcat([xs[indices_R, ind_T, ind_alfsig, 1, 1, ttt] for ttt in length(times)-last_n+1:length(times)]...)
# yss = vcat([ys[indices_R, ind_T, ind_alfsig, 1, 1, ttt] for ttt in length(times)-last_n+1:length(times)]...)

# # # # # Positions # # # # #
xs_plus, ys_plus, xs_minus, ys_minus = split_plus_minus(qss, xss, yss)

NN = 50
# NN = round(Int, sum(Int, number_of_defects[:, 1, 1, 1, 1, end])/R) # to have approx the same number of defects in total compared to the real defect number 
# xs_plus, ys_plus, xs_minus, ys_minus = [Lx * rand(NN) for r in 1:R], [Ly * rand(NN) for r in 1:R], [Lx * rand(NN) for r in 1:R], [Ly * rand(NN) for r in 1:R]
# using Sobol
# s = SobolSeq(2)
# p = Lx * reduce(hcat, next!(s) for i = 1:2NN)'
# xs_plus, ys_plus, xs_minus, ys_minus = [p[1:NN, 1] for r in 1:R], [p[1:NN, 2] for r in 1:R], [p[NN+1:end, 1] for r in 1:R], [p[NN+1:end, 2] for r in 1:R]

fig = Figure() 
ax = Axis(fig[1, 1], aspect=DataAspect(),
    limits = (0, Lx, 0, Ly)
    )
for r in 1:length(indices_R)
    scatter!(ax, xs_plus[r], ys_plus[r], color=:blue)
    scatter!(ax, xs_minus[r], ys_minus[r], color=:red)
end
fig
# #


# # # # # g(r) and S(q) # # # #
rs = collect(1:dr:Lx/2)
ks = reverse(2pi ./ collect(1:dr:Lx))


dico_GR = gr(xs_plus, ys_plus, xs_minus, ys_minus, Lx, length(indices_R), dr=dr)
GR_all_avg = dico_GR["all_avg"]
GR_plus_minus_avg = dico_GR["plus_minus_avg"]
GR_plus_plus_avg = dico_GR["plus_plus_avg"]
GR_minus_minus_avg = dico_GR["minus_minus_avg"]

q_vals_minus_minus, S_q_vals_minus_minus = compute_structure_factor(xs_minus, ys_minus, Lx, grid_size)
q_vals_plus_plus, S_q_vals_plus_plus = compute_structure_factor(xs_plus, ys_plus, Lx, grid_size)
# q_vals_plus_minus, S_q_vals_plus_minus = compute_structure_factor(xs_plus, ys_plus, Lx, grid_size)
q_vals_all, S_q_vals_all = compute_structure_factor(vcat(xs_plus, xs_minus), vcat(ys_plus, ys_minus), Lx, grid_size)
S_q_vals_avg = mean(S_q_vals, dims=2)[:,1]
# S_q_vals_plus_minus_avg = mean(S_q_vals_plus_minus, dims=2)[:,1]
S_q_vals_plus_plus_avg = mean(S_q_vals_plus_plus, dims=2)[:,1]
S_q_vals_minus_minus_avg = mean(S_q_vals_minus_minus, dims=2)[:,1]


# dico_SF = structure_factor(xs_plus, ys_plus, xs_minus, ys_minus, Lx, length(indices_R), dr=dr)
# SF_all_avg = dico_SF["all_avg"]
# SF_plus_minus_avg = dico_SF["plus_minus_avg"]
# SF_plus_plus_avg = dico_SF["plus_plus_avg"]
# SF_minus_minus_avg = dico_SF["minus_minus_avg"]

# dico_SF_gr = structure_factor_via_gr(xs_plus, ys_plus, xs_minus, ys_minus, Lx, length(indices_R), dr=dr)
# SF_all_gr_avg = dico_SF_gr["all_avg"]
# SF_plus_minus_gr_avg = dico_SF_gr["plus_minus_avg"]
# SF_plus_plus_gr_avg = dico_SF_gr["plus_plus_avg"]
# SF_minus_minus_gr_avg = dico_SF_gr["minus_minus_avg"]


# # # # # Figure # # # # #
fig = Figure(size=(500, 1200))

# Pair correlation function g(r)
ax = Axis(
    fig[1, 1],
    xlabel=L"r",
    ylabel=L"g(r)",
    title="Pair correlation function",
    # limits=(nothing, Lx/2, 0, 4.5),  yticks=0:0.5:4, 
    xscale=log10, yscale=log10
)
lines!(ax, rs, GR_all_avg, label="All")
lines!(ax, rs, GR_plus_minus_avg, label="+/-")
lines!(ax, rs, GR_plus_plus_avg, label="+/+")
lines!(ax, rs, GR_minus_minus_avg, label="-/-")
hlines!(ax, [1], color=:grey80)
axislegend(ax)
fig


# Pair correlation function S(q)
ax = Axis(fig[2, 1], xlabel=L"k", ylabel=L"S(k)", title="Structure factor",
limits = (nothing, nothing, nothing, nothing),
    xscale=log10)#, yscale = log10)
# lines!(ax, q_vals_all, S_q_vals_all_avg, label="All")
scatterlines!(ax, q_vals_plus_plus, S_q_vals_plus_plus_avg, label="+/+")
scatterlines!(ax, q_vals_minus_minus, S_q_vals_minus_minus_avg, label="-/-")
# lines!(ax, q_vals_plus_minus, S_q_vals_plus_minus_avg, label="+/-")
# lines!(ax, q_vals_plus_plus[1:5], x->2E2x^0.5,color=:black)
axislegend(ax, position=:rb)




# lines!(ax, ks, SF_all_avg, label="All")
# lines!(ax, ks, SF_plus_minus_avg, label="+/-")
# lines!(ax, ks, SF_plus_plus_avg, label="+/+")
# lines!(ax, ks, SF_minus_minus_avg, label="-/-")
# axislegend(ax, position=:rt)


# # Pair correlation function S(q) from g(r)
# ax = Axis(fig[3, 1], xlabel=L"k", ylabel=L"S(k)", title="Structure factor via g(r)",
#     xscale=log10) #, yscale = log10)

# lines!(ax, ks[round(Int, end / 2 + 1):end], SF_all_gr_avg, label="All")
# lines!(ax, ks[round(Int, end / 2 + 1):end], SF_plus_minus_gr_avg, label="+/-")
# lines!(ax, ks[round(Int, end / 2 + 1):end], SF_plus_plus_gr_avg, label="+/+")
# lines!(ax, ks[round(Int, end / 2 + 1):end], SF_minus_minus_gr_avg, label="-/-")
# axislegend(ax, position=:lt)


resize_to_layout!(fig)
fig

##


