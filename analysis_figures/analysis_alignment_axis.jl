# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include(pwd()*"/src/load_everything.jl")

function plot_circle!(fig, ax; radius=1)
    θ = range(0, 2π, length=100)
    x_circle = radius*cos.(θ)
    y_circle = radius*sin.(θ)
    lines!(ax, x_circle, y_circle, color=:black, linewidth=2)
    return fig
end
function plot_axes!(fig, ax, lattice_type_string)
    if lattice_type_string == "square"
        hlines!(ax, [0], color=:grey, linewidth=1)
        vlines!(ax, [0], color=:grey, linewidth=1)
    elseif lattice_type_string == "triangular"
        lines!(ax, [-2, +2], sqrt(3) * [-2, +2], color=:grey, linewidth=1)
        lines!(ax, [-2, +2], sqrt(3) * [+2, -2], color=:grey, linewidth=1)
        hlines!(ax, [0], color=:grey, linewidth=1)

    end

    return fig  
end

## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##
## ------------------ Load data ------------------ ##

pwd()
filepath = pwd() * "/models/nrxy/data/"
filename = "alignment_axis_Lx128_Ly128_Rtot384_ratios_TKT[0.1]_tmax100000.0"
filename = "alignment_axis_Lx128_Ly128_Rtot32_ratios_TKT[0.1]_tmax100000.0"
filename = "comparison_soft_sharp_Lx128_Ly128_Rtot480_ratios_TKT[0.2]_tmax1.0e6"
@load filepath * filename * ".jld2" number_of_defects magnetisations Cs xis Lx Ly R R_simus R_tot ratios_TKT sigmas times tmax dt init comments runtime 
# @save filepath * filename * "_with_thetas_final_time.jld2" number_of_defects magnetisations Cs xis thetas_saved_cpu Lx Ly R R_simus R_tot ratios_TKT sigmas tmax times dt init comments runtime = z
# lattice_type = "square"
sharp_vision_cones = Tf.(2pi .- 16 * pi * sigmas ./ (8sigmas .+ pi * (4 .+ sigmas .^ 2))) # equivalence between sharp and smooth vision cones

##

filename = "triangular_dynamics_comparison_soft_sharp_Lx128_Ly128_Rtot64_ratios_TKT[0.2]_tmax100000.0"
filename = "smalldt_square_dynamics_comparison_soft_sharp_Lx128_Ly128_Rtot128_ratios_TKT[0.2]_tmax40000.0"
@load filepath * filename * ".jld2" magnetisations Lx Ly R R_simus R_tot ratios_TKT sigmas times tmax dt init comments runtime lattice_type
sharp_vision_cones = Tf.(2pi .- 16 * pi * sigmas ./ (8sigmas .+ pi * (4 .+ sigmas .^ 2))) # equivalence between sharp and smooth vision cones
number_of_defects_avg = mean(number_of_defects, dims=2)
xis_avg = mean(xis, dims=2)
prinz(runtime)



## ------------------ Scatter Plots comparison ------------------ ##
## ------------------ Scatter Plots comparison ------------------ ##   
## ------------------ Scatter Plots comparison ------------------ ##
## ------------------ Scatter Plots comparison ------------------ ##
using CairoMakie
CairoMakie.activate!()
markers = [:circle, :cross]

ind_sig = 15
ind_sig = length(sigmas)

filename = "comparison_soft_sharp_Lx128_Ly128_Rtot480_ratios_TKT[0.2]_tmax1.0e6"
@load pwd() * "/models/nrxy/data/" * filename * ".jld2" number_of_defects magnetisations Cs xis Lx Ly R R_simus R_tot ratios_TKT sigmas times tmax dt init comments runtime
angles_square = magnetisations[:, 1, :, :, :, end]
magnitudes_square = magnetisations[:, 2, :, :, :, end]


filename = "triangular_comparison_soft_sharp_Lx128_Ly128_Rtot480_ratios_TKT[0.2]_tmax1.0e6"
# @load pwd() * "/models/nrxy/data/" * filename * ".jld2" number_of_defects magnetisations Cs xis Lx Ly R R_simus R_tot ratios_TKT sigmas times tmax dt init comments runtime
# angles_triangular = magnetisations[:, 1, :, :, :, end]
# magnitudes_triangular = magnetisations[:, 2, :, :, :, end]

ind_sig = 2
fig = Figure(size=(400, 400)) 
ax = Axis(fig[1, 1]; limits=(-1., 1., -1., 1.), 
    xlabel=L"m_x", ylabel=L"m_y", 
        aspect=1, xticks=-1:1, yticks=-1:1, 
    title=L"σ = %$(round(sigmas[ind_sig], digits=2)) \ ,\ Θ = %$(round(sharp_vision_cones[ind_sig], digits=2))")
for j in 1#:length(ratios_TKT)
    for k in ind_sig#:length(sigmas)
        for tt in length(times)-3
            for i in 1:2 # soft / sharp
                mag_x = magnetisations[i, 2, :, j, k, tt] .* cos.(magnetisations[i, 1, :, j, k, tt])
                mag_y = magnetisations[i, 2, :, j, k, tt] .* sin.(magnetisations[i, 1, :, j, k, tt])
                if i == 1 
                    scatter!(ax, mag_x, mag_y, markersize = 7, marker = :circle, color = :green, label="Smooth kernel")
                elseif i == 2
                    scatter!(ax, mag_x, mag_y, markersize = 10, marker = :cross, color = :tomato, label="Sharp kernel")
                end
            end
        end
    end
end
plot_circle!(fig, ax, radius=1)


ind_sig = length(sigmas)
ax2 = Axis(fig[1,2]; limits=(-1.0, 1.0, -1.0, 1.0),
    xlabel=L"m_x", ylabel=L"m_y",
    aspect=1, xticks=-1:1, yticks=-1:1,
    title=L"σ = %$(round(sigmas[ind_sig], digits=2)) \ ,\ Θ = %$(round(sharp_vision_cones[ind_sig], digits=2))")
for j in 1#:length(ratios_TKT)
    for k in ind_sig#:length(sigmas)
        for tt in length(times) - 3
            for i in 1:2 # soft / sharp
                mag_x = magnetisations[i, 2, :, j, k, tt] .* cos.(magnetisations[i, 1, :, j, k, tt])
                mag_y = magnetisations[i, 2, :, j, k, tt] .* sin.(magnetisations[i, 1, :, j, k, tt])
                if i == 1
                    scatter!(ax2, mag_x, mag_y, markersize=7, marker=:circle, color=:green, label="Smooth kernel")
                elseif i == 2
                    scatter!(ax2, mag_x, mag_y, markersize=10, marker=:cross, color=:tomato, label="Sharp kernel")
                end
            end
        end
    end
end
plot_circle!(fig, ax2, radius=1)
axislegend(ax2, position=:cc, framevisible=false)


# magnetisations is : soft/sharp, angle/abs, R, T, sigma, time
zetas_square = (abs.(sum(exp.(4im * angles_square) .* magnitudes_square, dims=2))./sum(magnitudes_square, dims=2))[:, 1, :, :] # soft / sharp, T, sigma
# zetas_triangular = (abs.(sum(exp.(6im * angles_triangular) .* magnitudes_triangular, dims=2))./sum(magnitudes_triangular, dims=2))[:, 1, :, :] # soft / sharp, T, sigma
ax3 = Axis(fig[1, 3]; xlabel=L"σ", ylabel="ζ [%]", aspect=1, ylabelsize=20, title="Alignment")
for j in 1#:length(ratios_TKT)
    data = 100zetas_square[1, j, :]
    scatterlines!(ax3, sigmas, data, marker=:rect, color=:green)
    # data = 100zetas_triangular[1, j, :]
    # scatterlines!(ax3, sigmas, data, marker=:utriangle, color=:green)
    data = 100zetas_square[2, j, :]
    scatterlines!(ax3, sigmas, data, marker=:rect, color=:tomato)
    # data = 100zetas_triangular[2, j, :]
    # scatterlines!(ax3, sigmas, data, marker=:utriangle, color=:tomato)
end
for i in 1:3 
    colsize!(fig.layout, i, Aspect(1, 1))
end
resize_to_layout!(fig)
fig




## ------------------ Scatter Movies comparison ------------------ ##
## ------------------ Scatter Movies comparison ------------------ ##   
## ------------------ Scatter Movies comparison ------------------ ##
## ------------------ Scatter Movies comparison ------------------ ##
using GLMakie
GLMakie.activate!()
markers = [:circle, :cross]
 
mag_x = magnetisations[:, 2, :, :, :, :] .* cos.(magnetisations[:, 1, :, :, :, :])
mag_y = magnetisations[:, 2, :, :, :, :] .* sin.(magnetisations[:, 1, :, :, :, :])


for ind_sig in 1:length(sigmas)
    println("σ = $(sigmas[ind_sig])")

        # initial figure
    fig = Figure(size=(400, 400))
    ax = Axis(fig[1, 1]; limits=(-1.1, 1.1, -1.1, 1.1),
        # xlabel=L"m_x", ylabel=L"m_y",
        aspect=1, xticks=-1:1, yticks=-1:1,
        title=L"σ = %$(round(sigmas[ind_sig], digits=2)) \ ,\ Θ = %$(round(sharp_vision_cones[ind_sig], digits=2))\ ,\ t = %$(round(Int,times[1]/1e4))⋅10^4")
    data_x_soft = Observable(mag_x[1, :, 1, ind_sig, 1])
    data_y_soft = Observable(mag_y[1, :, 1, ind_sig, 1])
    data_x_sharp = Observable(mag_x[2, :, 1, ind_sig, 1])
    data_y_sharp = Observable(mag_y[2, :, 1, ind_sig, 1])
    scatter!(ax, data_x_soft, data_y_soft, markersize=7, marker=:circle, color=:green, label="Smooth kernel")
    scatter!(ax, data_x_sharp, data_y_sharp, markersize=10, marker=:cross, color=:tomato, label="Sharp kernel")
    plot_circle!(fig, ax, radius=1)
    plot_axes!(fig, ax, lattice_type)

    # text!(ax, 1, 1, text=L"t = %$(round(Int,times[4]/1e4))⋅10^4",
    #     fontsize=24, font=:bold_italic,
    #     align=(:right, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)


    text!(ax,1, 0, text=L"m_x",
        fontsize=30, font=:bold_italic,
        align=(:right, :bottom), offset=(-8, 4), space=:relative) # in the coordinate system of the axis)

    text!(ax, 0, 1, text=L"m_y",
        fontsize=30, font=:bold_italic,
        align=(:left, :top), offset=(8, -4), space=:relative) # in the coordinate system of the axis)

    filename = "small_dt_$(lattice_type)_σ_$(round(sigmas[ind_sig], digits=2))_Θ$(round(sharp_vision_cones[ind_sig], digits=2))"
    GLMakie.record(fig, pwd() * "/models/nrxy/movies/comparison_scatter_plot_magnetisation/" * filename * ".mp4", 1:length(times), framerate=15) do tt
        data_x_soft[] = mag_x[1, :, 1, ind_sig, tt]
        data_y_soft[] = mag_y[1, :, 1, ind_sig, tt]
        data_x_sharp[] = mag_x[2, :, 1, ind_sig, tt]
        data_y_sharp[] = mag_y[2, :, 1, ind_sig, tt]
        ax.title = L"σ = %$(round(sigmas[ind_sig], digits=2)) \ ,\ Θ = %$(round(sharp_vision_cones[ind_sig], digits=2))\ ,\ t = %$(round(Int,times[tt]/1e4))⋅10^4"
    end
end

## Coarsening number of defects
fig = Figure()
ax = Axis(fig[1, 1]; xlabel=L"t", ylabel=L"n", aspect=1,
    xscale = log10, yscale = log10)
for j in 1:length(ratios_TKT)
    for k in 3#:length(sigmas)
        for i in 1:2 # soft / sharp
            data = remove_negative(number_of_defects_avg[i, 1, j, k, :])
            scatterlines!(ax, times, data)
        end
    end
end
fig




## ------------------  Comparison sigma Theta soft sharp vision cones  ------------------ ##
## ------------------  Comparison sigma Theta soft sharp vision cones  ------------------ ##   
## ------------------  Comparison sigma Theta soft sharp vision cones  ------------------ ##
## ------------------  Comparison sigma Theta soft sharp vision cones  ------------------ ##
pwd()
filepath = pwd() * "/models/nrxy/data/"
filename = "comparison_soft_sharp_Lx128_Ly128_Rtot480_ratios_TKT[0.2]_tmax1.0e6"
@load filepath * filename * ".jld2" number_of_defects magnetisations Cs xis Lx Ly R R_simus R_tot ratios_TKT sigmas times tmax dt init comments runtime

using CairoMakie
CairoMakie.activate!()
markers = [:circle, :cross]

##
fig = Figure()
ax = Axis(fig[1, 1]; xlabel=L"t", ylabel=L"n", aspect=1,
    xscale = log10#, yscale = log10
)
for j in 1#:length(ratios_TKT)
    for k in [1, 5, 10, 15, 20]#length(sigmas)
        for i in 1:2 # soft / sharp
            data = remove_negative(number_of_defects_avg[i, 1, j, k, :])
            coll = cgrad(:rainbow)[sigmas[k]/sigmas[end]]
            if i == 1 
                lines!(ax, times, data, color=coll, linestyle=:solid)
            elseif i == 2
                lines!(ax, times, data, color=coll, linestyle=:dash)
            end
        end
    end
end
Colorbar(fig[1, 2], colormap=:rainbow, label="σ", colorrange=(0, sigmas[end]))

ax2 = Axis(fig[1, 3]; xlabel=L"σ", ylabel=L"n", aspect=1,)


for j in 1:length(ratios_TKT)
        for indt in round.(Int, length(times) * [1, 3, 5, 10] / 100)
        coll_t = cols_times[(times[indt]-times[1])/(times[Int(end / 100 * 10)]-times[1])]
        data = remove_negative(number_of_defects_avg[1, 1, j, :, indt])
        lines!(ax2, sigmas, data, color=coll_t, linestyle=:solid)
        data = remove_negative(number_of_defects_avg[2, 1, j, :, indt])
        lines!(ax2, sigmas, data, color=coll_t, linestyle=:dash)
    end
end

arrows!(ax2, [0.3], [8], [0.2], [-8], linecolor=:black, arrowcolor=:black)
text!(ax2, 0.3, 8.3, text="time", fontsize=20, color=:black, align=(:center, :bottom))
Colorbar(fig[1, 4], colormap=cols_times, label=L"t\ \  [\times 10^4]", labelrotation=pi/2, colorrange=(times[1] / 1e4, times[Int(end / 100 *10)] / 1e4), minorticksvisible=false)

for i in [1, 3]
    colsize!(fig.layout, i, Aspect(1, 1))
end
resize_to_layout!(fig)
fig
##
fig = Figure()
ax = Axis(fig[1, 1]; 
)
lines!(ax, 0:0.01:1, x -> pi/( (2pi - (16pi*x) / (4pi  + 4x^2 + 8x) ) -pi) -1)
lines!(ax, 0:0.01:1, x -> (pi + 2x) / (pi - 2x) - 1)
lines!(ax, 0:0.01:1, x -> (pi + 2x + pi * x^2 / 4) / (pi - 2x + pi * x^2 / 4) - 1)
fig 

##
