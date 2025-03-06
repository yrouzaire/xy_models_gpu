# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

include("src/load_everything.jl")

## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
## ------------------ Parameters and data generation ------------------ ##
lattice_type = "triangular"

L = 128 # For this script, Lx must be equal to Ly
Lx = Ly = L
R = 2 # number of independent realisations
R_simus = 1
R_tot = R * R_simus

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)

tmax = Tf(1E5)
dt = Tf(1E-1)
time_spacing = "linear" # choose among "linear", "log", "quadratic"
if time_spacing == "linear"
    every = tmax / 100
    times = Tf.(collect(every:every:tmax)) # linear
elseif time_spacing == "log"
    times = logspace(10dt, tmax, nb_frames) # log
elseif time_spacing == "quadratic"
    times = Tf.(range(sqrt(dt), sqrt(tmax), nb_frames)) .^ 2  # quadratic
else
    error("time_spacing should be 'linear', 'log' or 'quadratic'")
end

init = "hightemp"
params_init = ("dummy") # only useful for the "pair" init (a pair of defects)

sigmas = [0, 0.15, 0.3]
# sigmas = collect(0:0.05:0.30)
sharp_vision_cones = Tf.(2pi .- 16 * pi * sigmas ./ (8sigmas .+ pi * (4 .+ sigmas .^ 2))) # equivalence between sharp and smooth vision cones

ratios_TKT = [0.2]
threshold_magnetisation = 1.9 # if larger than 1, it's like if it was desactivated. 

# the first index is the smooth vision cone model, the second is the sharp vision cone model
magnetisations = zeros(2, 2, R_tot, length(ratios_TKT), length(sigmas), length(times))
# thetas_saved_cpu = zeros(Float16, 2, Lx, Ly, R_tot, length(ratios_TKT), length(sigmas), length(times)) # in Float16 for storage reasons
# number_of_defects = zeros(2, R_tot, length(ratios_TKT), length(sigmas), length(times))
# Cs = zeros(2, round(Int, L / 2 - 1), R_tot, length(ratios_TKT), length(sigmas), length(times)) # correlation function C(r,t)
# xis = zeros(2, R_tot, length(ratios_TKT), length(sigmas), length(times)) # correlation length ξ(t)

m = 0
M = length(ratios_TKT) * length(sigmas) * R_simus * 2 # 2 is the number of models

z = @elapsed for i in each(ratios_TKT), j in each(sigmas), r in 1:R_simus
    for ind_model in 1:2 # 1 is soft, 2 is sharp

        T = Tf(round(ratios_TKT[i] * 0.22, digits=4))
        thetas = create_thetas(Lx, Ly, R, init, params_init)
        thetas_new = similar(thetas)
        t = Float64(0)

        global m += 1

        if ind_model == 1 # # # # #  Smooth Vision Cone # # # # # 

            sigma = Tf(sigmas[j])
            println("Simulation $m/$M : T = $T, σ = $(sigma)")

            for tt in eachindex(times)
                thetas, t = evolve_NRXY!(thetas, thetas_new, Lx, Ly, R, T, sigma, t, dt, times[tt], lattice_type)

                # thetas_saved_cpu[1, :, :, 1+(r-1)*R:(r*R), i, j, tt] = Array(thetas)
                
                tmp1 = OP(thetas)
                magnetisations[ind_model, 1, 1+(r-1)*R:(r*R), i, j, tt] = Array(tmp1[1])
                magnetisations[ind_model, 2, 1+(r-1)*R:(r*R), i, j, tt] = Array(tmp1[2])

                # dico = spot_defects(thetas)
                # number_of_defects[ind_model, 1+(r-1)*R:(r*R), i, j, tt] = dico["number_defects"]

                # tmp2 = Array(corr_fft(thetas))
                # Cs[ind_model, :, 1+(r-1)*R:(r*R), i, j, tt] = tmp2
                # xis[ind_model, 1+(r-1)*R:(r*R), i, j, tt] = correlation_length(tmp2)

                if mean(tmp1[2] .> threshold_magnetisation) == 1
                    println("Soft model has reached the threshold magnetisation at time $(times[tt])")
                    magnetisations[ind_model, 1, 1+(r-1)*R:(r*R), i, j, tt+1:end] .= Array(tmp1[1])
                    magnetisations[ind_model, 2, 1+(r-1)*R:(r*R), i, j, tt+1:end] .= Array(tmp1[2])

                    # Cs[ind_model, :, 1+(r-1)*R:(r*R), i, j, tt+1:end] .= tmp2
                    # xis[ind_model, 1+(r-1)*R:(r*R), i, j, tt+1:end] .= correlation_length(tmp2)

                    break
                end
            end
        elseif ind_model == 2 # # # # #  Sharp Vision Cone # # # # # 

            sigma = sharp_vision_cones[j]
            println("Simulation $m/$M : T = $T, Θ = $(sigma)")

            for tt in eachindex(times)
                thetas, t = evolve_NRXY_sharp!(thetas, thetas_new, Lx, Ly, R, T, sigma, t, dt, times[tt], lattice_type)

                # thetas_saved_cpu[ind_model, :, :, 1+(r-1)*R:(r*R), i, j, tt] = Array(thetas)
                
                tmp1 = OP(thetas)
                magnetisations[ind_model, 1, 1+(r-1)*R:(r*R), i, j, tt] = Array(tmp1[1])
                magnetisations[ind_model, 2, 1+(r-1)*R:(r*R), i, j, tt] = Array(tmp1[2])

                # dico = spot_defects(thetas)
                # number_of_defects[ind_model, 1+(r-1)*R:(r*R), i, j, tt] = dico["number_defects"]

                # tmp2 = Array(corr_fft(thetas))
                # Cs[ind_model, :, 1+(r-1)*R:(r*R), i, j, tt] = tmp2
                # xis[ind_model, 1+(r-1)*R:(r*R), i, j, tt] = correlation_length(tmp2)

                if mean(tmp1[2] .> threshold_magnetisation) == 1
                    println("Sharp model has reached the threshold magnetisation at time $(times[tt])")
                    magnetisations[ind_model, 1, 1+(r-1)*R:(r*R), i, j, tt+1:end] .= Array(tmp1[1])
                    magnetisations[ind_model, 2, 1+(r-1)*R:(r*R), i, j, tt+1:end] .= Array(tmp1[2])

                    # Cs[ind_model, :, 1+(r-1)*R:(r*R), i, j, tt+1:end] .= tmp2
                    # xis[ind_model, 1+(r-1)*R:(r*R), i, j, tt+1:end] .= correlation_length(tmp2)

                    break
                end
            end
        end
    end
end
println("Simulation over.")
prinz(z) # prints the runtime in seconds, minutes and hours. Defined in `auxiliary.jl` file

# plot_thetas(thetas_saved_cpu[:, :, 1, 1, 2,end], defects=true) # for a quick visual check if needed 

## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##
## ------------------ Save data ------------------ ##

comments = "$(time_spacing)ly spaced times. $(lattice_type) lattice."

pwd()
filepath = pwd() #* "/models/nrxy/data/"
filename = "/$(lattice_type)_dynamics_comparison_soft_sharp_Lx$(Lx)_Ly$(Ly)_Rtot$(R_tot)_ratios_TKT$(ratios_TKT)_tmax$(tmax)"

@save filepath * filename * ".jld2" magnetisations Lx Ly R R_simus R_tot ratios_TKT sigmas times tmax dt init comments lattice_type runtime = z
# @save filepath * filename * "_with_thetas_final_time.jld2" number_of_defects magnetisations Cs xis thetas_saved_cpu Lx Ly R R_simus R_tot ratios_TKT sigmas tmax times dt init comments runtime = z

println("Data saved in $(filepath * filename).jld2")


## ------------------ Plots ------------------ ##
## ------------------ Plots ------------------ ##   
## ------------------ Plots ------------------ ##
## ------------------ Plots ------------------ ##
# fig = Figure(size=(400, 400)) 
# ax = Axis(fig[1, 1]; limits = (-1.1, 1.1, -1.1, 1.1), xlabel = "m_x", ylabel = "m_y", aspect=1)
# for j in 1:length(ratios_TKT)
#     for k in 2#:length(sigmas)
#         for tt in length(times)
#             for i in 1:2 # soft / sharp
#                 mag_x = magnetisations[i, 2, :, j, k, tt] .* cos.(magnetisations[i, 1, :, j, k, tt])
#                 mag_y = magnetisations[i, 2, :, j, k, tt] .* sin.(magnetisations[i, 1, :, j, k, tt])
#                 scatter!(ax, mag_x, mag_y, markersize = 10, 
#                     # color=,
#                     # color=cols_sigmas[sigmas[k]/sigmas[end]],
#                     )
#             end
#         end
#     end
# end
# Colorbar(fig[1, 2], colormap=cols_sigmas, label="σ", colorrange=(0, sigmas[end]))
# fig


# ##
# indr = rand(1:R_tot)
# indsig = 1
# titre = string(round(magnetisations[1, 2, indr, 1, indsig, end], digits=2)) * " " * string(round(magnetisations[1, 1, indr, 1, indsig, end], digits=2))
# plot_thetas(thetas_saved_cpu[1, :, :, indr, 1, indsig, end], defects=true, title=titre) # for a quick visual check if needed 
