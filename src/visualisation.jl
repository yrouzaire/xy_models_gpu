# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

# to easily check on the theta field 

function plot_thetas(thetas; defects=false, lattice_type="square")
    fig_without_defects, ax_without_defects = plot_thetas_without_defects(thetas)
    if defects == false
        return fig_without_defects
    else
        fig_with_defects = circle_defects(thetas, fig_without_defects, ax_without_defects, lattice_type)
        return fig_with_defects
    end
end

function plot_thetas_without_defects(thetas_format_unknown)
    cols_thetas = cgrad([:black, :blue, :green, :orange, :red, :black])

    thetas = Array(thetas_format_unknown[:,:,1]) # convert to Array if it's a CuArray and select only one realisation 
    Lx, Ly = size(thetas)
    
    thetas = mod.(thetas, 2pi)

    fig = Figure() 
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", 
        )
    h = heatmap!(ax, thetas[:, :, 1], 
        colormap=cols_thetas, 
        colorrange=(0, 2pi), )
    
    Colorbar(fig[1, 2], h, label="θ",
        ticklabelsize=20, labelsize=50, tickalign=1, 
        ticks=([0, pi/2, pi, 3pi/2, 2pi], [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),)
        
    colsize!(fig.layout, 1, Aspect(1, Lx / Ly))
    resize_to_layout!(fig)

    return fig, ax
end


function circle_defects(thetas, fig, ax, lattice_type)
    there_is_a_plus_defect, there_is_a_minus_defect = spot_defects(thetas, lattice_type)
    plus_defects = Array(there_is_a_plus_defect[:, :, 1])
    minus_defects = Array(there_is_a_minus_defect[:, :, 1])

    thetas = Array(thetas)
    Lx, Ly, R = size(thetas)
    

    nb_defects_plus = sum(plus_defects)
    nb_defects_minus = sum(minus_defects)
    if nb_defects_plus !== nb_defects_minus 
        ax.title = "Defects missing"
    end

    xys_plus = zeros(2, nb_defects_plus)
    xys_minus = zeros(2, nb_defects_minus)

    mm_plus = 1
    mm_minus = 1
    for i in 1:Lx, j in 1:Ly
        if plus_defects[i, j] == 1
            xys_plus[1, mm_plus] = i
            xys_plus[2, mm_plus] = j
            mm_plus += 1
        end

        if minus_defects[i, j] == 1
            xys_minus[1, mm_minus] = i
            xys_minus[2, mm_minus] = j
            mm_minus += 1
        end
    end
            
    xys_plus = xys_plus .+ 0.5

    scatter!(ax, xys_plus[1, :], xys_plus[2, :], marker=:utriangle, strokecolor=:white, color=:transparent, markersize=25, strokewidth=2)
    scatter!(ax, xys_minus[1, :], xys_minus[2, :], marker=:dtriangle, strokecolor=:white, color=:transparent, markersize=25, strokewidth=2)

    return fig 
end
# plot_thetas(thetas, defects=true)
