# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

# to easily check on the theta field 
function plot_thetas(thetas; defects=false, title="")
    @assert !(isa(thetas, CuArray)) "thetas should not be a CuArray"
    @assert length(size(thetas)) == 2 "thetas should be a 2D Array"

    thetas = mod.(thetas, 2pi)
    fig_without_defects, ax_without_defects = plot_thetas_without_defects(thetas, title=title)
    if defects == false
        return fig_without_defects
    else
        fig_with_defects = circle_defects(thetas, fig_without_defects, ax_without_defects)
        return fig_with_defects
    end
end

function plot_thetas_without_defects(thetas; title="")
    cols_thetas = cgrad([:black, :blue, :green, :orange, :red, :black])

    Lx, Ly = size(thetas)

    fig = Figure() 
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", title = title
        )
    h = heatmap!(ax, thetas[:, :, 1], 
        colormap=cols_thetas, 
        colorrange=(0, 2pi))
    
    Colorbar(fig[1, 2], h, label="θ",
        ticklabelsize=20, labelsize=50, tickalign=1, 
        ticks=([0, pi/2, pi, 3pi/2, 2pi], [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),)
        
    colsize!(fig.layout, 1, Aspect(1, Lx / Ly))
    colgap!(fig.layout, 1, 6)
    resize_to_layout!(fig)

    return fig, ax
end

function circle_defects(thetas, fig, ax)
    Lx, Ly = size(thetas)
    
    thetas = CuArray(reshape(thetas, Lx, Ly, 1)) # to use spot_defects
    dico = spot_defects(thetas)
    qs, xs, ys = dico["qs"], dico["xs"], dico["ys"]

    xs_plus, ys_plus, xs_minus, ys_minus = split_plus_minus(qs, xs, ys)
    xs_plus, ys_plus, xs_minus, ys_minus = xs_plus[1], ys_plus[1], xs_minus[1], ys_minus[1] # to get rid of the R=1 dimension added artificially

    scatter!(ax, xs_plus, ys_plus, marker=:utriangle, strokecolor=:white, color=:transparent, markersize=25, strokewidth=2)
    scatter!(ax, xs_minus, ys_minus, marker=:dtriangle, strokecolor=:white, color=:transparent, markersize=25, strokewidth=2)

    return fig 
end
# plot_thetas(thetas, defects=true)
