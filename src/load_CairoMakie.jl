# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine

cols_thetas = cgrad([:black, :blue, :green, :orange, :red, :black]);
cols_sigmas = cgrad(ColorSchemes.viridis)
cols_times = reverse(cgrad(ColorSchemes.Spectral))
cols_py= cgrad(ColorSchemes.tab10)

CairoMakie.activate!()
custom_theme = Theme(
    Axis=(
        xgridvisible=false, ygridvisible=false,
        xtickalign=1, ytickalign=1,
        xticklabelpad=2, yticklabelpad=4,
        xlabelsize = 40, ylabelsize = 40,
    ),
    fontsize=24,
    Colorbar=(ticklabelsize=20, labelrotation=0, labelsize=40,
        tickalign=1,),
    # Legend=(framevisible = false)
    fonts=(; regular=texfont(:text),
        bold=texfont(:bold),
        italic=texfont(:italic),
        bold_italic=texfont(:bolditalic))
)
set_theme!(custom_theme)