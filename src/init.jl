# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

function create_thetas(Lx, Ly, R, init::String, params_init)

    if init in ["random", "hightemp", "disorder", "disordered"]
        thetas = Float32(2pi) * rand(Float32, Lx, Ly, R)
    elseif init in ["order", "ordered", "lowtemp"]
        thetas = zeros(Float32, Lx, Ly, R)
    elseif init in ["pair"]
        @unpack r0, mu_plus, phi = params_init # phi is the angle that the +/- line makes with the horizontal.

        @assert 0 < r0 <= Lx / 2 "The radius r0 should be smaller than Lx/2"
        @assert 0 < r0 <= Ly / 2 "The radius r0 should be smaller than Ly/2"

        x_plus = Lx / 2 + r0 / 2 * cos(phi)
        y_plus = Ly / 2 + r0 / 2 * sin(phi)

        x_minus = Lx / 2 - r0 / 2 * cos(phi)
        y_minus = Ly / 2 - r0 / 2 * sin(phi)

        thetas = zeros(Float32, Lx, Ly, R)
        for i in 1:Lx, j in 1:Ly

            dx_plus = i - x_plus
            # dx_plus = argmin(abs, [dx_plus, Lx - dx_plus, Lx + dx_plus])
            dy_plus = j - y_plus
            # dy_plus = argmin(abs, [dy_plus, Ly - dy_plus, Ly + dy_plus])

            dx_minus = i - x_minus
            # dx_minus = argmin(abs, [dx_minus, Lx - dx_minus, Lx + dx_minus])
            dy_minus = j - y_minus
            # dy_minus = argmin(abs, [dy_minus, Ly - dy_minus, Ly + dy_minus])

            thetas[i, j, :] .= mu_plus + phi +
                               atan(dy_plus, dx_plus) - atan(dy_minus, dx_minus) +
                               atan(dy_plus + Ly, dx_plus) - atan(dy_minus + Ly, dx_minus) +
                               atan(dy_plus - Ly, dx_plus) - atan(dy_minus - Ly, dx_minus) +
                               atan(dy_plus, dx_plus + Lx) - atan(dy_minus, dx_minus + Lx) +
                               atan(dy_plus, dx_plus - Lx) - atan(dy_minus, dx_minus - Lx)
        end
    elseif init in ["single"]
        @unpack q, mu0 = params_init
        x_center = Lx / 2
        y_center = Ly / 2
        thetas = zeros(Float32, Lx, Ly, R)

        for i in 1:Lx, j in 1:Ly
            dx = i - x_center
            dy = j - y_center
            thetas[i, j, :] .= mu0 + q * atan(dy, dx)
        end
    else
        accepted_inits = ["random", "hightemp", "disorder", "disordered",
            "ordered", "order", "lowtemp",
            "pair", "single"]
        error("init should be one of : $(join(accepted_inits, ", "))")
    end

    return CuArray(thetas)
end
initialisation_thetas = create_thetas # alias

function create_omegas(Lx, Ly, R, σ, distribution_type)
    distribution_type = lowercase(distribution_type)
    if distribution_type in ["gaussian", "normal", "norm"]
        sigmas = σ * randn(Lx, Ly, R)
    elseif distribution_type in ["uniform", "unif"]
        # The variance of a uniform distribution between -a and a is a^2/3. 
        # If we want the same variance as the gaussian, we should have a = σ*sqrt(3)
        sigmas = σ * sqrt(3) * (2 * rand(Lx, Ly, R) .- 1)
    elseif distribution_type in ["laplace", "exponential"]
        #= Here we implement the symmetrised exponential distribution, aka the Laplace distribution : https://en.wikipedia.org/wiki/Laplace_distribution
        The pdf is f(x|μ,b) = 1/(2b) * exp(-|x-μ|/b)
        The variance is 2b^2, so we should have b = σ/sqrt(2) =#
        if σ > 0
            b = σ / sqrt(2)
            sigmas = rand(Laplace(0, b), Lx, Ly, R)
        else
            sigmas = zeros(Lx, Ly, R)
        end
    elseif distribution_type in ["lorentzian", "cauchy"]
        #= Here we implement the Cauchy distribution : https://en.wikipedia.org/wiki/Cauchy_distribution
        The pdf is f(x|x0,γ) = 1/(πγ * (1 + ((x-x0)/γ)^2))
        The variance is ∞, so we should have γ = σ =#
        if σ > 0
            sigmas = rand(Cauchy(0, σ), Lx, Ly, R)
        else
            sigmas = zeros(Lx, Ly, R)
        end
    elseif distribution_type in ["trunc_lorentzian", "truncated_lorentzian", "truncated_cauchy", "trunc_cauchy"]
        if σ > 0
            distrib = Truncated(Cauchy(0, σ), -5σ, 5σ) # the bounds 5σ are arbitrary
            sigmas = rand(distrib, Lx, Ly, R)
        else
            sigmas = zeros(Lx, Ly, R)
        end
    else
        println("Error : distribution_type should be \"gaussian\" or \"uniform\" or \"laplace\" or \"lorentzian\".")
    end
    return CuArray(Float32.(sigmas))
end
initialisation_omegas = create_omegas # alias