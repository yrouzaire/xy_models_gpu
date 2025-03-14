# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

## ------------------- Non Reciprocal Kuramoto Model ------------------- ##
## ------------------- Non Reciprocal Kuramoto Model ------------------- ##
## ------------------- Non Reciprocal Kuramoto Model ------------------- ##
## ------------------- Non Reciprocal Kuramoto Model ------------------- ##

function kernel_update_NRKuramoto_square!(thetas, thetas_new, omegas, Lx::Int, Ly::Int, R::Int, T::Tf, sigma::Tf, dt::Tf) where {Tf}
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    check = i <= Lx && j <= Ly && k <= R # if L is not a power of 2, we need to check that we are not out of bounds
    if check
        # not slower than differentiating the bulk and the boundaries to avoid unnecessary mod1 usage
        ip = mod1(i + 1, Lx)
        i_ = mod1(i - 1, Lx)
        jp = mod1(j + 1, Ly)
        jm = mod1(j - 1, Ly)
        theta = thetas[i, j, k]
        omega = omegas[i, j, k]
        
        pi_by_2 = Float32(pi / 2)

        #     3
        # 4   X   2   X is the spin of interest. The x direction is ↓, the y direction is → . 
        #     1

        force =
            sin(thetas[ip, j, k] - theta) * g_kernel(theta, 0pi_by_2, sigma) +
            sin(thetas[i, jp, k] - theta) * g_kernel(theta, 1pi_by_2, sigma) +
            sin(thetas[i_, j, k] - theta) * g_kernel(theta, 2pi_by_2, sigma) +
            sin(thetas[i, jm, k] - theta) * g_kernel(theta, 3pi_by_2, sigma)


        thetas_new[i, j, k] = theta + (force/4 + omega) * dt + sqrt(2T * dt) * randn(Tf)
    end

    return nothing
end

function kernel_update_NRKuramoto_triangular!(thetas, thetas_new, omegas, Lx::Int, Ly::Int, R::Int, T::Tf, sigma::Tf, dt::Tf) where {Tf}
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    check = i <= Lx && j <= Ly && k <= R # if L is not a power of 2, we need to check that we are not out of bounds
    if check
        # not slower than differentiating the bulk and the boundaries to avoid unnecessary mod1 usage
        ip = mod1(i + 1, Lx)
        i_ = mod1(i - 1, Lx)
        jp = mod1(j + 1, Ly)
        jm = mod1(j - 1, Ly)
        theta = thetas[i, j, k]
        omega = omegas[i, j, k]

        pi_by_3 = Float32(pi / 3)

        if iseven(i)
            #     4   3
            # 5   X   2   X is the spin of interest. The x direction is ↓, the y direction is → . 
            #     0   1

            force =
                sin(thetas[ip, j, k] - theta) * g_kernel(theta, 0pi_by_3, sigma) +
                sin(thetas[ip, jp, k] - theta) * g_kernel(theta, pi_by_3, sigma) +
                sin(thetas[i, jp, k] - theta) * g_kernel(theta, 2pi_by_3, sigma) +
                sin(thetas[i_, jp, k] - theta) * g_kernel(theta, 3pi_by_3, sigma) +
                sin(thetas[i_, j, k] - theta) * g_kernel(theta, 4pi_by_3, sigma) +
                sin(thetas[i, jm, k] - theta) * g_kernel(theta, 5pi_by_3, sigma)

        elseif isodd(i)
            # 3   2   
            # 4   X   1   X is the spin of interest. The x direction is ↓, the y direction is → . 
            # 5   0   
            force =
                sin(thetas[ip, j, k] - theta) * g_kernel(theta, 0pi_by_3, sigma) +
                sin(thetas[i, jp, k] - theta) * g_kernel(theta, pi_by_3, sigma) +
                sin(thetas[i_, j, k] - theta) * g_kernel(theta, 2pi_by_3, sigma) +
                sin(thetas[i_, jm, k] - theta) * g_kernel(theta, 3pi_by_3, sigma) +
                sin(thetas[i, jm, k] - theta) * g_kernel(theta, 4pi_by_3, sigma) +
                sin(thetas[ip, jm, k] - theta) * g_kernel(theta, 5pi_by_3, sigma)
        end

        thetas_new[i, j, k] = theta + (force/6 + omega) * dt + sqrt(2T * dt) * randn(Tf)
    end

    return nothing
end

function update_NRKuramoto_square!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, dt)
    @cuda threads = block3D blocks = grid3D kernel_update_NRKuramoto_square!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, dt)
    @. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end
function update_NRKuramoto_triangular!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, dt)
    @cuda threads = block3D blocks = grid3D kernel_update_NRKuramoto_triangular!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, dt)
    @. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end

function evolve_NRKuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, t, dt, tmax, lattice_type)
    lattice_type = lowercase(lattice_type)
    @assert lattice_type in ["square", "triangular"]
    if lattice_type == "square"
    while t < tmax
        t += dt
        update_NRKuramoto_square!(thetas, thetas_new, omegas, Lx, Ly, R, T,sigma, dt)
    end

    elseif lattice_type == "triangular"
        while t < tmax
            t += dt
            update_NRKuramoto_triangular!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, dt)
        end
    end
    return thetas_new, t
end
