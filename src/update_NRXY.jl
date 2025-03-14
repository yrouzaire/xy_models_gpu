# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

## -------------------  NRXY Model ------------------- ##
## -------------------  NRXY Model ------------------- ##
## -------------------  NRXY Model ------------------- ##
## -------------------  NRXY Model ------------------- ##

function g_kernel(theta, uij, sigma)
    return exp(sigma * cos(theta - uij)) # soft vision cone kernel
    # return 1.0 # isotropic = XY
end


function kernel_update_NRXY_square!(thetas, thetas_new, Lx::Int, Ly::Int, R::Int, T::Tf, sigma::Tf, dt::Tf) where {Tf}
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

        pi_by_2 = Float32(pi / 2)

        #     3
        # 4   X   2   X is the spin of interest. The x direction is ↓, the y direction is → . 
        #     1

        force =
            sin(thetas[ip, j, k] - theta) * g_kernel(theta, 0pi_by_2, sigma) +
            sin(thetas[i, jp, k] - theta) * g_kernel(theta, 1pi_by_2, sigma) +
            sin(thetas[i_, j, k] - theta) * g_kernel(theta, 2pi_by_2, sigma) +
            sin(thetas[i, jm, k] - theta) * g_kernel(theta, 3pi_by_2, sigma)

        thetas_new[i, j, k] = theta + force * dt / 4 + sqrt(2T * dt) * randn(Tf)
    end

    return nothing
end

function kernel_update_NRXY_triangular!(thetas, thetas_new, Lx::Int, Ly::Int, R::Int, T::Tf, sigma::Tf, dt::Tf) where {Tf}
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


        thetas_new[i, j, k] = theta + force * dt / 6 + sqrt(2T * dt) * randn(Tf)
    end

    return nothing
end

function update_NRXY_square!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
    @cuda threads = block3D blocks = grid3D kernel_update_NRXY_square!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
    @. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end
function update_NRXY_triangular!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
    @cuda threads = block3D blocks = grid3D kernel_update_NRXY_triangular!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
    @. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end

function evolve_NRXY!(thetas, thetas_new, Lx, Ly, R, T, sigma, t, dt, tmax, lattice_type)
    lattice_type = lowercase(lattice_type)
    @assert lattice_type in ["square", "triangular"]
    if lattice_type == "square"
        while t < tmax
            t += dt
            update_NRXY_square!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
        end
    elseif lattice_type == "triangular"
        while t < tmax
            t += dt
            update_NRXY_triangular!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
        end
    end

    return thetas_new, t
end


## -------------------  Sharp Vision Cone (Langevin Dynamics) ------------------- ##
## -------------------  Sharp Vision Cone (Langevin Dynamics) ------------------- ##
## -------------------  Sharp Vision Cone (Langevin Dynamics) ------------------- ##
## -------------------  Sharp Vision Cone (Langevin Dynamics) ------------------- ##

function get_vision_cone(sigma)
    return Float32(2pi .- 16 * pi * sigma ./ (8sigma .+ pi * (4 .+ sigma .^ 2)))
end

function g_kernel_sharp(theta, uij, vision_cone_aperture)
    pi_2 = 2Float32(pi)
    theta = mod(theta, pi_2)
    uij = mod(uij, pi_2)
    dist = abs(theta - uij)
    dist = min(dist, pi_2 - dist) # Periodic boundary conditions
    return Float32(dist < vision_cone_aperture / 2) # returns a bool, but in julia, true*x = x and false*x = 0 for any number x
end



function kernel_update_NRXY_sharp_square!(thetas, thetas_new, Lx::Int, Ly::Int, R::Int, T::Tf, sigma::Tf, dt::Tf) where {Tf}
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

        pi_by_2 = Float32(pi / 2)

        #     3
        # 4   X   2   X is the spin of interest. The x direction is ↓, the y direction is → . 
        #     1

        force =
            sin(thetas[ip, j, k] - theta) * g_kernel_sharp(theta, 0pi_by_2, sigma) +
            sin(thetas[i, jp, k] - theta) * g_kernel_sharp(theta, 1pi_by_2, sigma) +
            sin(thetas[i_, j, k] - theta) * g_kernel_sharp(theta, 2pi_by_2, sigma) +
            sin(thetas[i, jm, k] - theta) * g_kernel_sharp(theta, 3pi_by_2, sigma)

        thetas_new[i, j, k] = theta + force * dt / 4 + sqrt(2T * dt) * randn(Tf)
    end

    return nothing
end

function kernel_update_NRXY_sharp_triangular!(thetas, thetas_new, Lx::Int, Ly::Int, R::Int, T::Tf, sigma::Tf, dt::Tf) where {Tf}
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

        pi_by_3 = Float32(pi / 3)

        if iseven(i)
            #     4   3
            # 5   X   2   X is the spin of interest. The x direction is ↓, the y direction is → . 
            #     0   1
            force =
                sin(thetas[ip, j, k] - theta) * g_kernel_sharp(theta, 0pi_by_3, sigma) +
                sin(thetas[ip, jp, k] - theta) * g_kernel_sharp(theta, pi_by_3, sigma) +
                sin(thetas[i, jp, k] - theta) * g_kernel_sharp(theta, 2pi_by_3, sigma) +
                sin(thetas[i_, jp, k] - theta) * g_kernel_sharp(theta, 3pi_by_3, sigma) +
                sin(thetas[i_, j, k] - theta) * g_kernel_sharp(theta, 4pi_by_3, sigma) +
                sin(thetas[i, jm, k] - theta) * g_kernel_sharp(theta, 5pi_by_3, sigma)

        elseif isodd(i)
            # 3   2   
            # 4   X   1   X is the spin of interest. The x direction is ↓, the y direction is → . 
            # 5   0   
            force =
                sin(thetas[ip, j, k] - theta) * g_kernel_sharp(theta, 0pi_by_3, sigma) +
                sin(thetas[i, jp, k] - theta) * g_kernel_sharp(theta, pi_by_3, sigma) +
                sin(thetas[i_, j, k] - theta) * g_kernel_sharp(theta, 2pi_by_3, sigma) +
                sin(thetas[i_, jm, k] - theta) * g_kernel_sharp(theta, 3pi_by_3, sigma) +
                sin(thetas[i, jm, k] - theta) * g_kernel_sharp(theta, 4pi_by_3, sigma) +
                sin(thetas[ip, jm, k] - theta) * g_kernel_sharp(theta, 5pi_by_3, sigma)
        end

        thetas_new[i, j, k] = theta + force * dt / 6 + sqrt(2T * dt) * randn(Tf)
    end

    return nothing
end

function update_NRXY_sharp_square!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
    @cuda threads = block3D blocks = grid3D kernel_update_NRXY_sharp_square!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
    @. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end
function update_NRXY_sharp_triangular!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
    @cuda threads = block3D blocks = grid3D kernel_update_NRXY_sharp_triangular!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
    @. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end

function evolve_NRXY_sharp!(thetas, thetas_new, Lx, Ly, R, T, sigma, t, dt, tmax, lattice_type)
    vision_cone_aperture = get_vision_cone(sigma)
    return evolve_NRXY!(thetas, thetas_new, Lx, Ly, R, T, vision_cone_aperture, t, dt, tmax, lattice_type)
end
