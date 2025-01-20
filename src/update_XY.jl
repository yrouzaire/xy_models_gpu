# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

## ------------------- Vanilla XY Model ------------------- ##
## ------------------- Vanilla XY Model ------------------- ##
## ------------------- Vanilla XY Model ------------------- ##
## ------------------- Vanilla XY Model ------------------- ##

function kernel_update_XY_square!(thetas, thetas_new, Lx::Int, Ly::Int, R::Int, T::Tf, dt::Tf) where {Tf}
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

        force = sin(thetas[i, jp, k] - theta) +
                sin(thetas[ip, j, k] - theta) +
                sin(thetas[i, jm, k] - theta) +
                sin(thetas[i_, j, k] - theta)

        thetas_new[i, j, k] = theta + force * dt / 4 + sqrt(2T * dt) * randn(Tf)
    end

    return nothing
end

function kernel_update_XY_triangular!(thetas, thetas_new, Lx::Int, Ly::Int, R::Int, T::Tf, dt::Tf) where {Tf}
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

        if iseven(i)
            force = sin(thetas[i, jp, k] - theta) +
                    sin(thetas[i_, jp, k] - theta) +
                    sin(thetas[i_, j, k] - theta) +
                    sin(thetas[i, jm, k] - theta) +
                    sin(thetas[ip, j, k] - theta) +
                    sin(thetas[ip, jp, k] - theta)
        elseif isodd(i)
            force = sin(thetas[i, jp, k] - theta) +
                    sin(thetas[i_, j, k] - theta) +
                    sin(thetas[i_, jm, k] - theta) +
                    sin(thetas[i, jm, k] - theta) +
                    sin(thetas[ip, jm, k] - theta) +
                    sin(thetas[ip, j, k] - theta)
        end

        thetas_new[i, j, k] = theta + force * dt / 6 + sqrt(2T * dt) * randn(Tf)
    end

    return nothing
end

function update_XY_square!(thetas, thetas_new, Lx, Ly, R, T, dt)
    @cuda threads = block3D blocks = grid3D kernel_update_XY_square!(thetas, thetas_new, Lx, Ly, R, T, dt)
    @. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end
function update_XY_triangular!(thetas, thetas_new, Lx, Ly, R, T, dt)
    @cuda threads = block3D blocks = grid3D kernel_update_XY_triangular!(thetas, thetas_new, Lx, Ly, R, T, dt)
    @. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end

function evolve_XY!(thetas, thetas_new, Lx, Ly, R, T, t, dt, tmax, lattice_type)
    @assert lattice_type in ["square", "triangular"]
    while t < tmax
        t += dt

        if lattice_type == "square"
            update_XY_square!(thetas, thetas_new, Lx, Ly, R, T, dt)
        elseif lattice_type == "triangular"
            update_XY_triangular!(thetas, thetas_new, Lx, Ly, R, T, dt)
        end
    end
    return thetas_new, t
end
