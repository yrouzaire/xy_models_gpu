# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`


## -------------- Polarisation ----------------- ##
## -------------- Polarisation ----------------- ##
## -------------- Polarisation ----------------- ##
## -------------- Polarisation ----------------- ##

function OP(thetas) # here thetas can be a 3d CuArray
    tmp = mean(exp.(im * thetas), dims=(1, 2))
    return vec(angle.(tmp)), vec(abs.(tmp))
end

## -------------- Spatial correlation function ----------------- ##
## -------------- Spatial correlation function ----------------- ##
## -------------- Spatial correlation function ----------------- ##
## -------------- Spatial correlation function ----------------- ##

function corr_fft(thetas)
    Lx, Ly, R = size(thetas)
    Lmin_over_2 = round(Int, min(Lx, Ly) / 2)

    C = cos.(thetas)
    S = sin.(thetas)
    c_fft_2d = abs.((ifft(abs.(fft(C) .^ 2 + fft(S) .^ 2))))
    c_fft_2d = Array(c_fft_2d)
    c_fft = zeros(Lmin_over_2 - 1, R)
    for r in 1:R
        c_fft[:, r] = (c_fft_2d[2:Lmin_over_2, 1, r] + c_fft_2d[1, 2:Lmin_over_2, r]) / (2c_fft_2d[1, 1, r])
    end
    # does the mean of the first row and first column and renormalises by the first element
    return c_fft
end

function correlation_length(C, threshold=exp(-1))
    L, R = size(C)
    xi = zeros(R)
    for r in 1:R
        tmp = findfirst(x -> x < threshold, C[:, r])
        if isnothing(tmp)
            xi[r] = NaN
        else
            xi[r] = tmp
        end
    end
    return xi
end


## -------------- Topological defects ----------------- ##
## -------------- Topological defects ----------------- ##
## -------------- Topological defects ----------------- ##
## -------------- Topological defects ----------------- ##

#= Note: arclength(theta1, theta2) is defined in src/auxiliary.jl. 
It returns the signed arclength (in radians) from theta1 to theta2 on the unit trigonometric circle. =#


function spot_defects(thetas, lattice_type="square")
    @assert length(size(thetas)) == 3 "thetas should be a 3D array"
    thetas = CuArray(thetas)
    Lx, Ly, R = size(thetas)
    there_is_a_plus_defect = CUDA.zeros(Int, Lx, Ly, R) # predeclaration
    there_is_a_minus_defect = CUDA.zeros(Int, Lx, Ly, R) # predeclaration
    if lattice_type == "square"
        @cuda threads = block3D blocks = grid3D kernel_spot_defects_square!(mod.(thetas, 2Float32(pi)), there_is_a_plus_defect, there_is_a_minus_defect, Lx, Ly, R)
    elseif lattice_type == "triangular"
        @cuda threads = block3D blocks = grid3D kernel_spot_defects_triangular!(mod.(thetas, 2Float32(pi)), there_is_a_plus_defect, there_is_a_minus_defect, Lx, Ly, R)
    else
        error("lattice_type should be one of : square, triangular")
    end

    return there_is_a_plus_defect, there_is_a_minus_defect
end

function number_defects(thetas, lattice_type="square")
    there_is_a_plus_defect, there_is_a_minus_defect = spot_defects(thetas, lattice_type)
    return reduce(+, there_is_a_plus_defect, dims=(1, 2))[1, 1, :] + reduce(+, there_is_a_minus_defect, dims=(1, 2))[1, 1, :]
end


function number_defects_plus_minus(thetas, lattice_type="square")
    there_is_a_plus_defect, there_is_a_minus_defect = spot_defects(thetas, lattice_type)
    n_plus = reduce(+, there_is_a_plus_defect, dims=(1, 2))[1, 1, :]
    n_minus = reduce(+, there_is_a_minus_defect, dims=(1, 2))[1, 1, :]
    n_total = n_plus + n_minus
    return n_plus, n_minus, n_total
end


function kernel_spot_defects_square!(thetas, there_is_a_plus_defect, there_is_a_minus_defect, Lx::Int, Ly::Int, R::Int)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    # check = i <= Lx && j <= Ly && k <= R # if L is not a power of 2, we need to check that we are not out of bounds
    check_bulk = i > 1 && j > 1 && i < Lx && j < Ly && k <= R
    if check_bulk

        #= The plaquette is defined as follows (X is the defect)
        i,j+1 -------i+1,j+1
        |           |
        |     X     |
        |           |
        i,j-------i+1,j

        =#

        # not slower than differentiating the bulk and the boundaries to avoid unnecessary mod1 usage
        ip = mod1(i + 1, Lx)
        jp = mod1(j + 1, Ly)

        # Compute the charge q of the topological defect
        q = Float32(0)
        q += arclength(thetas[i, j, k], thetas[ip, j, k])
        q += arclength(thetas[ip, j, k], thetas[ip, jp, k])
        q += arclength(thetas[ip, jp, k], thetas[i, jp, k])
        q += arclength(thetas[i, jp, k], thetas[i, j, k])
        q = round(q / 2 / pi, digits=1) # to avoid numerical errors such as 0.99999 ≠ 1

        if q == 1
            there_is_a_plus_defect[i, j, k] = true
        elseif q == -1
            there_is_a_minus_defect[i, j, k] = true
        end

    end

    return nothing
end



function kernel_spot_defects_triangular!(thetas, there_is_a_plus_defect, there_is_a_minus_defect, Lx::Int, Ly::Int, R::Int)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    # check = i <= Lx && j <= Ly && k <= R # if L is not a power of 2, we need to check that we are not out of bounds
    check_bulk = i > 1 && j > 1 && i < Lx && j < Ly && k <= R
    if check_bulk

        # not slower than differentiating the bulk and the boundaries to avoid unnecessary mod1 usage 
        imm = mod1(i - 1, Lx)
        ip = mod1(i + 1, Lx)
        jm = mod1(j - 1, Ly)
        jp = mod1(j + 1, Ly)


        if iseven(i)
            #= The up and down plaquettes are defined as follows (X is the defect)
               i+1,j+1               
                /   \       
               /  X  \                 
              /       \                 
            i,j-------i+1,j             
              \       /
               \  X  /
                \   /
               i-1,j+1
            =#

            # up plaquette
            q = Float32(0)
            q += arclength(thetas[i, j, k], thetas[ip, j, k])
            q += arclength(thetas[ip, j, k], thetas[ip, jp, k])
            q += arclength(thetas[ip, jp, k], thetas[i, j, k])
            q = round(q / 2 / pi, digits=1) # to avoid numerical errors such as 0.99999 ≠ 1

            if q == 1
                there_is_a_plus_defect[i, j, k] = true
            elseif q == -1
                there_is_a_minus_defect[i, j, k] = true
            end

            # down plaquette
            q = Float32(0)
            q += arclength(thetas[i, j, k], thetas[imm, jp, k])
            q += arclength(thetas[imm, jp, k], thetas[ip, j, k])
            q += arclength(thetas[ip, j, k], thetas[i, j, k])
            q = round(q / 2 / pi, digits=1) # to avoid numerical errors such as 0.99999 ≠ 1

            if q == 1
                there_is_a_plus_defect[i, j, k] = true
            elseif q == -1
                there_is_a_minus_defect[i, j, k] = true
            end

        elseif isodd(i)
            #= The up and down plaquettes are defined as follows (X is the defect)

                i-1,j+1 ----- i,j+1            
                    \       /
                     \  X  /
                      \   /
                       i,j
                      /   \       
                     /  X  \                 
                    /       \ 
                i-1,j-1 ----- i,j-1
            =#


            # down plaquette (the triangle pointing down, above)
            q = Float32(0)
            q += arclength(thetas[imm, jp, k], thetas[i, j, k])
            q += arclength(thetas[i, j, k], thetas[ip, jp, k])
            q += arclength(thetas[i, jp, k], thetas[imm, jp, k])
            q = round(q / 2 / pi, digits=1) # to avoid numerical errors such as 0.99999 ≠ 1

            if q == 1
                there_is_a_plus_defect[i, j, k] = true
            elseif q == -1
                there_is_a_minus_defect[i, j, k] = true
            end

            # up plaquette (the triangle pointing up, below)
            q = Float32(0)
            q += arclength(thetas[i, j, k], thetas[imm, jm, k])
            q += arclength(thetas[imm, jm, k], thetas[i, jm, k])
            q += arclength(thetas[i, jm, k], thetas[i, j, k])
            q = round(q / 2 / pi, digits=1) # to avoid numerical errors such as 0.99999 ≠ 1

            if q == 1
                there_is_a_plus_defect[i, j, k] = true
            elseif q == -1
                there_is_a_minus_defect[i, j, k] = true
            end


        end

    end


    return nothing
end

