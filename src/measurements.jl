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

function spot_defects(thetas)
    @assert length(size(thetas)) == 3 "thetas should be a 3D array"
    thetas = CuArray(thetas)
    thetas = mod.(thetas, 2Float32(pi))
    Lx, Ly, R = size(thetas)

    qs_matrix = CUDA.zeros(Int, Lx, Ly, R) # topological charges
    #= 
    One could define  
        xs_matrix = CUDA.zeros(Int, Lx, Ly, R) # position x
        ys_matrix = CUDA.zeros(Int, Lx, Ly, R) # position y
    but, because then one would juste have to write 
        xs_matrix[i, j, k] = i
        ys_matrix[i, j, k] = j 
    let's just create those matrices manually. 
    =#
    xs_matrix = repeat(reshape(1:Lx, Lx, 1, 1), 1, Ly, R)
    ys_matrix = repeat(reshape(1:Ly, 1, Ly, 1), Lx, 1, R)
    xs_vector_tmp = reshape(xs_matrix, Lx * Ly, R)
    ys_vector_tmp = reshape(ys_matrix, Lx * Ly, R)

    @cuda threads = block3D blocks = grid3D kernel_spot_defects(thetas, qs_matrix, Lx, Ly, R)
    
    qs_vector_tmp = Array(reshape(qs_matrix, Lx * Ly, R))

    qs = Vector{Vector{Float64}}(undef, R)
    xs = Vector{Vector{Float64}}(undef, R)
    ys = Vector{Vector{Float64}}(undef, R)

    for r in 1:R
        indices_defects = findall(qs_vector_tmp[:, r] .!= 0)
        qs[r] = qs_vector_tmp[indices_defects, r]
        xs[r] = xs_vector_tmp[indices_defects, r]
        ys[r] = ys_vector_tmp[indices_defects, r]
    end
    number_defects_plus = [sum(el .== 1) for el in qs]
    number_defects_minus = [sum(el .== -1) for el in qs]
    number_defects = number_defects_plus .+ number_defects_minus

    dico = Dict()
    dico["qs"] = qs
    dico["xs"] = xs
    dico["ys"] = ys
    dico["number_defects_plus"] = number_defects_plus
    dico["number_defects_minus"] = number_defects_minus
    dico["number_defects"] = number_defects

    return dico
end


split_plus_minus(dico) = split_plus_minus(dico["qs"], dico["xs"], dico["ys"])
function split_plus_minus(qs, xs, ys)
    R = length(qs)

    indices_plus = [findall(q .== 1) for q in qs]
    indices_minus = [findall(q .== -1) for q in qs]

    xs_plus = [xs[r][indices_plus[r]] for r in 1:R]
    ys_plus = [ys[r][indices_plus[r]] for r in 1:R]
    xs_minus = [xs[r][indices_minus[r]] for r in 1:R]
    ys_minus = [ys[r][indices_minus[r]] for r in 1:R]

    return xs_plus, ys_plus, xs_minus, ys_minus
end


function kernel_spot_defects(thetas, qs_matrix, Lx::Int, Ly::Int, R::Int)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    check = i <= Lx && j <= Ly && k <= R # if L is not a power of 2, we need to check that we are not out of bounds
    if check

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

        qs_matrix[i, j, k] = q
    end

    return nothing
end


## -------------- Structure Factor ----------------- ##
## -------------- Structure Factor ----------------- ##
## -------------- Structure Factor ----------------- ##
## -------------- Structure Factor ----------------- ##

distance_matrix_defects(dico, L) = distance_matrix_defects(split_plus_minus(dico)..., L, R)
function distance_matrix_defects(xs_plus, ys_plus, xs_minus, ys_minus, L, R)
    Ns = length.(xs_plus)
    @assert Ns == length.(xs_minus) "The number of +/- defects should be the same for each realisation."

    distance_matrices_all = Vector{Matrix{Float64}}(undef, R)
    distance_matrices_plus_minus = Vector{Matrix{Float64}}(undef, R)
    distance_matrices_plus_plus = Vector{Matrix{Float64}}(undef, R)
    distance_matrices_minus_minus = Vector{Matrix{Float64}}(undef, R)
    for r in 1:R
        N = Ns[r]

        # all defects (no matter the charge)
        distance_matrices_all[r] = construct_distance_matrix(2N, vcat(xs_plus[r], xs_minus[r]), vcat(ys_plus[r], ys_minus[r]), L)

        # +/+ defects (only between positive charges)
        distance_matrices_plus_plus[r] = construct_distance_matrix(N, xs_plus[r], ys_plus[r], L)

        # -/- defects (only between negative charges)
        distance_matrices_minus_minus[r] = construct_distance_matrix(N, xs_minus[r], ys_minus[r], L)

        # +/- defects (opposite charges)
        distance_matrices_plus_minus[r] = zeros(N, N)
        for i in 1:N, j in 1:N # we have to do it manually because we mix infos between 2 separate arrays and it's therefore non-symmetric
            dx = abs(xs_plus[r][i] - xs_minus[r][j])
            dx = min(dx, L - dx)

            dy = abs(ys_plus[r][i] - ys_minus[r][j])
            dy = min(dy, L - dy)

            distance_matrices_plus_minus[r][i, j] = sqrt(dx^2 + dy^2)
        end
    end

    return distance_matrices_all, distance_matrices_plus_minus, distance_matrices_plus_plus, distance_matrices_minus_minus
end


function delta_matrix_defects(xs_plus, ys_plus, xs_minus, ys_minus, L, R)
    Ns = length.(xs_plus)
    @assert Ns == length.(xs_minus) "The number of +/- defects should be the same for each realisation."

    delta_matrices_all = Vector{Matrix{Float64}}(undef, R)
    delta_matrices_plus_minus = Vector{Matrix{Float64}}(undef, R)
    delta_matrices_plus_plus = Vector{Matrix{Float64}}(undef, R)
    delta_matrices_minus_minus = Vector{Matrix{Float64}}(undef, R)
    for r in 1:R
        N = Ns[r]

        # all defects (no matter the charge)
        delta_matrices_all[r] = construct_delta_matrix(2N, vcat(xs_plus[r], xs_minus[r]), vcat(ys_plus[r], ys_minus[r]), L)

        # +/+ defects (only between positive charges)
        delta_matrices_plus_plus[r] = construct_delta_matrix(N, xs_plus[r], ys_plus[r], L)

        # -/- defects (only between negative charges)
        delta_matrices_minus_minus[r] = construct_delta_matrix(N, xs_minus[r], ys_minus[r], L)

        # +/- defects (opposite charges)
        delta_matrices_plus_minus[r] = zeros(N, N)
        for i in 1:N, j in 1:N # we have to do it manually because we mix infos between 2 separate arrays and it's therefore non-symmetric
            dx = (xs_plus[r][i] - xs_minus[r][j])
            dx = argmin(abs, [dx, dx + L, dx - L])

            dy = (ys_plus[r][i] - ys_minus[r][j])
            dy = argmin(abs, [dy, dy + L, dy - L])

            delta_matrices_plus_minus[r][i, j] = (dx + dy)
        end
    end

    return delta_matrices_all, delta_matrices_plus_minus, delta_matrices_plus_plus, delta_matrices_minus_minus
end

function construct_delta_matrix(N, xs, ys, L) 
    # element_ij = (x_i - x_j) + (y_i - y_j)
    # WARNING ! It's not a distance ! (it is signed)

    delta_matrix = zeros(N, N)
    for i in 1:N, j in i+1:N # only the upper triangular part. The matrix is symmetric. By definition, the diagonal is 0.
        dx = (xs[i] - xs[j])
        dx = argmin(abs, [dx, dx + L, dx - L])

        dy = (ys[i] - ys[j])
        dy = argmin(abs, [dy, dy + L, dy - L])

        delta = dx + dy
        delta_matrix[i, j] = delta
        delta_matrix[j, i] = -delta
    end

    return delta_matrix
end


function construct_distance_matrix(N, xs, ys, L)
    distance_matrix = zeros(N, N)
    for i in 1:N, j in i+1:N # only the upper triangular part. The matrix is symmetric. By definition, the diagonal is 0.
        dx = abs(xs[i] - xs[j])
        dx = min(dx, L - dx)

        dy = abs(ys[i] - ys[j])
        dy = min(dy, L - dy)

        dist = sqrt(dx^2 + dy^2)
        distance_matrix[i, j] = dist
        distance_matrix[j, i] = dist
    end

    return distance_matrix
end

function structure_factor(xs_plus, ys_plus, xs_minus, ys_minus, L, R; dr = 1)
    ks = reverse(2pi ./ collect(1:dr:L)) # fourier space 

    distance_matrices_all,
    distance_matrices_plus_minus,
    distance_matrices_plus_plus,
    distance_matrices_minus_minus = delta_matrix_defects(xs_plus, ys_plus, xs_minus, ys_minus, L, R)


    all = zeros(length(ks), R)
    plus_minus = zeros(length(ks), R)
    plus_plus = zeros(length(ks), R)
    minus_minus = zeros(length(ks), R)

    for r in 1:R
        all[:, r] = structure_factor(distance_matrices_all[r], ks)
        plus_minus[:, r] = structure_factor(distance_matrices_plus_minus[r], ks)
        plus_plus[:, r] = structure_factor(distance_matrices_plus_plus[r], ks)
        minus_minus[:, r] = structure_factor(distance_matrices_minus_minus[r], ks)
    end

    dico = Dict()
    dico["all"] = all
    dico["plus_minus"] = plus_minus
    dico["plus_plus"] = plus_plus
    dico["minus_minus"] = minus_minus
    dico["all_avg"] = mean(all, dims=2)[:,1]
    dico["plus_minus_avg"] = mean(plus_minus, dims=2)[:,1]
    dico["plus_plus_avg"] = mean(plus_plus, dims=2)[:,1]
    dico["minus_minus_avg"] = mean(minus_minus, dims=2)[:,1]

    return dico
end

function structure_factor(distance_matrix, ks)
    N = size(distance_matrix, 1)
    structure_factor = zeros(length(ks))
    for (i, k) in enumerate(ks)
        structure_factor[i] = abs(sum(exp.(-im * k * distance_matrix))) / N^2
        # structure_factor[i] = 1 + 2/N * sum(cos.(k * distance_matrix))
    end
    return structure_factor
end


function gr(xs_plus, ys_plus, xs_minus, ys_minus, L, R; dr = 1)
    distance_matrices_all,
    distance_matrices_plus_minus,
    distance_matrices_plus_plus,
    distance_matrices_minus_minus = distance_matrix_defects(xs_plus, ys_plus, xs_minus, ys_minus, L, R)

    rs = collect(1:dr:L/2)
    gr_all = zeros(length(rs), R)
    gr_plus_minus = zeros(length(rs), R)
    gr_plus_plus = zeros(length(rs), R)
    gr_minus_minus = zeros(length(rs), R)

    for r in 1:R
        gr_all[:, r] = gr(distance_matrices_all[r], rs, dr, L)
        gr_plus_minus[:, r] = gr(distance_matrices_plus_minus[r], rs, dr, L)
        gr_plus_plus[:, r] = gr(distance_matrices_plus_plus[r], rs, dr, L)
        gr_minus_minus[:, r] = gr(distance_matrices_minus_minus[r], rs, dr, L)
    end

    dico = Dict()
    dico["all"] = gr_all
    dico["plus_minus"] = gr_plus_minus
    dico["plus_plus"] = gr_plus_plus
    dico["minus_minus"] = gr_minus_minus
    dico["all_avg"] = mean(gr_all, dims=2)[:,1]
    dico["plus_minus_avg"] = mean(gr_plus_minus, dims=2)[:,1]
    dico["plus_plus_avg"] = mean(gr_plus_plus, dims=2)[:,1]
    dico["minus_minus_avg"] = mean(gr_minus_minus, dims=2)[:,1]

    return dico
end

function gr(distance_matrix, rs, dr, L)
    N = size(distance_matrix, 1)
    grr = zeros(length(rs))
    for (i, r) in enumerate(rs)
        indices = findall(r .<= distance_matrix .< r + dr)
        normalisation = N^2 / L^2 * 2pi * dr * r
        grr[i] = length(indices) / normalisation
    end
    return grr
end

function structure_factor_via_gr(xs_plus, ys_plus, xs_minus, ys_minus, L, R; dr = 1)
    rs = collect(1:dr:Lx/2)
    ks = reverse(2pi ./ rs)
    N = length(xs_plus)

    dico_GR = gr(xs_plus, ys_plus, xs_minus, ys_minus, L, R, dr = dr)
    GR_all_avg = dico_GR["all_avg"]
    GR_plus_minus_avg = dico_GR["plus_minus_avg"]
    GR_plus_plus_avg = dico_GR["plus_plus_avg"]
    GR_minus_minus_avg = dico_GR["minus_minus_avg"]

    SF_all_avg = zeros(length(ks))
    SF_plus_minus_avg = zeros(length(ks))
    SF_plus_plus_avg = zeros(length(ks))
    SF_minus_minus_avg = zeros(length(ks))

    for (i, k) in enumerate(ks)
        # S(k) = 1 + N/L^2 * integral_0^∞ of (g(r) - 1) * J_0(kr) * r * dr) 
        SF_all_avg[i] = 1 + N/L^2 * sum( (GR_all_avg .- 1) .* besselj0.(k * rs) .* rs * dr)
        SF_plus_minus_avg[i] = 1 + N/L^2 * sum( (GR_plus_minus_avg .- 1) .* besselj0.(k * rs) .* rs * dr)
        SF_plus_plus_avg[i] = 1 + N/L^2 * sum( (GR_plus_plus_avg .- 1) .* besselj0.(k * rs) .* rs * dr)
        SF_minus_minus_avg[i] = 1 + N/L^2 * sum( (GR_minus_minus_avg .- 1) .* besselj0.(k * rs) .* rs * dr)
    end

    dico = Dict()
    dico["all_avg"] = SF_all_avg
    dico["plus_minus_avg"] = SF_plus_minus_avg
    dico["plus_plus_avg"] = SF_plus_plus_avg
    dico["minus_minus_avg"] = SF_minus_minus_avg

    return dico
end
