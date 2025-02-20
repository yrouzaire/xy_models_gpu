#= 
for the parameters 
    Lx = 256
    Ly = 256
    R = 32
    
    lattice_type = "triangular"
    
    dt = 0.1
    tmax = 1000.0
    
    --------------- Successive improvements of the models (Square//Triangular) ---------------
                                               (base)    (wrt square)
    runtime XY         = 18.5 // 26 seconds |   100 %  |   140 % 
    runtime Kuramoto   = 18.9 // 25 seconds |   102 %  |   130 %
    runtime NRXY       = 21.6 // 33 seconds |   116 %  |   150 %
    runtime NRKuramoto = 22.1 // 33 seconds |   120 %  |   150 %
=#




include("../src/load_everything.jl")
Tf = Float32

T = Tf(0.01)
alpha = Tf(0.15) # non-reciprocal parameter (α=0 is fully reciprocal. The vision cone narrows as α increases) 
sigma = Tf(0.15) # standard deviation of the distribution of natural frequencies (ω) of the Kuramoto oscillators

Lx = 256
Ly = 256
R = 32

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)


init = "hightemp"
params_init = ("dummy") # only useful for the defect pair initialisation. Useless here. 

lattice_type = "triangular"
# lattice_type = "square"

dt = Tf(0.1)
tmax = Tf(1000.0)


## -------------------- XY model -------------------- ##
## -------------------- XY model -------------------- ##
## -------------------- XY model -------------------- ##
## -------------------- XY model -------------------- ##
thetas = create_thetas(Lx, Ly, R, init, params_init)
thetas_new = copy(thetas)

update_XY_triangular!(thetas, thetas_new, Lx, Ly, R, T, dt)
update_XY_square!(thetas, thetas_new, Lx, Ly, R, T, dt)


t = Float64(0.0)
z = @elapsed thetas, t = evolve_XY!(thetas, thetas_new, Lx, Ly, R, T, t, dt, tmax, lattice_type);
println("XY model:")
prinz(z)

## -------------------- NRXY model -------------------- ##
## -------------------- NRXY model -------------------- ##
## -------------------- NRXY model -------------------- ##
## -------------------- NRXY model -------------------- ##
thetas = create_thetas(Lx, Ly, R, init, params_init)
thetas_new = copy(thetas)

update_NRXY_triangular!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
update_NRXY_square!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)

t = Float64(0.0)
z = @elapsed thetas, t = evolve_NRXY!(thetas, thetas_new, Lx, Ly, R, T, sigma, t, dt, tmax, lattice_type);
println("NRXY model:")
prinz(z)

## -------------------- Kuramoto model -------------------- ##
## -------------------- Kuramoto model -------------------- ##
## -------------------- Kuramoto model -------------------- ##
## -------------------- Kuramoto model -------------------- ##
thetas = create_thetas(Lx, Ly, R, init, params_init)
thetas_new = copy(thetas)
omegas = create_omegas(Lx, Ly, R, alpha, "gaussian")

update_Kuramoto_triangular!(thetas, thetas_new, omegas, Lx, Ly, R, T, dt)
update_Kuramoto_square!(thetas, thetas_new, omegas, Lx, Ly, R, T, dt)

t = Float64(0.0)
z = @elapsed thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, tmax, lattice_type);
println("Kuramoto model:")
prinz(z)

## -------------------- Non-reciprocal Kuramoto model -------------------- ##
## -------------------- Non-reciprocal Kuramoto model -------------------- ##
## -------------------- Non-reciprocal Kuramoto model -------------------- ##
## -------------------- Non-reciprocal Kuramoto model -------------------- ##
thetas = create_thetas(Lx, Ly, R, init, params_init)
thetas_new = copy(thetas)
omegas = create_omegas(Lx, Ly, R, alpha, "gaussian")

update_NRKuramoto_triangular!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, dt)
update_NRKuramoto_square!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, dt)

t = Float64(0.0)
z = @elapsed thetas, t = evolve_NRKuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, sigma, t, dt, tmax, lattice_type);
println("NR Kuramoto model:")
prinz(z)


