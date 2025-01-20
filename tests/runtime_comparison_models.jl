#= 
for the parameters 
    Lx = 256
    Ly = 256
    R = 16
    
    lattice_type = "triangular"
    
    dt = 0.1
    tmax = 1000.0
    
    --------------- Successive improvements of the models (Square//Triangular) ---------------

    runtime XY       = 10//14 seconds |
    runtime NRXY     = 11//17 seconds |
    runtime Kuramoto = 10//13 seconds |
=#




include("../src/load_everything.jl")
Tf = Float32

T = Tf(0.01)
sigma = Tf(0.15)

Lx = 256
Ly = 256
R = 16

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)


init = "hightemp"
lattice_type = "triangular"
lattice_type = "square"

dt = Tf(0.1)
tmax = Tf(1000.0)


## -------------------- XY model -------------------- ##
## -------------------- XY model -------------------- ##
## -------------------- XY model -------------------- ##
## -------------------- XY model -------------------- ##
t = Float64(0.0)
thetas = create_thetas(Lx, Ly, R, init)
thetas_new = copy(thetas)

update_XY_triangular!(thetas, thetas_new, Lx, Ly, R, T, dt)
update_XY_square!(thetas, thetas_new, Lx, Ly, R, T, dt)


t = Float64(0.0)
z = @elapsed thetas, t = evolve_XY!(thetas, thetas_new, Lx, Ly, R, T, t, dt, tmax, lattice_type);
prinz(z)

## -------------------- NRXY model -------------------- ##
## -------------------- NRXY model -------------------- ##
## -------------------- NRXY model -------------------- ##
## -------------------- NRXY model -------------------- ##
t = Float64(0.0)
thetas = create_thetas(Lx, Ly, R, init)
thetas_new = copy(thetas)

update_NRXY_triangular!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
update_NRXY_square!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)

t = Float64(0.0)
z = @elapsed thetas, t = evolve_NRXY!(thetas, thetas_new, Lx, Ly, R, T, sigma, t, dt, tmax, lattice_type);
prinz(z)

## -------------------- Kuramoto model -------------------- ##
## -------------------- Kuramoto model -------------------- ##
## -------------------- Kuramoto model -------------------- ##
## -------------------- Kuramoto model -------------------- ##
t = Float64(0.0)
thetas = create_thetas(Lx, Ly, R, init)
thetas_new = copy(thetas)
omegas = create_omegas(Lx, Ly, R, sigma, "gaussian")

update_Kuramoto_triangular!(thetas, thetas_new, omegas, Lx, Ly, R, T, dt)
update_Kuramoto_square!(thetas, thetas_new, omegas, Lx, Ly, R, T, dt)

t = Float64(0.0)
z = @elapsed thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, tmax, lattice_type);
prinz(z)


