include("../src/load_everything.jl")
Tf = Float32

## ----------------- Test Init Pair ----------------- ##
## ----------------- Test Init Pair ----------------- ##
## ----------------- Test Init Pair ----------------- ##
## ----------------- Test Init Pair ----------------- ##
T = Tf(0.01)
Lx = 512
Ly = 512
R = 1

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)

init = "hightemp"
params_init = (x_plus=Lx/4, y_plus=Ly/2, r0=Lx / 2, mu_plus=0.0)
lattice_type = "square"

dt = Tf(0.1)
tmax = Tf(100)

thetas = create_thetas(Lx, Ly, R, init , params_init)
thetas_new = copy(thetas)
omegas = create_omegas(Lx, Ly, R, 0.1, "normal")
# update_XY_square!(thetas, thetas_new, Lx, Ly, R, T, dt)
# update_XY_triangular!(thetas, thetas_new, Lx, Ly, R, T, dt)

t = Float64(0.0)
thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, 5000.0, lattice_type)
# plot_thetas(circshift(Array(thetas), (0,0Ly/2, 0)))
plot_thetas(thetas, defects=true, lattice_type="square")
number_defects_plus_minus(thetas, lattice_type)

## ----------------- Kuramoto Model ----------------- ##
## ----------------- Kuramoto Model ----------------- ##
## ----------------- Kuramoto Model ----------------- ##
## ----------------- Kuramoto Model ----------------- ##
T = Tf(0.01)
sigma = 0.15

Lx = 512
Ly = 512
R = 2

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)


init = "hightemp"
lattice_type = "triangular"

dt = Tf(0.1)
tmax = Tf(2000.0)

thetas = create_thetas(Lx, Ly, R, init)
omegas = create_omegas(Lx, Ly, R, sigma, "gaussian")
thetas_new = copy(thetas)

update_Kuramoto_square!(thetas, thetas_new, omegas, Lx, Ly, R, T, dt)
update_Kuramoto_triangular!(thetas, thetas_new, omegas, Lx, Ly, R, T, dt)

t = Float64(0.0)
thetas, t = evolve_Kuramoto!(thetas, thetas_new, omegas, Lx, Ly, R, T, t, dt, tmax, lattice_type)
plot_thetas(thetas)






## ----------------- NRYY Model ----------------- ##
## ----------------- NRYY Model ----------------- ##
## ----------------- NRYY Model ----------------- ##
## ----------------- NRYY Model ----------------- ##
T = Tf(0.01)
sigma = Tf(0.15)

Lx = 128
Ly = 128
R = 2

wrapsT = 16
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)


init = "pair"
lattice_type = "triangular"
# lattice_type = "square"

dt = Tf(0.1)
tmax = Tf(1000.0)

thetas = create_thetas(Lx, Ly, R, init)
# thetas_new = copy(thetas)

# update_NRXY_square!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)
# update_NRXY_triangular!(thetas, thetas_new, Lx, Ly, R, T, sigma, dt)

# t = Float64(0.0)
# thetas, t = evolve_NRXY!(thetas, thetas_new, Lx, Ly, R, T, sigma, t, dt, tmax, lattice_type)
plot_thetas(thetas)


