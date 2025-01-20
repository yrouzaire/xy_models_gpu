include("../src/load_everything.jl")
Tf = Float32

T = Tf(0.01)

Lx = 256
Ly = 256
R = 16

wrapsT = 16 # tune this parameter to optimize the runtime
block3D = (wrapsT, wrapsT, 1)
grid3D = (Int(ceil(Lx / wrapsT)), Int(ceil(Ly / wrapsT)), R)


init = "hightemp"
lattice_type = "triangular"

dt = Tf(0.1)
tmax = Tf(100.0)

thetas = create_thetas(Lx, Ly, R, init)
thetas_new = copy(thetas)

update_XY_square!(thetas, thetas_new, Lx, Ly, R, T, dt)
update_XY_triangular!(thetas, thetas_new, Lx, Ly, R, T, dt)

t = Float64(0.0)
CUDA.@time thetas, t = evolve_XY!(thetas, thetas_new, Lx, Ly, R, T, t, dt, tmax, lattice_type)






