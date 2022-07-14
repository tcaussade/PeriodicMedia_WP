using devPeriodicMedia

# set physical params
k = [4.5, 8.]
θ = [π/12.,0.]
L = [1.0, 1.0]
P = Problem(k,θ,L; ambdim = 3, geodim = 2)

# Set Windowed Green function parameters
WGF = Window(0.7,5*(2π/k[1]))

#create geometry
devPeriodicMedia.clear_entities!()
Sphere = devPeriodicMedia.ParametricSurfaces.Sphere
Fig = devPeriodicMedia.Obstacle(Sphere,maximum(L)/4)

ppw = 3
dim = 1

Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
Γt = devPeriodicMedia.extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)

# Find the unkwnown densities
ϕ = solver(P,Γs,Γt,WGF)

# Asses method accuracy (note it uses triple cell configuration)
@show devPeriodicMedia.energytest(P,Γt,WGF, ϕ; FRO = false)

# Plot the solution over the desired cells
X,Y,U = devPeriodicMedia.cellsolution(P,Γt,WGF,ϕ; ppw = 20)
devPeriodicMedia.viewsolution(P,X,Y,U,Fig; ncell = -2:2)




