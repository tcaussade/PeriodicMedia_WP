using PeriodicMedia

# set physical params
k = [10.68, 20.]
θ = π/4.
L = 2.0
P = Problem(k,θ,L; ambdim = 2, geodim = 1)

# Set Windowed Green function parameters
WGF = Window(0.5,10*(2π/k[1]))

#create geometry
PeriodicMedia.clear_entities!()
Kite = PeriodicMedia.ParametricSurfaces.Kite
Disk = PeriodicMedia.ParametricSurfaces.Disk
Fig = Obstacle(Disk,L/4)

ppw = 3
dim = 2

Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
Γt = extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)

fro = true
# Find the unkwnown densities
ϕ = solver(P,Γs,Γt,WGF; FRO = fro)
# Asses method accuracy (note it uses triple cell configuration)
@show energytest(P,Γt,WGF, ϕ; FRO = fro, H=1.0)

# Plot the solution over the desired cells
X,Y,U = cellsolution(P,Γt,WGF,ϕ; ppw = 20, FRO = fro)
viewsolution(P,X,Y,-U,Fig; ncell = -1:1)

