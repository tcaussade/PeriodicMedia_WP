using PeriodicMedia

# set physical params
k = [10.68, 20.]
θ = π/4.
L = 2.0
P = Problem(k,θ,L; ambdim = 2, geodim = 1)

# Set Windowed Green function parameters
WGF = Window(0.5,15*(2π/k[1]))

#create geometry
PeriodicMedia.clear_entities!()
Kite = PeriodicMedia.ParametricSurfaces.Kite
Disk = PeriodicMedia.ParametricSurfaces.Disk
Fig = Obstacle(Kite,L/4)

ppw = 4
dim = 2

Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
Γt = extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)

# Find the unkwnown densoties
ϕ = solver(P,Γs,Γt,WGF; FRO = false)

# Asses method accuracy (note it uses triple cell configuration)
@show energytest(P,Γt,WGF, ϕ; FRO = false)

# Plot the solution over the desired cells
X,Y,U = cellsolution(P,Γt,WGF,ϕ; ppw = 20)
viewsolution(P,X,Y,U,Fig; ncell = -0:0)