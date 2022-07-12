using devPeriodicMedia

# set physical params
k = [10.68, 20.]
θ = π/4.
L = 2.0
P = Problem(k,θ,L; ambdim = 2, geodim = 1)

# Set Windowed Green function parameters
WGF = Window(0.5,15*(2π/k[1]))

#create geometry
Kite = devPeriodicMedia.ParametricSurfaces.Kite
Disk = devPeriodicMedia.ParametricSurfaces.Disk
Fig = devPeriodicMedia.Obstacle(Kite,L/4)

ppw = 4
dim = 2

Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
Γt = devPeriodicMedia.triplecell(P,Fig, WGF; ppw = ppw, dimorder = dim)

# Find the unkwnown densoties
ϕ = solver(P,Γs,Γt,WGF; FRO = true)

# Asses method accuracy (note it uses triple cell configuration)
@show devPeriodicMedia.energytest(P,Γt,WGF, ϕ)

# Plot the solution over the desired cells
X,Y,U = devPeriodicMedia.cellsolution(P,Γt,WGF,ϕ; ppw = 20)
devPeriodicMedia.viewsolution(P,X,Y,U,Fig; ncell = -1:1)