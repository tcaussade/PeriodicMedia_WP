using PeriodicMedia

# set physical params
θ = π/4.
L = 2.0
# k = [10.68, 20.]
# k = [10.76, 20.]
k = [2π/(L*(1-sin(θ))), 20.]
P = Problem(k,θ,L; ambdim = 2, geodim = 1)

# Set Windowed Green function parameters
WGF = Window(0.5,30*(2π/k[1]))

#create geometry
PeriodicMedia.clear_entities!()
Kite = PeriodicMedia.ParametricSurfaces.Kite
Disk = PeriodicMedia.ParametricSurfaces.Disk
Fig = Obstacle(Kite,L/4)

ppw = 8
dim = 5

Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
Γt = extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)

fro = true
# Find the unkwnown densities
ϕt = solver(P,Γs,Γt,WGF; FRO = fro)
ϕ = solver(P,Γs,Γt,WGF; FRO = ~fro)
# Asses method accuracy (note it uses triple cell configuration)
@show energytest(P,Γt,WGF, ϕt; FRO = fro, H=1.0)
@show energytest(P,Γt,WGF, ϕ; FRO = ~fro, H=1.0)

# Plot the solution over the desired cells
X,Y,U = cellsolution(P,Γt,WGF,ϕ; ppw = 20, FRO = fro)
viewsolution(P,X,Y,-U,Fig; ncell = -1:1)

import Plots
Plots.heatmap(X,Y,log10.(abs.(Ut-U)), aspect_ratio = 1, clim = (-6, 0))
Plots.plot(log10.(abs.(ϕ-ϕt)))