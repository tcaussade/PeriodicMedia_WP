using devPeriodicMedia

# set physical params
k = [4.5, 8.]
θ = [π/4.,π/4.]
L = [1.0, 1.0]
P = Problem(k,θ,L; ambdim = 3, geodim = 2)

# Set Windowed Green function parameters
WGF = Window(0.7,5*(2π/k[1]))

#create geometry
devPeriodicMedia.clear_entities!()
Sphere = devPeriodicMedia.ParametricSurfaces.Sphere
Fig = devPeriodicMedia.Obstacle(Sphere,maximum(L)/4)

ppw = 5
dim = 3

Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2+length(Γs[4].dofs)*2
Γt = devPeriodicMedia.extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)

# Find the unkwnown densities
@time ϕ = solver(P,Γs,Γt,WGF)

# Asses method accuracy (note it uses triple cell configuration)
@show devPeriodicMedia.energytest(P,Γt,WGF, ϕ; FRO = false)

# Plot the solution over the desired cells
X,Y,Z, U = devPeriodicMedia.cellsolution(P,Γt,WGF,ϕ; ppw = 20)

import Plots
ncell = -1:1
p1 = devPeriodicMedia.XYviewsolution(P,X,Y,U.XY; ncell = ncell)
p2 = devPeriodicMedia.YZviewsolution(P,Y,Z,U.YZ; ncell = ncell)
p3 = devPeriodicMedia.XZviewsolution(P,X,Z,U.XZ; ncell = ncell)
Plots.plot(p1,p2,p3, layout=(1,3), size = (1500,700) )





