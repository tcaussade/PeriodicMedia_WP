# /home/thomas/Downloads/julia-1.7.2/bin/julia
# import Pkg; Pkg.activate("/home/thomas/PeriodicMedia_WP")
using PeriodicMedia

# set physical params
# k1 = 4.59961087822572 #RW-anomaly
k1 = 9.0 # only m = (0,0) propagates
# k = [k1, 8.0]
k = [k1,15.]
θ = [π/4.,π/4.]
# θ = [π/4.,π/2.]
L = [1.0, 1.0]
P = Problem(k,θ,L; ambdim = 3, geodim = 2)

# Set Windowed Green function parameters
WGF = Window(0.6,15*(2π/k[1]))

#create geometry
PeriodicMedia.clear_entities!()
Sphere = PeriodicMedia.ParametricSurfaces.Sphere
Fig = Obstacle(Sphere,minimum(L)/4)

ppw = 4
dim = 2

Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
@show smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2+length(Γs[4].dofs)*2
Γt = extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)

PeriodicMedia.meshplot(Γs[1])
PeriodicMedia.meshplot(get(Γt[2],-1,""))
# PeriodicMedia.meshplot!(get(Γt[2],0,""))
# PeriodicMedia.meshplot!(get(Γt[5],-1,""))
# PeriodicMedia.meshplot!(get(Γt[5],+1,""))

# PeriodicMedia.meshplot!(Γs[1])

# PeriodicMedia.meshplot(get(Γt[4],0,""))

# Find the unkwnown densities
# @time ϕt = solver(P,Γs,Γt,WGF; FRO = true)
@time ϕf = solver(P,Γs,Γt,WGF; FRO = false)

# Asses method accuracy (note it uses triple cell configuration)
# @show energytest(P,Γt,WGF, ϕt; FRO = true, H = 1.0)
@show eb = energytest(P,Γt,WGF, ϕf; FRO = false, H = 1.0)

# Plot the solution over the desired cells
# X,Y,Z, U = cellsolution(P,Γt,WGF,ϕt; ppw = 20, zlims = [-2.0,2.0], FRO = true)
X,Y,Z, U = cellsolution(P,Γt,WGF,ϕf; ppw = 20, zlims = [-2.0,4.0], FRO = false)

import Plots
ncell = -1:0
p1 = XYviewsolution(P,X,Y,U.XY; ncell = ncell)
p2 = YZviewsolution(P,Y,Z,U.YZ; ncell = ncell)
p3 = XZviewsolution(P,X,Z,U.XZ; ncell = ncell)
λ  = 2π/k[1]
lb = "(c,A/λ)="*string((WGF.c,WGF.A/λ))*", k= "*string(k)*", θ= "*string(θ./π)*"π, L="*string(L)
println( lb )
println( "EB: "*string(round(eb, digits = 8)) )
p = Plots.plot(p1,p2,p3, layout=(1,3), size = (1500,700), plot_title = lb*" EB: "*string(round(eb, digits = 8)))
#Plots.savefig(p,"bigexperiment.png")
