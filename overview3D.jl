# /home/thomas/Downloads/julia-1.7.2/bin/julia
# import Pkg; Pkg.activate("/home/thomas/PeriodicMedia_WP")

using PeriodicMedia

correct  = false
test     = true
plotting = true

# set physical params
θ   = [π/6.,π/4.] 
L   = [0.5, 0.5]
k1  = 9.2 
k2  = 15.0
pol = "TE"

P = Problem([k1,k2],θ,L,pol; ambdim = 3, geodim = 2)

# set discretization parameters
λ    = 2π/k1
c    = 0.5
A    = 5*λ
Wpar = Window(c,A)

ppw = 3
dim = 1

#create geometry
Shape = PeriodicMedia.ParametricSurfaces.Sphere
PeriodicMedia.clear_entities!()

Fig = Obstacle(Shape,minimum(L)/4)
Γs = unitcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
Γt = extendedcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
G  = (Γs=Γs, Γt=Γt)
smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2+length(Γs[4].dofs)*2

@info "Geometry created" Wpar c*A smat
@info "Assembling..."

# Solving for the scattered field
MB,Wa,E = matrixcreator(P,G,Wpar)
b       = rightside(P,Γs)

δ       = 0.75*k1
hcorr   = Wpar.c * Wpar.A
if correct  
    @info "Adding corrections..." δ hcorr
    @assert Γt[1] isa Dict{Tuple{Int,Int}}
    MB += finiterankoperator(P,Γs,Γt; δ = δ , h = hcorr)
end

# GMRES parameters
res = size(MB,2)
vbs = false
tol = 1e-6
# Solve linear system
ϕ,niter = gmres(E+MB*Wa,b; restart = res, verbose = vbs, reltol = tol, log = true)
@info niter

he = 0.9 * hcorr
eb = 0.0
if test  
    @info "Computing energy" he
    @assert he < Wpar.c * Wpar.A
    R,T = energytest(P,Γt,Wpar, ϕ; FRO = correct, h = he)
    eb = abs(R+T-1)
    @info "Test results" eb R T 
end

import Plots
include("matlab_export.jl")
p = Plots.Plot()
if plotting
    
    ncell = -3:3
    len   = length(ncell) 
    zlims = [-3*maximum(L),3*maximum(L)]
    res   = 120
    @show h = 2π/max(P.pde[1].k,P.pde[2].k) / res
    @info "Plotting solution" len res
    X,Y,Z, U = cellsolution(P,Γt,Wpar,ϕ; ppw = res, zlims = zlims, FRO = correct)
    # p1 = XYviewsolution(P,X,Y,U.XY; ncell = ncell)
    # p2 = YZviewsolution(P,Y,Z,U.YZ; ncell = ncell)
    # p3 = XZviewsolution(P,X,Z,U.XZ; ncell = ncell)

    # lb = "(c,A/λ)="*string((Wpar.c,Wpar.A/λ))*", k= "*string([k1 k2])*", θ= "*string(θ./π)*"π, L="*string(L)
    # if test
    #     lb *=" EB: "*string(round(eb, digits = 8))
    # end
    # p = Plots.plot(p1,p2,p3, layout=(1,3), size = (1500,700), plot_title = lb)
    # # Plots.savefig(p,"bigexperiment.png")

    matlab_export(X,Y,Z,U,"3d_data")
end


# # Find the unkwnown densities
# @time ϕt = solver(P,Γs,Γt,WGF; FRO = true)
# @time ϕf = solver(P,Γs,Γt,WGF; FRO = false)

# # Asses method accuracy (note it uses triple cell configuration)
# @show energytest(P,Γt,WGF, ϕt; FRO = true, H = 1.0)
# @show eb = energytest(P,Γt,WGF, ϕf; FRO = false, H = 1.0)

# # Plot the solution over the desired cells
# X,Y,Z, U = cellsolution(P,Γt,WGF,ϕt; ppw = 20, zlims = [-2.0,2.0], FRO = true)
# X,Y,Z, U = cellsolution(P,Γt,WGF,ϕf; ppw = 20, zlims = [-2.0,4.0], FRO = false)

# import Plots
# ncell = -1:0
# p1 = XYviewsolution(P,X,Y,U.XY; ncell = ncell)
# p2 = YZviewsolution(P,Y,Z,U.YZ; ncell = ncell)
# p3 = XZviewsolution(P,X,Z,U.XZ; ncell = ncell)
# λ  = 2π/k[1]
# lb = "(c,A/λ)="*string((WGF.c,WGF.A/λ))*", k= "*string(k)*", θ= "*string(θ./π)*"π, L="*string(L)
# println( lb )
# println( "EB: "*string(round(eb, digits = 8)) )
# p = Plots.plot(p1,p2,p3, layout=(1,3), size = (1500,700), plot_title = lb*" EB: "*string(round(eb, digits = 8)))
# #Plots.savefig(p,"bigexperiment.png")



# PeriodicMedia.meshplot(Γs[1])
# PeriodicMedia.meshplot(get(Γt[2],-1,""))
# PeriodicMedia.meshplot!(get(Γt[2],0,""))
# PeriodicMedia.meshplot!(get(Γt[5],-1,""))
# PeriodicMedia.meshplot!(get(Γt[5],+1,""))
# PeriodicMedia.meshplot!(Γs[1])
# PeriodicMedia.meshplot(get(Γt[4],0,""))

