using PeriodicMedia

correct  = false
test     = false
plotting = true

# set physical parameters
θ  = [π/4.,π/4.]
L  = 1.0
k1 = 10.68 # 10.76
k2 = 20.0
pol = "TE"

P = Problem([k1,k2],θ,L,pol; ambdim = 3, geodim = 1)

# set discretization parameters
λ    = 2π/k1
c    = 0.5
A    = 5*λ
Wpar = Window(c,A)

ppw = 2
dim = 2

# create geometry

Shape = PeriodicMedia.ParametricSurfaces.Sphere
PeriodicMedia.clear_entities!()
Fig = Obstacle(Shape,L/4)
Γs = unitcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
Γt = extendedcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
G  = (Γs=Γs, Γt=Γt)
smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2

# PeriodicMedia.meshplot(Γs[1])
# PeriodicMedia.meshplot!(Γt[2])
# PeriodicMedia.meshplot!(Γt[3])
# PeriodicMedia.meshplot(get(Γt[1],-1,""))


@info "Geometry created" Wpar c*A smat
@info "Assembling..."

# Solving for the scattered field
MB,Wa,E = matrixcreator(P,G,Wpar)
b       = rightside(P,Γs)

δ       = 0.75*k1
hcorr   = Wpar.c * Wpar.A
if correct  
    @info "Adding corrections..." δ hcorr
    @assert Γt[1] isa Dict
    MB += finiterankoperator(P,Γs,Γt; δ = δ , h = hcorr)
    # Solving with corrections
else
    # Solving without corrections
end

# GMRES parameters
res = size(MB,2)
vbs = false
tol = 1e-6
# Solve linear system
ϕ,niter = gmres(E+MB*Wa,b; restart = res, verbose = vbs, reltol = tol, log = true)
@info niter


he = 0.9 * hcorr
if test  
    @info "Computing energy" he
    e,R,T = energytest(P,Γt,Wpar, ϕ; FRO = correct, h = he)
    # err = abs(R+T-1)
    @info "Test results" e R T 
end


if plotting
    import Plots
    ncell = -1:1
    len   = length(ncell)
    resolution = 10
    @info "Plotting solution" len resolution
    X,Y,Z,U = cellsolution(P,Γt,Wpar,ϕ; ppw = resolution, FRO = correct)
    p1 = XYviewsolution(P,X,Y,U.XY; ncell = ncell)
    p2 = YZviewsolution(P,Y,Z,U.YZ; ncell = ncell)
    p3 = XZviewsolution(P,X,Z,U.XZ; ncell = ncell)
    Plots.plot(p1,p2,p3, layout=(1,3), size = (1500,700), plot_title = "3D1D")
end

# import Plots
# Plots.heatmap(X,Y,log10.(abs.(Ut-U)), aspect_ratio = 1, clim = (-6, 0))
# Plots.plot(log10.(abs.(ϕ-ϕt)))