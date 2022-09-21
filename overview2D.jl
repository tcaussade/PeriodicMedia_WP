using PeriodicMedia

correct  = false
test     = true
plotting = false

# set physical parameters
θ  = π/4.
L  = 2.0
k1 = 10.76 # 10.68
k2 = 20.0
pol = "TE"

P = Problem([k1,k2],θ,L,pol; ambdim = 2, geodim = 1)

# set discretization parameters
λ    = 2π/k1
c    = 0.5
A    = 40*λ
Wpar = Window(c,A)

ppw = 7
dim = 5

# create geometry

Shape = PeriodicMedia.ParametricSurfaces.Kite
# Shape = PeriodicMedia.ParametricSurfaces.Disk
PeriodicMedia.clear_entities!()
Fig = Obstacle(Shape,L/4)
Γs = unitcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
Γt = extendedcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
G  = (Γs=Γs, Γt=Γt)

@info "Geometry created" Wpar c*A
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
    ϕ = (E+MB*Wa)\b
else
    # Solving without corrections
    ϕ = (E+MB*Wa)\b
end

he = 0.9 * hcorr
if test  
    @info "Computing energy" he
    R,T = energytest(P,Γt,Wpar, ϕ; FRO = correct, h = he)
    err = abs(R+T-1)
    @info "Test results" err R T 
end

if plotting
    ncell = -1:1
    len   = length(ncell)
    resolution = 20
    @info "Plotting solution" len resolution
    X,Y,U = cellsolution(P,Γt,Wpar,ϕ; ppw = resolution, FRO = correct)
    viewsolution(P,X,Y,U,Fig; ncell = ncell)
end

# import Plots
# Plots.heatmap(X,Y,log10.(abs.(Ut-U)), aspect_ratio = 1, clim = (-6, 0))
# Plots.plot(log10.(abs.(ϕ-ϕt)))