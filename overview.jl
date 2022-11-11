using PeriodicMedia

""" 
    ambient_dimension  : spatial dimension Rᵈ for d=2,3
    geometric_dimension: d=1 for line-arrays or d=2 for surface-arrays
"""

ambient_dimension   = 2
geometric_dimension = 1

"""
    correct = add finite-rank corrections
    etest   = compute the energy balance, reflectance and transmittance
    plots   = show scattered field plot
"""

correct = true
etest   = true
plots   = true

### Physical parameters ###
k1  = 8.8    # exterior wavenumber
k2  = 14.0   # interior wavenumber
pol = "TE"   # polarization

λ   = 2π/k1  # wavelength

""" 
    Relevant parameters
    - θ : incidence angle
    - L : periodic array spatial period

    Obstacle parameters
    - Shape : shape
    - r     : radius
"""

if ambient_dimension == 2 && geometric_dimension == 1
    θ  = π/4.
    L  = 2.0
    Shape = PeriodicMedia.ParametricSurfaces.Kite
    r     = L/4

elseif ambient_dimension == 3 && geometric_dimension == 1
    θ  = [π/4.,π/4.]
    L  = 1.0
    Shape = PeriodicMedia.ParametricSurfaces.Sphere
    r     = L/4

elseif ambient_dimension == 3 && geometric_dimension == 2
    θ   = [π/6.,π/4.] 
    L   = [0.5, 0.5]
    Shape = PeriodicMedia.ParametricSurfaces.Sphere
    r = minimum(L)/4

end

# stores the physical data
P = Problem([k1,k2],θ,L,pol; ambdim = ambient_dimension, geodim = geometric_dimension)

### Discretization parameters ###
c    = 0.5 # parameter that controls the window function decayment
A    = 10λ  # window sisze
Wpar = Window(c,A) # stores relevant windowing parameters

ppw = 7 # Points per wavelength to build meshes
dim = 3 # element integration order

δ       = 0.75*k1 # parameter to select the correction modes in Cδ    
hcorr   = c*A     # height used to compute corrective terms
he      = 0.9 * hcorr # height used to compute the energy test
@assert r < he 

### Geometry ###
PeriodicMedia.clear_entities!()

"""
    This version supports currently only two geometry hyperparameters: shape and radius.
    If more hyper-parameters were desired such as rotation, number of layers, etc 
    it is encouraged to add such functionalities into the Obstacle struct in "src/Configuration.jl"
"""

Fig = Obstacle(Shape, r)                                  # stores obstacle data
Γs = unitcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)     # single-cell mesh construction
Γt = extendedcell(P,Fig, Wpar; ppw = ppw, dimorder = dim) # extended-cell mesh construction
G  = (Γs=Γs, Γt=Γt)                                       # stores meshed geometry

@info "Geometry created"

### Linear system matrices ###
MB, Wa, E = matrixcreator(P,G,Wpar) # main system matrices
b         = rightside(P,Γs)         # right-hand side of the system

if correct
    if ambient_dimension == 3
        @warn "Corrections are not (yet) validated, disabling corrections" 
    else
        @info "Adding corrections..." δ hcorr
        MB += finiterankoperator(P,Γs,Γt; δ = δ , h = hcorr)
    end
end

### Linear system solution ###

"""
    We use GMRES iterative solver if the problem is in R³
    We use gaussian elimination exact solver if the problem is in R²
"""

if ambient_dimension == 3
    res = size(MB,2) # GMRES iterations restart parameters
    vbs = false      # verbose
    tol = 1e-6       # GMRES tolerance
    ϕ,niter = gmres(E+MB*Wa,b; restart = res, verbose = vbs, reltol = tol, log = true)
    @info niter
elseif ambient_dimension == 2
    ϕ = (E+MB*Wa)\b
end

### Energy conservation error ###
if etest
    if ambient_dimension == 3
        @warn "Energy balance is not (yet) validated"
    end
    if ambient_dimension-geometric_dimension == 1
        @info "Computing energy..." he
        eb,R,T = energytest(P,Γt,Wpar, ϕ; FRO = correct, h = he)
        @info "Test results:" eb R T 
    else
        @warn "Not (yet) implemented"
        # @info "Computing energy..." he
        # eb = energytest(P,Γt,Wpar, ϕ; FRO = correct, h = he)
        # @info "Test results:" eb
    end
end

### Plotting the solution ###

"""
    1. Evaulate potential for a single cell
        - cellsolution()
    2. Retrieve over a range "ncell" using quasi-periodicity
        - viewsolution()
"""

ncell      = -1:1 # range of cells to be plotted
resolution = 10   # points per wavelength used for plotting

if plots
    if ambient_dimension == 2
        X,Y,U = cellsolution(P,Γt,Wpar,ϕ; ppw = resolution, FRO = correct)
        field = viewsolution(P,X,Y,U,Fig; ncell = ncell)
    elseif ambient_dimension == 3
        X,Y,Z,U = cellsolution(P,Γt,Wpar,ϕ; ppw = resolution, FRO = false)
        p1      = XYviewsolution(P,X,Y,U.XY; ncell = ncell)
        p2      = YZviewsolution(P,Y,Z,U.YZ; ncell = ncell)
        p3      = XZviewsolution(P,X,Z,U.XZ; ncell = ncell)
        field   = PeriodicMedia.Plots.plot(p1,p2,p3, layout=(1,3), size = (1500,700))
    end
    field
end
