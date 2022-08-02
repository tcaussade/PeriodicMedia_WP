using PeriodicMedia

setup = "2D1D"

if setup == "2D1D"
    θ = π/4.
    L = 2.0
    # k = [10.68, 20.]
    # k = [10.76, 20.]
    k = [2π/(L*(1-sin(θ))), 20.]
    P = Problem(k,θ,L; ambdim = 2, geodim = 1)
    Shape = PeriodicMedia.ParametricSurfaces.Kite
    # Shape = PeriodicMedia.ParametricSurfaces.Disk
    Fig = Obstacle(Shape,L/4)
elseif setup == "3D2D"
    k = [4.5, 8.0]
    θ = [π/4.,π/4.]
    L = [1.0, 1.0]
    P = Problem(k,θ,L; ambdim = 3, geodim = 2)
    Shape = PeriodicMedia.ParametricSurfaces.Sphere
    Fig = Obstacle(Sphere,minimum(L)/4)
end

# Mesh parameters
ppw = 8
dim = 5

# Experiment
λ = 2π/k[1]
windowsizes = λ * collect(10:1:30)
ef = Vector{Float64}(undef,0)
et = Vector{Float64}(undef,0)
for w in windowsizes
    println("Solving with A = "*string(round(w*k[1]/2π,digits=1))*"λ")
    WGF = Window(0.5, w)
    Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
    Γt = extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
    # Non-corrected solution
    ϕf = solver(P,Γs,Γt,WGF; FRO = false)
    @show ebf = energytest(P,Γt,WGF, ϕf; FRO = false, H = 1.0)
    push!(ef,ebf)
    # Corrected solution
    ϕt = solver(P,Γs,Γt,WGF; FRO = true)
    @show ebt = energytest(P,Γt,WGF, ϕt; FRO = true, H = 1.0)
    push!(et,ebt)  
end

import Plots
lbs     = ["without correction" "with correction"]
semilog = Plots.plot(windowsizes/λ, log10.([ef et]); title = "Semilog", label = lbs)
loglog  = Plots.plot(log10.(windowsizes/λ), log10.([ef et]), title = "Log-Log", label = lbs)
Plots.plot(semilog,loglog, legend = true, xlabel = "A/λ", ylabel = "EB", ylims = (-6,-1))