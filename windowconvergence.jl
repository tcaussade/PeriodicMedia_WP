using PeriodicMedia
import Plots

function windowconvergence(P::Problem, Fig::Obstacle, windowsizes)
    ef = Vector{Float64}(undef,0)
    et = Vector{Float64}(undef,0)
    Threads.@threads for w in windowsizes
        println("Solving with A = "*string(round(w/λ, digits=1))*"λ")
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
    return ef,et
end
function viewresults(windowsizes,ef,et; semilog)
    lbs     = ["without correction" "with correction"]
    if semilog
        return Plots.plot(windowsizes/λ, log10.([ef et]); title = "k="*string(round(2π/λ,digits = 2)), label = lbs)
    else
        return Plots.plot(log10.(windowsizes/λ), log10.([ef et]), "k="*string(round(2π/λ,digits = 2)), label = lbs)
    end
    # Plots.plot(semilog,loglog, legend = true, xlabel = "A/λ", ylabel = "EB", ylims = (-6,-1))
end

""" 
    set desired experiment
    - setup = "2D1D"
    - setup = "3D2D"
"""
setup = "2D1D"

if setup == "2D1D"
    θ = π/4.
    L = 2.0
    k1 = [10.68, 2π/(L*(1-sin(θ))), 10.76] # k2 = 20.
    # Shape = PeriodicMedia.ParametricSurfaces.Kite
    Shape = PeriodicMedia.ParametricSurfaces.Disk
    Fig = Obstacle(Shape,L/4)

    ####################################
    global ppw = 8
    global dim = 5
    ####################################

    p = Vector{Plots.Plot}(undef,3)
    Threads.@threads for i = 1:3
        P = Problem([k1[i],20.],θ,L; ambdim = 2, geodim = 1)

        ####################################
        global λ = 2π/k1[i]
        windowsizes = λ * collect(10:5:40)
        ####################################

        ef,et = windowconvergence(P,Fig,windowsizes)
        p[i] = viewresults(windowsizes,ef,et; semilog = true) 
    end 
    Plots.savefig(  Plots.plot(p[1],p[2],p[3], layout = Plots.@layout([A B C]), 
                    xlabel = "A/λ", ylabel = "EB", # ylims = (-6,-1),
                    size = (1500,800)), "2d1d.png")

elseif setup == "3D2D"
    k = [4.5, 8.0]
    θ = [π/4.,π/4.]
    L = [1.0, 1.0]
    P = Problem(k,θ,L; ambdim = 3, geodim = 2)
    Shape = PeriodicMedia.ParametricSurfaces.Sphere
    Fig = Obstacle(Sphere,minimum(L)/4)
    ef,et = windowconvergence(P,Fig)
end


