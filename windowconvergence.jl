# /home/thomas/Downloads/julia-1.7.2/bin/julia
# import Pkg; Pkg.activate("/home/thomas/PeriodicMedia_WP")

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
        # ϕt = solver(P,Γs,Γt,WGF; FRO = true)
        # @show ebt = energytest(P,Γt,WGF, ϕt; FRO = true, H = 1.0)
        ebt = 1.0
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
setup = "3D2D"
p = Vector{Plots.Plot}(undef,3)

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
    k1 = [9.0, 9.199221756451442, 9.5]
    θ = [π/4.,π/4.]
    L = [0.5, 0.5]
    Shape = PeriodicMedia.ParametricSurfaces.Sphere
    Fig = Obstacle(Shape,minimum(L)/4)

    ####################################
    global ppw = 1
    global dim = 1
    ####################################

    
    Threads.@threads for i = 1:3
        P = Problem([k1[i],17.],θ,L; ambdim = 3, geodim = 2)

        ####################################
        global λ = 2π/k1[i]
        windowsizes = λ * collect(10:10:30)
        ####################################

        ef,et = windowconvergence(P,Fig,windowsizes)
        p[i] = viewresults(windowsizes,ef,et; semilog = true) 
    end
    Plots.savefig(  Plots.plot(p[1],p[2],p[3], layout = Plots.@layout([A B C]), 
                    xlabel = "A/λ", ylabel = "EB", ylims = (-6,0),
                    size = (1500,800)), "3d2d.png") 
end


