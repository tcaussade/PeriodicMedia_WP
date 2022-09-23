# /home/thomas/Downloads/julia-1.7.2/bin/julia
# import Pkg; Pkg.activate("/home/thomas/PeriodicMedia_WP")

using PeriodicMedia
import Plots

error("pending update...")

function meshconvergence(P::Problem,WGF::Window, meshsizes, dim)

    Shape = PeriodicMedia.ParametricSurfaces.Sphere
    Fig = Obstacle(Shape,minimum(P.L)/4)

    err = Vector{Float64}(undef,0)
    for ppw in meshsizes
        println("Solving with ppw = "*string(ppw))
        Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
        Γt = extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
        ϕ  = solver(P,Γs,Γt,WGF; FRO = false)
        @show eb = energytest(P,Γt,WGF, ϕ; FRO = false, H = WGF.c*WGF.A)
        push!(err,eb)
    end
    return err
end


setup = "3D2D"
# p = Vector{Plots.Plot}(undef,3)

if setup == "3D2D"
    # k1 = [9.0, 9.199221756451442, 9.5]
    k1 = [9.0]
    θ = [π/4.,π/4.]
    L = [0.5, 0.5]
    
    ####################################
    P = Problem([k1[1],15.], θ, L; ambdim = 3, geodim = 2)
    meshsizes = collect(1:3)
    dim       = 1
    wdwsizes  = collect(5:5:15) 
    ####################################

    p = Plots.plot()
    λ = (2π/k1[1])
    for w in wdwsizes
        WGF = Window(0.3, w*λ)
        e   = meshconvergence(P,WGF,meshsizes,dim)
        p   = Plots.plot!(meshsizes, log10.(e), ylims = (-4.5,-0.5), label = "A/λ="*string(WGF.A/λ))
    end
    Plots.savefig(p, "meshconv.png") 
end


