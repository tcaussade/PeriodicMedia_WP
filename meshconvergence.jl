# /home/thomas/Downloads/julia-1.7.2/bin/julia
# import Pkg; Pkg.activate("/home/thomas/PeriodicMedia_WP")

using PeriodicMedia
import Plots


function meshconvergence(P::Problem,wgf::Window, meshsizes, dim)

    Shape = PeriodicMedia.ParametricSurfaces.Sphere
    Fig = Obstacle(Shape,minimum(P.L)/4)

    err = Vector{Float64}(undef,0)
    for ppw in meshsizes
        println("Solving with ppw = "*string(ppw))
        Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
        Γt = extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
        ϕ  = solver(P,Γs,Γt,WGF; FRO = false)
        @show eb = energytest(P,Γt,WGF, ϕ; FRO = false, H = wgf.c*wgf.A)
        push!(err,eb)
    end
    return err
end


setup = "3D2D"
p = Vector{Plots.Plot}(undef,3)

if setup == "3D2D"
    # k1 = [9.0, 9.199221756451442, 9.5]
    k1 = [9.0]
    θ = [π/4.,π/4.]
    L = [0.5, 0.5]
    
    ####################################
    P = Problem([k1[1],15.], θ, L; ambdim = 3, geodim = 2)
    meshsizes = collect(3:5)
    dim       = 2
    WGF       = Window(0.2, 15*(2π/k1[1]))
    ####################################

    e = meshconvergence(P,WGF,meshsizes,dim)
    p = Plots.plot(meshsizes, log10.(e), ylims = (-4.5,-0.5))

    Plots.savefig(p, "meshconv.png") 
end


