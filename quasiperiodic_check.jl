# /home/thomas/Downloads/julia-1.7.2/bin/julia
# import Pkg; Pkg.activate("/home/thomas/PeriodicMedia_WP")
using PeriodicMedia

function solver(E,MB,Wa,b)
    # GMRES parameters
    res = size(MB,2)
    vbs = false
    tol = 1e-6
    # Solve linear system
    ϕ,niter = gmres(E+MB*Wa,b; restart = res, verbose = vbs, reltol = tol, log = true)
    @info niter
    return ϕ
end

function evalpoints(P, h)
    top = Matrix{Matrix}(undef,3,3)
    bot = Matrix{Matrix}(undef,3,3)
    shf = 5/12
    for n =-1:1, m=-1:1
        top[n+2,m+2] = [PeriodicMedia.SVector(shf*i + n*P.L[1],shf*j + m*P.L[2], +h) for i=[0],j=[0]]
        bot[n+2,m+2] = [PeriodicMedia.SVector(shf*i + n*P.L[1],shf*j + m*P.L[2], -h) for i=[0],j=[0]]
    end
    return top,bot
end

function qptest(P,G,σw; h)
    ev = evalpoints(P, h)
    e = [];
    for pts in ev
        ref = PeriodicMedia.scatpotential(P,vec(pts[2,2]),G)*σw
        err = zeros(3,3)
        for n1=-1:1,n2=-1:1
            scat = PeriodicMedia.scatpotential(P,vec(pts[n1+2,n2+2]),G)*σw
            γ    = P.γ[1]^n1 * P.γ[2]^n2
            err[n1+2,n2+2] = maximum(abs.(γ * ref - scat)) /maximum(abs.(ref))

        end
        push!(e,maximum(err))
    end
    return maximum(e)
end

function run_experiment(P,Wpar)
    
    Γs = unitcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
    Γt = extendedcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
    G  = (Γs=Γs, Γt=Γt)

    smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2+length(Γs[4].dofs)*2
    @show smat

    MB,Wa,E = matrixcreator(P,G,Wpar)
    b       = rightside(P,Γs)
    ϕ       = solver(E,MB,Wa,b)

    qptest(P,Γt,Wa*ϕ; h = he)
end


PeriodicMedia.clear_entities!()

θ   = [0., 0.] 
L   = [0.5, 0.5]
# θ = 0.0
# L = 0.5
k1  = 8.8 
k2  = 15.0
pol = "TE"
P = Problem([k1,k2],θ,L,pol; ambdim = 3, geodim = 2)
# P = Problem([k1,k2],θ,L,pol; ambdim = 2, geodim = 1)

Shape = PeriodicMedia.ParametricSurfaces.Sphere
# Shape = PeriodicMedia.ParametricSurfaces.Kite
Fig = Obstacle(Shape,minimum(L)/4)

# correction parameters
#        δ    = 1.5*k1
global he   = .25

# global Cδ = PeriodicMedia.fixdelta(P,δ)

# Window sizes (normalized to λ)
Asizes = collect(3:.15:12)
e = Vector{Float64}(undef,length(Asizes))

global ppw = 8
dim_orders = [2,3,4]
name_exp = "qp_k"*string(k1)*"_ppw"*string(ppw)*".txt"
for dim_ord in dim_orders
    global dim = dim_ord
    for (i,Ap) in enumerate(Asizes)
        # set discretization parameters
        @info dim Ap
        global λ = 2π/k1
        c    = 0.3
        A    = Ap * λ
        Wpar = Window(c,A)

        e[i] = run_experiment(P,Wpar)
        open(name_exp*"_p"*string(dim),"a") do io
            tail = Ap==Asizes[end] ? "" : "\n"
            write(io, string(e[i]) * tail)
        end
    end
end

################# Plot convergence #################
# import Plots
# Plots.plot(Asizes,log10.(e))
# namefig = "qp_convergence"*string(".png")
# Plots.savefig(p, namefig) 