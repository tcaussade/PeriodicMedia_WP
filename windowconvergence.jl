# /home/thomas/Downloads/julia-1.7.2/bin/julia
# import Pkg; Pkg.activate("/home/thomas/PeriodicMedia_WP")

using PeriodicMedia
import Plots

function solver(E,MB,Wa,b; setup)
    if setup == "2D1D"
        ϕ = (E + MB*Wa)\b
    elseif setup == "3D2D"
        # GMRES parameters
        res = size(MB,2)
        vbs = false
        tol = 1e-6
        # Solve linear system
        ϕ,niter = gmres(E+MB*Wa,b; restart = res, verbose = vbs, reltol = tol, log = true)
        @info niter
    end
    return ϕ
end
function accuracy(P,Γt,Wpar,ϕ; correct)
    # @info "Computing energy" he
    @assert he < Wpar.c * Wpar.A
    R,T = energytest(P,Γt,Wpar, ϕ; FRO = correct, h = he)
    eb = abs(R+T-1)
    # @info "Test results" eb R T 
    return eb
end
function run_experiment(P,Wpar; etest)
    
    Γs = unitcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
    Γt = extendedcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
    G  = (Γs=Γs, Γt=Γt)

    if setup == "2D1D"
        smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2
    elseif setup == "3D2D"
        smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2+length(Γs[4].dofs)*2
    end
    @show smat

    # @info "Geometry created" Wpar Wpar.c*Wpar.A smat
    # @info "Assembling..."

    # Solving for the scattered field
    MB,Wa,E = matrixcreator(P,G,Wpar)
    b       = rightside(P,Γs)

    if etest == true
        ϕ = solver(E,MB,Wa,b; setup)
        ebf = accuracy(P,Γt,Wpar,ϕ; correct = false)
    end

    # Compute again with corrections
    hcorr = Wpar.c * Wpar.A
    # @info "Adding corrections..." δ hcorr
    MB += finiterankoperator(P,Γs,Γt; δ = δ , h = hcorr)

    ϕ = solver(E,MB,Wa,b; setup)
    if etest == true
        ebt = accuracy(P,Γt,Wpar,ϕ; correct = true)
    else
        return ϕ
    end

    return ebf,ebt
end

""" 
    set desired experiment
    - setup = "2D1D"
    - setup = "3D2D"
"""

global setup = "3D2D"

# set physical params and geometry
PeriodicMedia.clear_entities!()
if setup == "2D1D"
    θ   = π/4.
    L   = 2.0
    k1  = 10.72 
    k2  = 20.0
    pol = "TE"
    P = Problem([k1,k2],θ,L,pol; ambdim = 2, geodim = 1)

    Shape = PeriodicMedia.ParametricSurfaces.Disk
    Fig = Obstacle(Shape,L/4)
elseif setup == "3D2D"
    θ   = [0.,0.] 
    L   = [0.5, 0.5]
    k1  = 9.2 
    k2  = 15.0
    pol = "TE"
    P = Problem([k1,k2],θ,L,pol; ambdim = 3, geodim = 2)

    Shape = PeriodicMedia.ParametricSurfaces.Sphere
    Fig = Obstacle(Shape,minimum(L)/4)
end

global ppw = 1
global dim = 1

# correction parameters
global δ    = 2*k1
global he   = 1.0

# Window sizes (normalized to λ)
Asizes = collect(10:10:30)
errors = []

iter = 1
for Ap in Asizes
    # set discretization parameters
    @info "Iteration: "*string(iter)*"/"*string(length(Asizes)) Ap
    global λ    = 2π/k1
    c    = 0.3
    A    = Ap * λ
    Wpar = Window(c,A)

    @show ebf,ebt = run_experiment(P,Wpar; etest = true)
    push!(errors,[ebf,ebt])

    iter += 1
end

lbs   = ["without correction" "with correction"]
ef,et = [e[1] for e in errors],[e[2] for e in errors]

p = Plots.plot(title = "k="*string(P.pde[1].k))
p = Plots.plot!(Asizes, log10.([ef et]); label = lbs)

if savefig == true
    namefig = setup*string(".png")
    Plots.savefig(p, namefig) 
end



