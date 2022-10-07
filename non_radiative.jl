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
function nonradiative(P,G,σw; h)
    dR⁺, R⁺, x⁺ = PeriodicMedia.planegradscat(P,G; H = +h)
    dR⁻, R⁻, x⁻ = PeriodicMedia.planegradscat(P,G; H = -h)
    C⁺ = []
    C⁻ = []
    for n ∈ Cδ
        _,_,βₙ = PeriodicMedia.seriesconstant(P,(n[1],n[2]))
        Cₙ⁺ = +exp(im*βₙ*h)/(2*im*βₙ) * PeriodicMedia.Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.0) * σw
        Cₙ⁻ = -exp(im*βₙ*h)/(2*im*βₙ) * PeriodicMedia.Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.0) * σw
        push!(C⁺,Cₙ⁺)
        push!(C⁻,Cₙ⁻)
    end
    return C⁺,C⁻
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

    C⁺,C⁻ = nonradiative(P,Γt,Wa*ϕ; h = he)
end

PeriodicMedia.clear_entities!()

θ   = [0.,0.] 
L   = [0.5, 0.5]
k1  = 8.8 
k2  = 15.0
pol = "TE"
P = Problem([k1,k2],θ,L,pol; ambdim = 3, geodim = 2)

Shape = PeriodicMedia.ParametricSurfaces.Sphere
Fig = Obstacle(Shape,minimum(L)/4)


global ppw = 2
global dim = 1

# correction parameters
       δ    = 1.5*k1
global he   = 1.0

global Cδ = PeriodicMedia.fixdelta(P,δ)

# Window sizes (normalized to λ)
Asizes = collect(3:1:5)
cp  = Matrix{Float64}(undef,length(Cδ),length(Asizes))
cm  = Matrix{Float64}(undef,length(Cδ),length(Asizes))

for (i,Ap) in enumerate(Asizes)
    # set discretization parameters
    @info Ap
    global λ = 2π/k1
    c    = 0.3
    A    = Ap * λ
    Wpar = Window(c,A)

    c⁺,c⁻ = run_experiment(P,Wpar)
    cp[:,i] .= abs.(c⁺)
    cm[:,i] .= abs.(c⁻)
end

################# Plot convergence #################
import Plots

pp = Plots.plot(title = "Cₙ⁺")
for (i,n) in enumerate(Cδ)
    pp = Plots.plot!(Asizes, log10.(cp[i,:]); label = "+"*string(n) )
end

pm = Plots.plot(title = "Cₙ⁻")
for (i,n) in enumerate(Cδ)
    pm = Plots.plot!(Asizes, log10.(cm[i,:]); label = "-"*string(n) )
end

p = Plots.plot(pp,pm, layout = (1,2))

namefig = "radiative_coefs"*string(".png")
Plots.savefig(p, namefig) 




