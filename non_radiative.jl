# /home/thomas/Downloads/julia-1.7.2/bin/julia
# import Pkg; Pkg.activate("/home/thomas/PeriodicMedia_WP")

using PeriodicMedia

function solver(E,MB,Wa,b)
    # return (E+MB*Wa)\b
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
    # dR⁺, R⁺, x⁺ = PeriodicMedia.linegradscat(P,G; H = +h)
    # dR⁻, R⁻, x⁻ = PeriodicMedia.linegradscat(P,G; H = -h)
    C⁺ = []
    C⁻ = []
    B⁺ = []
    B⁻ = []
    betas = []
    for n ∈ Cδ
        _,_,βₙ = PeriodicMedia.seriesconstant(P,(n[1],n[2]))

        Lₙ⁺ = PeriodicMedia.Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.0) * σw
        Lₙ⁻ = PeriodicMedia.Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.0) * σw

        Cₙ⁺ = -exp(im*βₙ*h)/(2*im*βₙ) * Lₙ⁺ 
        Cₙ⁻ = +exp(im*βₙ*h)/(2*im*βₙ) * Lₙ⁻ 
        push!(C⁺,Cₙ⁺)
        push!(C⁻,Cₙ⁻)

        Bₙ⁺ = +exp(-im*βₙ*h)/(2*im*βₙ) * Lₙ⁻ 
        Bₙ⁻ = -exp(-im*βₙ*h)/(2*im*βₙ) * Lₙ⁺ 
        push!(B⁺,Bₙ⁺)
        push!(B⁻,Bₙ⁻)

        push!(betas,βₙ)
    end
    return betas, C⁺,C⁻, B⁺,B⁻
end
function etest(bp,bm)
    α₁,α₂,β = P.pde[1].k*P.dir
    S = zero(ComplexF64)
    R = zero(ComplexF64)
    for (i,(n1,n2)) in enumerate(Cδ)
        αₙ₁ = α₁ + n1*2π/P.L[1]
        αₙ₂ = α₂ + n2*2π/P.L[2]
        βₙ = sqrt(complex(P.pde[1].k^2 - αₙ₁^2 - αₙ₂^2 ))
        if P.pde[1].k^2 ≥ abs(αₙ₁)^2 + abs(αₙ₂)^2
            βₙ = sqrt(complex(P.pde[1].k^2 - αₙ₁^2 - αₙ₂^2 ))
            # @show B⁺ = exp(-im*βₙ*h)/P.L[1]/P.L[2] * sum( u1.*exp.(-im*αₙ₁*xi-im*αₙ₂*yi) ) *P.L[1]*P.L[2]/length(u1)
            # @show B⁻ = exp(-im*βₙ*h)/P.L[1]/P.L[2] * sum( u2.*exp.(-im*αₙ₁*xi-im*αₙ₂*yi) ) *P.L[1]*P.L[2]/length(u1)
            B⁺ = bp[i]
            B⁻ = bm[i]
            if n1 == 0 && n2 == 0
                global B₀⁻ = B⁻
            end
            S += βₙ/β * (abs(B⁺)^2+abs(B⁻)^2)
            R += βₙ/β * abs(B⁺)^2
            # T += βₙ/β * abs(B⁻)^2
        end
    end
    e = abs(2*real(B₀⁻)+S)
    # T = 1 + 2*real(B₀⁻) + S - R
    # abs(2*real(B₀)+S)
    return e
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

    nonradiative(P,Γt,Wa*ϕ; h = he)
end

PeriodicMedia.clear_entities!()

# θ   = [π/4.,π/4.] 
θ = [0.,0.]
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


global ppw = 1
global dim = 1

# correction parameters
       δ    = 1.5*k1
global he   = 1.0

global Cδ = PeriodicMedia.fixdelta(P,δ)

# Window sizes (normalized to λ)
Asizes = collect(3:1:5)
cp  = Matrix{Float64}(undef,length(Cδ),length(Asizes))
cm  = Matrix{Float64}(undef,length(Cδ),length(Asizes))
e   = Vector{Float64}(undef,length(Asizes))
betas = Vector{ComplexF64}(undef,length(Cδ))

name_exp = "non_radiative"*"_ppw"*string(ppw)*"_dim"*string(dim)*"_k"*string(k1)*string(".txt")

for (i,Ap) in enumerate(Asizes)
    # set discretization parameters
    @info Ap
    global λ = 2π/k1
    c    = 0.3
    A    = Ap * λ
    Wpar = Window(c,A)

    betas, c⁺,c⁻,b⁺,b⁻ = run_experiment(P,Wpar)

    @show e[i]    = etest(b⁺,b⁻)
    cp[:,i] .= abs.(c⁺)
    cm[:,i] .= abs.(c⁻)
end

for (i,n) in enumerate(Cδ)
    @show (n,betas[i])
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

eb = Plots.plot(Asizes, log10.(e); title = "EB")

p = Plots.plot(pp,pm,eb, layout = (3,1), size = (600,800))

namefig = "radiative_coefs"*string(".png")
Plots.savefig(p, namefig) 




