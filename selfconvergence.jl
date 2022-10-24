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

function run_experiment(P,Wpar; ref)
    
    Γs = unitcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
    Γt = extendedcell(P,Fig, Wpar; ppw = ppw, dimorder = dim)
    G  = (Γs=Γs, Γt=Γt)

    smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2+length(Γs[4].dofs)*2
    @show smat

    # Solving for the scattered field
    MB,Wa,E = matrixcreator(P,G,Wpar)
    b       = rightside(P,Γs)
    # MB += finiterankoperator(P,Γs,Γt; δ = δ , h = hcorr)
    ϕ = solver(E,MB,Wa,b)

    # evaluate at reference point
    (PeriodicMedia.scatpotential(P,[ref],Γt) * ϕ)[1]
end

PeriodicMedia.clear_entities!()

θ   = [0.,0.] 
L   = [0.5, 0.5]
k1  = 8.8 
k2  = 14.0
pol = "TE"
P = Problem([k1,k2],θ,L,pol; ambdim = 3, geodim = 2)

Shape = PeriodicMedia.ParametricSurfaces.Sphere
Fig = Obstacle(Shape,minimum(L)/4)

# correction parameters
global δ    = 2*k1
global he   = 0.25

# Window sizes (normalized to λ)
Asizes = collect(3:.15:12)

global ppw = 8
global λ = 2π/k1
global c = 0.3
@assert he < c * Asizes[1]*λ

dim_orders = [2,3,4]
for dim_order in dim_orders
    global dim = dim_order
    name_exp = "selfconvergence"*"_ppw"*string(ppw)*"_dim"*string(dim_order)*"_k"*string(k1)*string(".txt")

    # the "most precise" solution
    ref  = [0.,0.,he]
    uend = run_experiment(P,Window(c,15.0 * λ); ref = ref) 

    for Ap in Asizes
        Wpar = Window(c,Ap*λ)
        @info dim Ap Wpar
        utp = run_experiment(P,Wpar; ref = ref)
        @show e = abs(uend-utp)

        open(name_exp,"a") do io
            tail = Ap==Asizes[end] ? "" : "\n"
            write(io, string(e) * tail)
        end
    end
end