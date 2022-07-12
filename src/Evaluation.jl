"""
    Functions related to potential evaluation are defined here
"""

Nystrom.dofs(vec) = vec
Nystrom.coords(vec) = vec

function transpotential(P::Problem{N,1}, eval, G1::NystromMesh{N,T,M,NM}) where {N,T,M,NM}
    [-DoubleLayerOperator(P.pde[2],eval,G1) SingleLayerOperator(P.pde[2],eval,G1)]
end

"""
    scatpotential(P,eval,G)
        P: parameters
        eval: where you want to evaluate
        G: Must be triplecell in order to avoid nearly singular in artificial boundaries
    returns potential operator as a matrix
"""

function scatpotential(P::Problem{N,1}, eval, G::Vector) where N
    @assert G[1] isa Dict
    D1 = sum([DoubleLayerOperator(P.pde[1], eval, get(G[1],i,"")) * P.γ^i for i=-1:1])
    S1 = -P.η * sum([SingleLayerOperator(P.pde[1], eval, get(G[1],i,"")) * P.γ^i for i=-1:1])
    D2 = DoubleLayerOperator(P.pde[1],eval,G[2]) / P.γ - DoubleLayerOperator(P.pde[1],eval,G[3]) * P.γ^2
    S2 = SingleLayerOperator(P.pde[1],eval,G[3]) * P.γ^2 - SingleLayerOperator(P.pde[1],eval,G[2]) / P.γ 
    [D1 S1 D2 S2]
end

function gradscatpotential(P::Problem{N,1}, eval, G::Vector) where N
    @assert G[1] isa Dict
    T1 = sum([HyperSingularOperator(P.pde[1], eval, get(G[1],i,"")) * P.γ^i for i=-1:1])
    K1 = -P.η * sum([AdjointDoubleLayerOperator(P.pde[1], eval, get(G[1],i,"")) * P.γ^i for i=-1:1])
    T2 = HyperSingularOperator(P.pde[1],eval,G[2]) / P.γ - HyperSingularOperator(P.pde[1],eval,G[3]) * P.γ^2
    K2 = AdjointDoubleLayerOperator(P.pde[1],eval,G[3]) * P.γ^2 - AdjointDoubleLayerOperator(P.pde[1],eval,G[2]) / P.γ 
    [T1 K1 T2 K2]
end

""" 
    Plots
    cellsolution(P,G,w; ppw)
        G must be triple cell
        ppw is points per wavelength
    returns X,Y,U for plotting
"""

function cellsolution(P::Problem{2,1},G::Vector,w::Window,densities::Vector{ComplexF64}; ppw, ylims = [-2P.L 2P.L])
    @assert G[1] isa Dict

    h = 2π/P.pde[1].k / ppw
    X = -0.5*P.L:h:0.5*P.L
    Y = ylims[1]:h:ylims[2]
    mshgrid = vec([[x,y]  for x in X, y in Y])

    W = wgfmatrix(G,w)
    us = scatpotential(P,mshgrid,G)*W*densities
    Γ₁ = get(G[1],0,"")
    ut = transpotential(P,mshgrid,Γ₁)*densities[1:2*length(Γ₁.dofs)]

    U = Matrix{ComplexF64}(undef,length(mshgrid),1)
    for n = 1:length(mshgrid)
        xp,yp = mshgrid[n]
        U[n] = Nystrom.isinside((xp,yp),Γ₁) ? ut[n] :  us[n] + uinc(P,[[xp,yp]])[1] # 
    end

    return X,Y,U

end

function viewsolution(P::Problem{2,1},X,Y,U::Matrix{ComplexF64},Fig::Obstacle; ncell, part = real)
    p = heatmap(X,Y,part.(U), clims = (-2,2), aspect_ratio = 1)
    for n in ncell
        Xn = X .+ n*P.L
        Un = U*P.γ^n
        p  = heatmap!(Xn,Y,part.(Un), color = :RdBu, clims = (-2,2))
        F  = Scatterer(50, Fig; c = (n*P.L,0.), dimorder = 2)
        p  = plot!(F, color = :black)
    end
    return p
end