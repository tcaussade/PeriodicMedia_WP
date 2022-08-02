"""
    Functions related to potential evaluation are defined here
"""

Nystrom.dofs(vec) = vec
Nystrom.coords(vec) = vec

""" transpotential computes transmitted field """

function transpotential(P::Problem{N,NP}, eval, G1::NystromMesh{N,T,M,NM}) where {N,T,M,NM,NP}
    [lmul!(-1.,Matrix(DoubleLayerOperator(P.pde[2],eval,G1))) Matrix(SingleLayerOperator(P.pde[2],eval,G1))]
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
function scatpotential(P::Problem{3,2}, eval, G::Vector)
    @assert G[1] isa Dict{Tuple{Int,Int}}
    D1 = sum([DoubleLayerOperator(P.pde[1], eval, get(G[1],(i,j),""))*P.γ[1]^i*P.γ[2]^j for i=-1:1,j=-1:1])
    S1 = lmul!(-P.η, sum([SingleLayerOperator(P.pde[1], eval, get(G[1],(i,j),""))*P.γ[1]^i*P.γ[2]^j for i=-1:1,j=-1:1]))

    D2 = lmul!(+1/P.γ[1], sum([DoubleLayerOperator(P.pde[1],eval,get(G[2],i,""))*P.γ[2]^i for i=-1:1]))
    S2 = lmul!(-1/P.γ[1], sum([SingleLayerOperator(P.pde[1],eval,get(G[2],i,""))*P.γ[2]^i for i=-1:1]))
    D3 = lmul!(-P.γ[1]^2, sum([DoubleLayerOperator(P.pde[1],eval,get(G[3],i,""))*P.γ[2]^i for i=-1:1]))
    S3 = lmul!(+P.γ[1]^2, sum([SingleLayerOperator(P.pde[1],eval,get(G[3],i,""))*P.γ[2]^i for i=-1:1]))
    
    D4 = lmul!(+1/P.γ[2], sum([DoubleLayerOperator(P.pde[1],eval,get(G[4],i,""))*P.γ[1]^i for i=-1:1]))
    S4 = lmul!(-1/P.γ[2], sum([SingleLayerOperator(P.pde[1],eval,get(G[4],i,""))*P.γ[1]^i for i=-1:1]))
    D5 = lmul!(-P.γ[2]^2, sum([DoubleLayerOperator(P.pde[1],eval,get(G[5],i,""))*P.γ[1]^i for i=-1:1]))
    S5 = lmul!(+P.γ[2]^2, sum([SingleLayerOperator(P.pde[1],eval,get(G[5],i,""))*P.γ[1]^i for i=-1:1]))

    [D1 S1 axpy!(1.,D2,D3) axpy!(1.,S2,S3) axpy!(1.,D4,D5) axpy!(1.,S4,S5)]
end

""" gradscatpotential computes the normal derivative over a given are """

function gradscatpotential(P::Problem{N,1}, eval, G::Vector) where N
    @assert G[1] isa Dict
    T1 = sum([HyperSingularOperator(P.pde[1], eval, get(G[1],i,"")) * P.γ^i for i=-1:1])
    K1 = -P.η * sum([AdjointDoubleLayerOperator(P.pde[1], eval, get(G[1],i,"")) * P.γ^i for i=-1:1])
    T2 = HyperSingularOperator(P.pde[1],eval,G[2]) / P.γ - HyperSingularOperator(P.pde[1],eval,G[3]) * P.γ^2
    K2 = AdjointDoubleLayerOperator(P.pde[1],eval,G[3]) * P.γ^2 - AdjointDoubleLayerOperator(P.pde[1],eval,G[2]) / P.γ 
    [T1 K1 T2 K2]
end

""" 
    cellsolution(P,G,w; ppw)
        G must be extended cell
        ppw is points per wavelength
        returns X,Y,U for plotting
    cutcellsolution
        is used for 3D to retrieve a given plane
"""

function cellsolution(P::Problem{2,1},G::Vector,w::Window,densities::Vector{ComplexF64}; ppw, ylims = [-2P.L 2P.L], FRO::Bool)
    @assert G[1] isa Dict

    h = 2π/P.pde[1].k / ppw
    X = -0.5*P.L:h:0.5*P.L
    Y = ylims[1]:h:ylims[2]
    mshgrid = vec([[x,y]  for x in X, y in Y])

    σw = lmul!(wgfmatrix(G,w),densities)
    us = scatpotential(P,mshgrid,G)*σw
    if FRO
        H = (P.L + w.A*w.c)*0.5
        us += scatcorrection(P,G,mshgrid,σw; H=H, δ = 0.75*P.pde[1].k)
    else
        @info "Non-corrected potential"
    end
    Γ₁ = get(G[1],0,"")
    ut = transpotential(P,mshgrid,Γ₁)*densities[1:2*length(Γ₁.dofs)]

    U = Matrix{ComplexF64}(undef,length(mshgrid),1)
    for n = 1:length(mshgrid)
        xp,yp = mshgrid[n]
        U[n] = Nystrom.isinside((xp,yp),Γ₁) ? ut[n] :  us[n] + uinc(P,[[xp,yp]])[1] # 
    end

    return X,Y,U

end
function cellsolution(P::Problem{3,2},G::Vector,w::Window,densities::Vector{ComplexF64}; ppw, zlims = [-w.c*w.A w.c*w.A])
    h = 2π/max(P.pde[1].k,P.pde[2].k) / ppw
    X = -0.5*P.L[1]:h:0.5*P.L[1]
    Y = -0.5*P.L[2]:h:0.5*P.L[2]
    Z = zlims[1]:h:zlims[2]

    W = wgfmatrix(G,w)
    
    Uxz = cutcellsolution(P,vec([SVector(x,0.,z)  for x in X, z in Z]),W,densities,G)
    Uyz = cutcellsolution(P,vec([SVector(0.,y,z)  for y in Y, z in Z]),W,densities,G)
    Uxy = cutcellsolution(P,vec([SVector(x,y,0.)  for x in X, y in Y]),W,densities,G)

    return X,Y,Z, (XZ = Uxz, YZ = Uyz, XY = Uxy)
end
function cutcellsolution(P::Problem{3,2},mshgrid,wgfmat::Diagonal,densities::Vector{ComplexF64},G::Vector)
    Γ₁ = get(G[1],(0,0),"")
    us = scatpotential(P,mshgrid,G)*wgfmat*densities
    ut = transpotential(P,mshgrid,Γ₁)*densities[1:2*length(Γ₁.dofs)]
    U = Matrix{ComplexF64}(undef,length(mshgrid),1)
    for n = 1:length(mshgrid)
        xp,yp,zp = mshgrid[n]
        U[n] = Nystrom.isinside((xp,yp,zp),Γ₁) ? ut[n] :  us[n] + uinc(P,[[xp,yp,zp]])[1] # 
    end
    return U
end

"""
    viewsolution
        heatmaps an arbitrary number of cells
"""
function viewsolution(P::Problem{2,1},X,Y,U::Matrix{ComplexF64},Fig::Obstacle; ncell, part = real)
    p = Plots.heatmap(X,Y,part.(U), clims = (-2,2), aspect_ratio = 1)
    for n in ncell
        Xn = X .+ n*P.L
        Un = U*P.γ^n
        p  = Plots.heatmap!(Xn,Y,part.(Un), color = :RdBu, clims = (-2,2))
        F  = Scatterer(50, Fig; c = (n*P.L,0.), dimorder = 2)
        p  = Plots.plot!(F, color = :black)
    end
    return p
end
function XZviewsolution(P::Problem{3,2},X,Z,U::Matrix{ComplexF64}; ncell, part = real)
    p = Plots.heatmap(X,Z,part.(U), clims = (-2,2), aspect_ratio = 1, title = "XZ view")
    for n in ncell
        Xn = X .+ n*P.L[1]
        Un = U*P.γ[1]^n
        p  = Plots.heatmap!(Xn,Z,part.(Un), color = :RdBu, clims = (-2,2))
    end
    return p
end
function YZviewsolution(P::Problem{3,2},Y,Z,U::Matrix{ComplexF64}; ncell, part = real)
    p = Plots.heatmap(Y,Z,part.(U), clims = (-2,2), aspect_ratio = 1, title = "YZ view")
    for n in ncell
        Yn = Y .+ n*P.L[2]
        Un = U*P.γ[2]^n
        p  = Plots.heatmap!(Yn,Z,part.(Un), color = :RdBu, clims = (-2,2))
    end
    return p
end
function XYviewsolution(P::Problem{3,2},X,Y,U::Matrix{ComplexF64}; ncell, part = real)
    p = Plots.heatmap(X,Y,part.(U), clims = (-2,2), aspect_ratio = 1, title = "XY view")
    for n1 in ncell, n2 in ncell
        Xn = X .+ n1*P.L[1]
        Yn = Y .+ n2*P.L[2]
        Un = U*P.γ[1]^n1*P.γ[2]^n2
        p  = Plots.heatmap!(Xn,Yn,part.(Un), color = :RdBu, clims = (-2,2))
    end
    return p
end