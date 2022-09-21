"""
    Base functionalities are defined 
"""

"""
    Problem(k,θ,L) computes physical parameters

    where k: [k₁,k₂] for inner and outer obstacle media, 
          θ: incidence angle,
          L: unit cell size.
          ambdim: ambient dimension (d=2 for problems in R² )
          geodim: periodic dimensions (e.g. d=1 for 1D photonic crystals)
"""

struct Problem{N,NP}
    pde :: Vector{Nystrom.Helmholtz{N, Float64}}
    dir :: Vector{Float64}
    η   :: Float64
    γ   :: Union{ComplexF64,Vector{ComplexF64}}
    L   :: Union{Float64,Vector{Float64}}
    function Problem(k::Vector{Float64},θ::Float64,L::Float64, pol::String; ambdim::Int, geodim::Int)
        @assert ambdim==2 && geodim==1 
        pde = [Helmholtz(dim=ambdim; k=k[1]), Helmholtz(dim=ambdim; k=k[2])]
        dir = [sin(θ),cos(θ)]
        γ   = exp(im*k[1]*dir[1]*L)
        pol == "TE" ? η = 1.0 : η = k[2]^2/k[1]^2
        @info "Created physical problem" (dir[1],dir[2]) γ η
        new{ambdim,geodim}(pde,dir,η,γ,L)
    end
    function Problem(k::Vector{Float64},θ::Vector{Float64},L::Vector{Float64}, pol::String; ambdim::Int, geodim::Int)
        @assert ambdim==3 && geodim==2 
        pde = [Helmholtz(dim=ambdim; k=k[1]), Helmholtz(dim=ambdim; k=k[2])]
        dir = [sin(θ[1])*cos(θ[2]),sin(θ[1])*sin(θ[2]), cos(θ[1])]
        γ   = [exp(im*k[1]*dir[1]*L[1]),exp(im*k[1]*dir[2]*L[2])]
        pol == "TE" ? η = 1.0 : η = k[2]^2/k[1]^2
        @info "Created physical problem" (dir[1],dir[2],dir[3]) (γ[1],γ[2]) η
        new{ambdim,geodim}(pde,dir,η,γ,L)
    end
end

"""
    uinc and ∂uinc return incident field and normal incidient field respectively
"""

function uinc(P::Problem{2,1},x::Vector)
    r = [P.dir[1]*xi[1] - P.dir[2]*xi[2] for xi in x]
    exp.(im*P.pde[1].k*r)
end

function uinc(P::Problem{3,NP},x::Vector) where NP
    r = [P.dir[1]*xi[1] + P.dir[2]*xi[2] - P.dir[3]*xi[3] for xi in x]
    exp.(im*P.pde[1].k*r)
end

function ∂uinc(P::Problem{2,1},x::Vector,ν::Vector)
    r = [P.dir[1]*xi[1] - P.dir[2]*xi[2] for xi in x]
    n = [P.dir[1]*νi[1] - P.dir[2]*νi[2] for νi in ν]
    im*P.pde[1].k*n.*exp.(im*P.pde[1].k*r)
end

function ∂uinc(P::Problem{3,NP},x::Vector,ν::Vector) where NP
    r = [P.dir[1]*xi[1] + P.dir[2]*xi[2] - P.dir[3]*xi[3] for xi in x]
    n = [P.dir[1]*νi[1] + P.dir[2]*νi[2] - P.dir[3]*νi[3] for νi in ν]
    im*P.pde[1].k*n.*exp.(im*P.pde[1].k*r)
end

"""
    Windowed Green function associated functionalities
"""

struct Window
    c :: Float64
    A :: Float64
    # function Window(c,A)
    #     abs(c)<1.0 & A>0 ? new(c,A) : error("Window not admitted")
    # end
end

function χ(x,w::Window)
    if abs(x) <= w.c*w.A
        return one(x)
    elseif w.c*w.A < abs(x) < w.A
        u = (abs(x)-w.c*w.A)/(w.A-w.c*w.A)
        return exp( 2*exp(-1/u)/(u-1) )
    else
        return zero(x)
    end
end
        
function wgfmatrix(G::Vector{NystromMesh{2,T,M,NM}}, w::Window) where {T,M,NM}
    d1 = ones(ComplexF64, length(G[1].dofs))
    d2 = [χ(q.coords[2],w) for q in G[2].dofs]
    Diagonal([d1;d1;d2;d2])
end

function wgfmatrix(G::Vector{NystromMesh{3,T,M,NM}}, w::Window) where {T,M,NM}
    d1 = ones(ComplexF64, length(G[1].dofs))
    d2 = [χ(q.coords[3],w) for q in G[2].dofs]
    d3 = [χ(q.coords[3],w) for q in G[4].dofs]
    Diagonal([d1;d1;d2;d2;d3;d3])
end

function wgfmatrix(G::Vector, w::Window)
    if G[1] isa Dict{Int}
        d1 = ones(ComplexF64, length(get(G[1],0,"").dofs))
        d2 = [χ(q.coords[2],w) for q in G[2].dofs]
        return Diagonal([d1;d1;d2;d2])
    elseif G[1] isa Dict{Tuple{Int,Int}}
        d1 = ones(ComplexF64, length(get(G[1],(0,0),"").dofs))
        d2 = [χ(q.coords[3],w) for q in get(G[2],0,"").dofs]
        d3 = [χ(q.coords[3],w) for q in get(G[4],0,"").dofs]
        return Diagonal([d1;d1;d2;d2;d3;d3])
    end
end