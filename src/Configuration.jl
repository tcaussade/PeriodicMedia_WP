"""
    All functions/structs related to geometry construction are defined here
"""

struct Obstacle
    shape  :: DataType
    radius :: Float64
    # center :: 
end

function StraightLine(xs,xe; M::Int, dimorder::Int, qrule = WavePropBase.Fejer)
    Γ   = Domain(Nystrom.line(xs,xe))
    msh = Nystrom.meshgen(Γ,M)
    NystromMesh(msh,Γ, qrule(dimorder))
end

function Scatterer(M::Int, Fig::Obstacle; c, dimorder::Int)
    Ω   = Fig.shape(;radius = Fig.radius, center = c)
    ∂Ω  = boundary(Domain(Ω))
    msh = Nystrom.meshgen(∂Ω, M)
    NystromMesh(msh,∂Ω; order = dimorder)
end

function unitcell(P::Problem{2,1}, Fig::Obstacle, w::Window ; ppw::Int,dimorder::Int)
    Mv = Int(ceil( 2*w.A*(P.pde[1].k/2π) * ppw ))
    Ms = Int(ceil( 2π*Fig.radius*(P.pde[1].k/2π) * ppw ))
    Γ₁ = Scatterer(Ms, Fig; c = (0.,0.), dimorder = dimorder)
    Γ₂ = StraightLine((-0.5*P.L,-w.A),(-0.5*P.L,w.A); M = Mv, dimorder = dimorder)
    Γ₃ = StraightLine((+0.5*P.L,-w.A),(+0.5*P.L,w.A); M = Mv, dimorder = dimorder)
    return [Γ₁,Γ₂,Γ₃]
end

function extendedcell(P::Problem{2,1}, Fig::Obstacle, w::Window ;ppw::Int,dimorder::Int)
    Mv = Int(ceil( 2*w.A*(P.pde[1].k/2π) * ppw ))
    Ms = Int(ceil( 2π*Fig.radius*(P.pde[1].k/2π) * ppw ))
    Γ₁ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = Scatterer(Ms, Fig; c = (i*P.L,0.), dimorder = dimorder)
        merge!(Γ₁,Dict(i => Γ))
    end
    Γ₂ = StraightLine((-1.5*P.L,-w.A),(-1.5*P.L,w.A); M = Mv, dimorder = dimorder)
    Γ₃ = StraightLine((+1.5*P.L,-w.A),(+1.5*P.L,w.A); M = Mv, dimorder = dimorder)
    return [Γ₁, Γ₂, Γ₃]
end

"""
    3D2D functions
"""

function meshplot!(Γ::NystromMesh{3,T,M,NM}) where {T,M,NM}
    x = [q.coords[1] for q in Γ.dofs]
    y = [q.coords[2] for q in Γ.dofs]
    z = [q.coords[3] for q in Γ.dofs]
    Plots.scatter!(x,y,z, legend = false)
end
function meshplot(Γ::NystromMesh{3,T,M,NM}) where {T,M,NM}
    x = [q.coords[1] for q in Γ.dofs]
    y = [q.coords[2] for q in Γ.dofs]
    z = [q.coords[3] for q in Γ.dofs]
    Plots.scatter(x,y,z, legend = false)
end

function StraightPlane(xe::SVector,xs::SVector; M::Tuple{Int,Int}, dimorder::Int64, qrule=WavePropBase.Fejer)
    f = (u) -> xe + [u[1], u[1], u[2]].*(xs - xe)
    d = Nystrom.HyperRectangle((0.,0.),(1.,1.))
    Γ = Domain(Nystrom.ParametricEntity(f,d))
    msh = Nystrom.meshgen(Γ,M)
    NystromMesh(msh,Γ, WavePropBase.TensorProductQuadrature(qrule(dimorder),qrule(dimorder)))
end
StraightPlane(xe,xs;M,dimorder,qrule=WavePropBase.Fejer) = StraightPlane(SVector(xe),SVector(xs);M,dimorder,qrule=WavePropBase.Fejer)

function HorizontalStraightPlane(xe::SVector,xs::SVector; M::Tuple{Int,Int}, dimorder::Int64, qrule=WavePropBase.Fejer)
    f = (u) -> xe + [u[1], u[2], xe[3]].*(xe - xs)
    d = Nystrom.HyperRectangle((0.,0.),(1.,1.))
    Γ = Domain(Nystrom.ParametricEntity(f,d))
    msh = Nystrom.meshgen(Γ,M)
    NystromMesh(msh,Γ, WavePropBase.TensorProductQuadrature(qrule(dimorder),qrule(dimorder)))
end
HorizontalStraightPlane(xe,xs;M,dimorder,qrule) = HorizontalStraightPlane(SVector(xe),SVector(xs);M,dimorder,qrule)

function unitcell(P::Problem{3,2},Fig::Obstacle, w::Window ; ppw::Int,dimorder::Int)
    N1 = Int(ceil( sqrt(4π*Fig.radius^2 /6) * (P.pde[1].k/2π) * ppw ))
    Nx = (Int(ceil(P.L[1]*(P.pde[1].k/2π) * ppw)), Int(ceil(2*w.A*(P.pde[1].k/2π) * ppw)))
    Ny = (Int(ceil(P.L[2]*(P.pde[1].k/2π) * ppw)), Int(ceil(2*w.A*(P.pde[1].k/2π) * ppw)))
    Γ₁ = Scatterer(N1,Fig; c = (0.,0.,0.), dimorder = dimorder)
    Γ₂ = StraightPlane((-0.5P.L[1],-0.5P.L[2],-w.A),(-0.5P.L[1],+0.5P.L[2],+w.A); M=Nx, dimorder = dimorder)
    Γ₃ = StraightPlane((+0.5P.L[1],-0.5P.L[2],-w.A),(+0.5P.L[1],+0.5P.L[2],+w.A); M=Nx, dimorder = dimorder)
    Γ₄ = StraightPlane((+0.5P.L[1],-0.5P.L[2],-w.A),(-0.5P.L[1],-0.5P.L[2],+w.A); M=Ny, dimorder = dimorder)
    Γ₅ = StraightPlane((+0.5P.L[1],+0.5P.L[2],-w.A),(-0.5P.L[1],+0.5P.L[2],+w.A); M=Ny, dimorder = dimorder)
    return [Γ₁,Γ₂,Γ₃,Γ₄,Γ₅]
end

function extendedcell(P::Problem{3,2},Fig::Obstacle, w::Window ; ppw::Int,dimorder::Int)
    N1 = Int(ceil( sqrt(4π*Fig.radius^2 /6) * (P.pde[1].k/2π) * ppw ))
    Nx = (Int(ceil(P.L[1]*(P.pde[1].k/2π) * ppw)), Int(ceil(2*w.A*(P.pde[1].k/2π) * ppw)))
    Ny = (Int(ceil(P.L[2]*(P.pde[1].k/2π) * ppw)), Int(ceil(2*w.A*(P.pde[1].k/2π) * ppw)))
    Γ₁ = Dict{Tuple{Int,Int},NystromMesh}()
    for i = -1:1, j=-1:1
        Γ = Scatterer(N1, Fig; c = (i*P.L[1],j*P.L[2],0.), dimorder = dimorder)
        merge!(Γ₁, Dict((i,j) => Γ))
    end
    Γ₂ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = StraightPlane((-1.5P.L[1],(-0.5+i)P.L[2],-w.A),(-1.5P.L[1],(+0.5+i)P.L[2],+w.A); M=Nx, dimorder = dimorder)
        merge!(Γ₂, Dict(i => Γ))
    end
    Γ₃ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = StraightPlane((+1.5P.L[1],(-0.5+i)P.L[2],-w.A),(+1.5P.L[1],(+0.5+i)P.L[2],+w.A); M=Nx, dimorder = dimorder)
        merge!(Γ₃, Dict(i => Γ))
    end
    Γ₄ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = StraightPlane(((+0.5+i)P.L[1],-1.5P.L[2],-w.A),((-0.5+i)P.L[1],-1.5P.L[2],+w.A); M=Ny, dimorder = dimorder)
        merge!(Γ₄, Dict(i => Γ))
    end
    Γ₅ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = StraightPlane(((+0.5+i)P.L[1],+1.5P.L[2],-w.A),((-0.5+i)P.L[1],+1.5P.L[2],+w.A); M=Ny, dimorder = dimorder)
        merge!(Γ₅, Dict(i => Γ))
    end
    return [Γ₁,Γ₂,Γ₃,Γ₄,Γ₅]
end