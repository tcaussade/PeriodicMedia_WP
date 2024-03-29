"""
    All functions/structs related to geometry construction are defined here
"""

""" functionalities for properly defining the obstacle 
    Obstacle is the struct containing the info
    Scatterer returns a NystromMesh given Obstacle data"""
struct Obstacle
    shape  :: DataType
    radius :: Float64
    # center :: 
end
function Scatterer(M::Int, Fig::Obstacle; c, dimorder::Int)
    Ω   = Fig.shape(;radius = Fig.radius, center = c)
    ∂Ω  = boundary(Domain(Ω))
    msh = Nystrom.meshgen(∂Ω, M)
    # NystromMesh(msh,∂Ω, WavePropBase.TrapezoidalOpen(dimorder))
    NystromMesh(msh,∂Ω; order= dimorder)
end

""" vertical/horizontal curves boundary (2D) """

function StraightLine(xs,xe; M::Int, dimorder::Int, qrule = WavePropBase.Fejer)
    Γ   = Domain(Nystrom.line(xs,xe))
    msh = Nystrom.meshgen(Γ,M)
    NystromMesh(msh,Γ, qrule(dimorder))
end

""" vertical surfaces boundary (3D)"""
function StraightPlane(xe::SVector,xs::SVector; M::Tuple{Int,Int}, dimorder::Int64, qrule=WavePropBase.Fejer)
    f = (u) -> xe + [u[1], u[1], u[2]].*(xs - xe)
    d = Nystrom.HyperRectangle((0.,0.),(1.,1.))
    Γ = Domain(Nystrom.ParametricEntity(f,d))
    msh = Nystrom.meshgen(Γ,M)
    NystromMesh(msh,Γ, WavePropBase.TensorProductQuadrature(qrule(dimorder),qrule(dimorder)))
end
StraightPlane(xe,xs;M,dimorder) = StraightPlane(SVector(xe),SVector(xs);M,dimorder)
""" horizontal surfaces boundary (3D)"""
function HorizontalStraightPlane(xe::SVector,xs::SVector; M::Tuple{Int,Int}, dimorder::Int64, qrule=WavePropBase.Fejer)
    f = (u) -> xe + [u[1], u[2], xe[3]].*(xe - xs)
    d = Nystrom.HyperRectangle((0.,0.),(1.,1.))
    Γ = Domain(Nystrom.ParametricEntity(f,d))
    msh = Nystrom.meshgen(Γ,M)
    NystromMesh(msh,Γ, WavePropBase.TensorProductQuadrature(qrule(dimorder),qrule(dimorder)))
end
HorizontalStraightPlane(xe,xs;M,dimorder,qrule) = HorizontalStraightPlane(SVector(xe),SVector(xs);M,dimorder,qrule)
HorizontalStraightPlane(xe,xs;M,dimorder) = HorizontalStraightPlane(SVector(xe),SVector(xs);M,dimorder)

""" ndofs returns the number of dofs to be used for each object to be discretized """

function ndofs(P::Problem{2,1},Fig::Obstacle, w::Window , ppw::Any)
    λ₁,λ₂ = 2π./[P.pde[1].k, P.pde[2].k]
    Mv = Int(ceil( 2*w.A/λ₁ * ppw ))
    Ms = Int(ceil( 2π*Fig.radius / max(λ₁,λ₂) * ppw ))
    Ms,Mv
end
function ndofs(P::Problem{3,2},Fig::Obstacle, w::Window , ppw::Any)
    λ₁,λ₂ = 2π./[P.pde[1].k, P.pde[2].k]
    N1 = Int(ceil( 4π*Fig.radius^2 / max(λ₁,λ₂)^2 * ppw ))
    Nx = (Int(ceil(P.L[1]/λ₁ * ppw)), Int(ceil(2*w.A/λ₁ * ppw)))
    Ny = (Int(ceil(P.L[2]/λ₁ * ppw)), Int(ceil(2*w.A/λ₁ * ppw)))
    N1,Nx,Ny
end
function ndofs(P::Problem{3,1},Fig::Obstacle, w::Window , ppw::Any)
    λ₁,λ₂ = 2π./[P.pde[1].k, P.pde[2].k]
    Ns = Int(ceil( 4π*Fig.radius^2 / max(λ₁,λ₂)^2 * ppw ))
    Nv = Int(ceil(2*w.A/λ₁ * ppw))
    Ns,Nv
end

""" unitcell creates a unit cell, thus contains no displaced curves and a single osbtacle """

function unitcell(P::Problem{2,1}, Fig::Obstacle, w::Window ; ppw::Any,dimorder::Int)
    Ms,Mv = ndofs(P,Fig,w,ppw)
    Γ₁ = Scatterer(Ms, Fig; c = (0.,0.), dimorder = dimorder)
    Γ₂ = StraightLine((-0.5*P.L,-w.A),(-0.5*P.L,w.A); M = Mv, dimorder = dimorder)
    Γ₃ = StraightLine((+0.5*P.L,-w.A),(+0.5*P.L,w.A); M = Mv, dimorder = dimorder)
    return [Γ₁,Γ₂,Γ₃]
end
function unitcell(P::Problem{3,2},Fig::Obstacle, w::Window ; ppw::Any,dimorder::Int)
    N1,Nx,Ny = ndofs(P,Fig,w,ppw)
    Γ₁ = Scatterer(N1,Fig; c = (0.,0.,0.), dimorder = dimorder)
    Γ₂ = StraightPlane((-0.5P.L[1],-0.5P.L[2],-w.A),(-0.5P.L[1],+0.5P.L[2],+w.A); M=Nx, dimorder = dimorder)  
    Γ₃ = StraightPlane((+0.5P.L[1],-0.5P.L[2],-w.A),(+0.5P.L[1],+0.5P.L[2],+w.A); M=Nx, dimorder = dimorder)
    Γ₄ = StraightPlane((+0.5P.L[1],-0.5P.L[2],-w.A),(-0.5P.L[1],-0.5P.L[2],+w.A); M=Ny, dimorder = dimorder)
    Γ₅ = StraightPlane((+0.5P.L[1],+0.5P.L[2],-w.A),(-0.5P.L[1],+0.5P.L[2],+w.A); M=Ny, dimorder = dimorder)
    @assert WavePropBase.normal(Γ₂.dofs[1]) == [1.,0.,0.]
    @assert WavePropBase.normal(Γ₃.dofs[1]) == [1.,0.,0.]
    @assert WavePropBase.normal(Γ₄.dofs[1]) == [0.,1.,0.]
    @assert WavePropBase.normal(Γ₅.dofs[1]) == [0.,1.,0.]
    return [Γ₁,Γ₂,Γ₃,Γ₄,Γ₅]
end
function unitcell(P::Problem{3,1},Fig::Obstacle, w::Window; ppw,dimorder::Int)
    Ns,Nv = ndofs(P,Fig,w,ppw)
    Γ₁ = Scatterer(Ns,Fig; c = (0.,0.,0.), dimorder = dimorder)
    Γ₂ = StraightPlane((-0.5w.A,-w.A,-w.A),(-0.5w.A,+w.A,+w.A); M=(Nv,Nv), dimorder = dimorder)  
    Γ₃ = StraightPlane((+0.5w.A,-w.A,-w.A),(+0.5w.A,+w.A,+w.A); M=(Nv,Nv), dimorder = dimorder)
    @assert WavePropBase.normal(Γ₂.dofs[1]) == [1.,0.,0.]
    @assert WavePropBase.normal(Γ₃.dofs[1]) == [1.,0.,0.]
    return [Γ₁,Γ₂,Γ₃]
end

""" extendedcell is used to avoid singularities
    - 2D1D returns a triple-cell setup
    - 3D2D returns a nine-cell setup """

function extendedcell(P::Problem{2,1}, Fig::Obstacle, w::Window ;ppw::Any,dimorder::Int)
    Ms,Mv = ndofs(P,Fig,w,ppw)
    Γ₁ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = Scatterer(Ms, Fig; c = (i*P.L,0.), dimorder = dimorder)
        merge!(Γ₁,Dict(i => Γ))
    end
    Γ₂ = StraightLine((-1.5*P.L,-w.A),(-1.5*P.L,w.A); M = Mv, dimorder = dimorder)
    @assert WavePropBase.normal(Γ₂.dofs[1]) == [1.,0.]
    Γ₃ = StraightLine((+1.5*P.L,-w.A),(+1.5*P.L,w.A); M = Mv, dimorder = dimorder)
    @assert WavePropBase.normal(Γ₃.dofs[1]) == [1.,0.]
    return [Γ₁, Γ₂, Γ₃]
end
function extendedcell(P::Problem{3,2},Fig::Obstacle, w::Window ; ppw::Any,dimorder::Int)
    N1,Nx,Ny = ndofs(P,Fig,w,ppw)
    Γ₁ = Dict{Tuple{Int,Int},NystromMesh}()
    for i = -1:1, j=-1:1
        Γ = Scatterer(N1, Fig; c = (i*P.L[1],j*P.L[2],0.), dimorder = dimorder)
        merge!(Γ₁, Dict((i,j) => Γ))
    end
    Γ₂ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = StraightPlane((-1.5P.L[1],(-0.5+i)P.L[2],-w.A),(-1.5P.L[1],(+0.5+i)P.L[2],+w.A); M=Nx, dimorder = dimorder)
        @assert WavePropBase.normal(Γ.dofs[1]) == [1.,0.,0.]
        merge!(Γ₂, Dict(i => Γ))
    end
    Γ₃ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = StraightPlane((+1.5P.L[1],(-0.5+i)P.L[2],-w.A),(+1.5P.L[1],(+0.5+i)P.L[2],+w.A); M=Nx, dimorder = dimorder)
        @assert WavePropBase.normal(Γ.dofs[1]) == [1.,0.,0.]
        merge!(Γ₃, Dict(i => Γ))
    end
    Γ₄ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = StraightPlane(((+0.5+i)P.L[1],-1.5P.L[2],-w.A),((-0.5+i)P.L[1],-1.5P.L[2],+w.A); M=Ny, dimorder = dimorder)
        @assert WavePropBase.normal(Γ.dofs[1]) == [0.,1.,0.]
        merge!(Γ₄, Dict(i => Γ))
    end
    Γ₅ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = StraightPlane(((+0.5+i)P.L[1],+1.5P.L[2],-w.A),((-0.5+i)P.L[1],+1.5P.L[2],+w.A); M=Ny, dimorder = dimorder)
        @assert WavePropBase.normal(Γ.dofs[1]) == [0.,1.,0.]
        merge!(Γ₅, Dict(i => Γ))
    end
    return [Γ₁,Γ₂,Γ₃,Γ₄,Γ₅]
end
function extendedcell(P::Problem{3,1}, Fig::Obstacle, w::Window ;ppw::Any,dimorder::Int)
    Ns,Nv = ndofs(P,Fig,w,ppw)
    Γ₁ = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ = Scatterer(Ns, Fig; c = (i*P.L,0.,0.), dimorder = dimorder)
        merge!(Γ₁,Dict(i => Γ))
    end
    Γ₂ = StraightPlane((-1.5*P.L,-w.A,-w.A),(-1.5*P.L,w.A,w.A); M = (Nv,Nv), dimorder = dimorder)
    @assert WavePropBase.normal(Γ₂.dofs[1]) == [1.,0.,0.]
    Γ₃ = StraightPlane((+1.5*P.L,-w.A,-w.A),(+1.5*P.L,w.A,w.A); M = (Nv,Nv), dimorder = dimorder)
    @assert WavePropBase.normal(Γ₃.dofs[1]) == [1.,0.,0.]
    return [Γ₁, Γ₂, Γ₃]
end

""" meshplot are used to plot NystromMesh dofs"""
function meshplot!(Γ::NystromMesh{3,T,M,NM}) where {T,M,NM}
    x = [q.coords[1] for q in Γ.dofs]
    y = [q.coords[2] for q in Γ.dofs]
    z = [q.coords[3] for q in Γ.dofs]
    Plots.scatter!(x,y,z, legend = false, aspect_ratio = 1)
end
function meshplot(Γ::NystromMesh{3,T,M,NM}) where {T,M,NM}
    x = [q.coords[1] for q in Γ.dofs]
    y = [q.coords[2] for q in Γ.dofs]
    z = [q.coords[3] for q in Γ.dofs]
    Plots.scatter(x,y,z, legend = false, aspect_ratio = 1)
end


