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

function triplecell(P::Problem{2,1}, Fig::Obstacle, w::Window ;ppw::Int,dimorder::Int)
    Mv = Int(ceil( 2*w.A*(P.pde[1].k/2π) * ppw ))
    Ms = Int(ceil( 2π*Fig.radius*(P.pde[1].k/2π) * ppw ))

    loctoscat = Dict{Int,NystromMesh}()
    for i = -1:1
        Γ₁ = Scatterer(Ms, Fig; c = (i*P.L,0.), dimorder = dimorder)
        merge!(loctoscat,Dict(i => Γ₁))
    end

    Γ₂ = StraightLine((-1.5*P.L,-w.A),(-1.5*P.L,w.A); M = Mv, dimorder = dimorder)
    Γ₃ = StraightLine((+1.5*P.L,-w.A),(+1.5*P.L,w.A); M = Mv, dimorder = dimorder)

    return [loctoscat, Γ₂, Γ₃]
end