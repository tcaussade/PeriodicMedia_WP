"""
    Field corrections are computed
"""

function fixdelta(P::Problem{N,1}, tol::Float64) where N
    δ = Int[]
    nC= 10
    α = P.dir[1]*P.pde[1].k
    for n = -nC:nC
        αₙ = α + 2π*n/P.L
        βₙ = sqrt(Complex(P.pde[1].k^2-αₙ^2))
        if abs(βₙ) < tol push!(δ,n) end
    end
    return δ
end

"""
    finiterankoperator() modifies the linear system to be solved accordingly
"""

function finiterankoperator(P::Problem{2,1},G::Vector{NystromMesh{N,T,M,NM}},Gt::Vector; δ::Float64,H::Float64) where {N,T,M,NM}
    Cδ = fixdelta(P, δ)
    r₁,n₁ = [q.coords for q in G[1].dofs], [WavePropBase.normal(q) for q in G[1].dofs]
    dR⁺, R⁺, x⁺ = linegradscat(P,Gt; H = +H)
    dR⁻, R⁻, x⁻ = linegradscat(P,Gt; H = -H)

    len = 2*length(G[1].dofs)+2*length(G[2].dofs)
    FRO = zeros(ComplexF64,len,len)
    for n in Cδ
        αₙ = P.dir[1]*P.pde[1].k + 2π*n/P.L
        βₙ = sqrt(Complex(P.pde[1].k^2-αₙ^2))

        Ψₙ⁺ = Ψₙ(P,G,r₁,n₁; n=n, sgn = +1.)
        Ψₙ⁻ = Ψₙ(P,G,r₁,n₁; n=n, sgn = -1.)
        Lₙ⁺ = Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.)
        Lₙ⁻ = Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.)

        if abs(βₙ) < 1e-4
            @warn "Near RW anomaly"
            continue #TO DO
        else
            FRO += exp(im*βₙ*H)/(2*im*βₙ) * (Ψₙ⁻ * Lₙ⁻ - Ψₙ⁺ * Lₙ⁺)
        end
    end
    return FRO
end


"""
    scatfiniterankoperator() modifies the scattered field evaluation
"""

function scatcorrection(P::Problem{2,1},G::Vector,eval,σw::Vector{ComplexF64}; H::Float64, δ::Float64)
    @assert G[1] isa Dict # @assert G isa extendedcell
    Cδ = fixdelta(P, δ)
    dR⁺, R⁺, x⁺ = linegradscat(P,G; H = +H)
    dR⁻, R⁻, x⁻ = linegradscat(P,G; H = -H)
    fix = zeros(length(eval))
    for n in Cδ
        αₙ  = P.dir[1]*P.pde[1].k + 2π*n/P.L
        βₙ  = sqrt(complex(P.pde[1].k^2-αₙ^2))
        uₙ⁺ = [exp.(im*αₙ*r[1] + im*βₙ*r[2]) for r in eval]
        uₙ⁻ = [exp.(im*αₙ*r[1] - im*βₙ*r[2]) for r in eval]
        Cₙ⁺ = +exp(im*βₙ*H)/(2*im*βₙ) * Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.0) * σw
        Cₙ⁻ = -exp(im*βₙ*H)/(2*im*βₙ) * Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.0) * σw
        fix += Cₙ⁺ * uₙ⁻ + Cₙ⁻ * uₙ⁺
    end
    return fix
end
# Finish

function linegradscat(P::Problem{2,1},G::Vector; H::Float64)
    lin = StraightLine((+0.5*P.L,H),(-0.5*P.L,H); M = 50, dimorder = 5, qrule = WavePropBase.TrapezoidalOpen)
    x   = [q.coords[1] for q in lin.dofs]
    return gradscatpotential(P,lin,G), scatpotential(P,lin,G), x
end

function Lₙ(P::Problem{2,1},x, R,dR; n, sgn) 
    # to integrate you must multiply by Lₙ
    αₙ = P.dir[1]*P.pde[1].k + 2π*n/P.L
    βₙ = sqrt(Complex(P.pde[1].k^2-αₙ^2))
    exp.(im*αₙ*x)' * (dR - sign(sgn) * im*βₙ* R) * P.L/length(x) /P.L
end

function Ψₙ(P,G,r1,n1; n, sgn)
    αₙ = P.dir[1]*P.pde[1].k + 2π*n/P.L
    βₙ = sqrt(Complex(P.pde[1].k^2-αₙ^2))
    uₙ = [exp(im*(αₙ*r[1] - sign(sgn)*βₙ*r[2])) for r in r1] # uₙ|Γ₁
    duₙ = im*[(αₙ*ν[1] - sign(sgn)*βₙ*ν[2]) for ν in n1] .* uₙ # duₙ|Γ₁
    [uₙ; duₙ; zeros(2*length(G[2].dofs))]
end



