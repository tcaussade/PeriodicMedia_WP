""" returns αₙ and βₙ for a given n"""

function seriesconstant(P::Problem{2,1}, n::Int)
    αₙ = P.dir[1]*P.pde[1].k + 2π*n/P.L
    βₙ = sqrt(Complex(P.pde[1].k^2-αₙ^2))
    return αₙ, βₙ
end
function seriesconstant(P::Problem{3,2}, n::Tuple{Int,Int})
    αₙ₁ = P.dir[1]*P.pde[1].k + 2π*n[1]/P.L[1]
    αₙ₂ = P.dir[2]*P.pde[1].k + 2π*n[2]/P.L[2]
    βₙ  = sqrt(complex(P.pde[1].k^2-αₙ₁^2-αₙ₂^2))
    return αₙ₁, αₙ₂, βₙ
end


"""
    Field corrections are computed
"""

function fixdelta(P::Problem{N,1}, tol::Float64) where N
    δ = Int[]
    nC= 10
    for n = -nC:nC
        _,βₙ = seriesconstant(P,n)
        if abs(βₙ) < tol push!(δ,n) end
    end
    return δ
end

function fixdelta(P::Problem{3,2}, tol::Float64)
    δ = Tuple{Int,Int}[]
    nC = 10
    for n1 = -nC:nC, n2 = -nC:nC
        _,_,βₙ = seriesconstant(P,(n1,n2))
        if abs(βₙ) < tol push!(δ,(n1,n2)) end
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
        _,βₙ = seriesconstant(P,n)

        Ψₙ⁺ = Ψₙ(P,G,r₁,n₁; n=n, sgn = +1.)
        Ψₙ⁻ = Ψₙ(P,G,r₁,n₁; n=n, sgn = -1.)
        Lₙ⁺ = Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.)
        Lₙ⁻ = Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.)

        if abs(βₙ) < 1e-8
            @info "Near RW anomaly" βₙ n
            dΨₙ⁺ = dΨₙ(P,G,r₁,n₁; n=n, sgn = +1.)
            dΨₙ⁻ = dΨₙ(P,G,r₁,n₁; n=n, sgn = -1.)
            dLₙ⁺ = dLₙ(P,x⁺,R⁺ ; n=n, sgn = +1.)
            dLₙ⁻ = dLₙ(P,x⁻,R⁻ ; n=n, sgn = -1.)
            FRO += 1/(2*im) * (dΨₙ⁻*Lₙ⁻ + Ψₙ⁻*dLₙ⁻ - dΨₙ⁺*Lₙ⁺ - Ψₙ⁺*dLₙ⁺)
        else
            FRO += exp(im*βₙ*H)/(2*im*βₙ) * (Ψₙ⁻ * Lₙ⁻ - Ψₙ⁺ * Lₙ⁺)
        end
    end
    return FRO
end

function finiterankoperator(P::Problem{3,2},G::Vector{NystromMesh{N,T,M,NM}},Gt::Vector; δ::Float64,H::Float64) where {N,T,M,NM}
    Cδ = fixdelta(P, δ)
    r₁,n₁ = [q.coords for q in G[1].dofs], [WavePropBase.normal(q) for q in G[1].dofs]
    dR⁺, R⁺, x⁺ = planegradscat(P,Gt; H = +H)
    dR⁻, R⁻, x⁻ = planegradscat(P,Gt; H = -H)
    len = 2*length(G[1].dofs)+2*length(G[2].dofs)+2*length(G[4].dofs)
    FRO = zeros(ComplexF64,len,len)
    for n in Cδ
        _,_,βₙ = seriesconstant(P,(n1,n2))

        Ψₙ⁺ = Ψₙ(P,G,r₁,n₁; n=n, sgn = +1.)
        Ψₙ⁻ = Ψₙ(P,G,r₁,n₁; n=n, sgn = -1.)
        Lₙ⁺ = Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.)
        Lₙ⁻ = Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.)

        if abs(βₙ) < 1e-8
            @info "Near RW anomaly" βₙ n
            # error("not implemented")
            @warn "not implemented"
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
        if abs(βₙ) < 1e-8
            Lₙ⁺  = Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.0) * σw
            Lₙ⁻  = Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.0) * σw
            dLₙ⁺ = dLₙ(P,x⁺,R⁺ ; n=n, sgn = +1.0) * σw
            dLₙ⁻ = dLₙ(P,x⁻,R⁻ ; n=n, sgn = -1.0) * σw
            y = [r[2] for r in eval]
            duₙ⁺ = +im*y.*uₙ⁺
            duₙ⁻ = -im*y.*uₙ⁻
            fix += 1/(2*im) * (duₙ⁻*Lₙ⁺ + uₙ⁻*dLₙ⁺ - duₙ⁺*Lₙ⁻ - uₙ⁺*dLₙ⁻)
        else
            Cₙ⁺ = +exp(im*βₙ*H)/(2*im*βₙ) * Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.0) * σw
            Cₙ⁻ = -exp(im*βₙ*H)/(2*im*βₙ) * Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.0) * σw
            fix += Cₙ⁺ * uₙ⁻ + Cₙ⁻ * uₙ⁺
        end
    end
    return fix
end

function scatcorrection(P::Problem{3,2},G::Vector,eval,σw::Vector{ComplexF64}; H::Float64, δ::Float64)
    @assert G[1] isa Dict{Tuple{Int,Int}} # @assert G isa extendedcell
    Cδ = fixdelta(P,δ)
    dR⁺, R⁺, x⁺ = planegradscat(P,G; H = +H)
    dR⁻, R⁻, x⁻ = planegradscat(P,G; H = -H)
    fix = zeros(length(eval))
    for n in Cδ
        # αₙ₁ = P.dir[1]*P.pde[1].k + 2π*n[1]/P.L[1]
        # αₙ₂ = P.dir[2]*P.pde[1].k + 2π*n[2]/P.L[2]
        # βₙ  = sqrt(complex(P.pde[1].k^2-αₙ₁^2-αₙ₂^2))
        αₙ₁,αₙ₂,βₙ = seriesconstant(P,(n1,n2))
        uₙ⁺ = [exp.(im*αₙ₁*r[1] + im*αₙ₂*r[2] + im*βₙ*r[3]) for r in eval]
        uₙ⁻ = [exp.(im*αₙ₂*r[1] + im*αₙ₂*r[2] - im*βₙ*r[3]) for r in eval]
        if abs(βₙ) < 1e-8
            error("RW anomaly")
        else
            # Check results
            Cₙ⁺ = +exp(im*βₙ*H)/(2*im*βₙ) * Lₙ(P,x⁺,R⁺,dR⁺; n=n, sgn = +1.0) * σw
            Cₙ⁻ = -exp(im*βₙ*H)/(2*im*βₙ) * Lₙ(P,x⁻,R⁻,dR⁻; n=n, sgn = -1.0) * σw
            fix += Cₙ⁺ * uₙ⁻ + Cₙ⁻ * uₙ⁺
        end
    end
    return fix
end

"""
    compute the gradscat and scat field over a horizontal curve/plane
"""

function linegradscat(P::Problem{2,1},G::Vector; H::Float64)
    lin = StraightLine((+0.5*P.L,H),(-0.5*P.L,H); M = 50, dimorder = 5, qrule = WavePropBase.TrapezoidalOpen)
    x   = [q.coords[1] for q in lin.dofs]
    return gradscatpotential(P,lin,G), scatpotential(P,lin,G), x
end

function planegradscat(P::Problem{3,2}, G::Vector; H::Float64)
    Trap = WavePropBase.TrapezoidalOpen
    Σ = HorizontalStraightPlane((-0.5*P.L[1],-0.5*P.L[2],H),(0.5*P.L[1],0.5*P.L[2],H); M = (20,20), dimorder = 3, qrule = Trap)
    x = [q.coords[1:2] for q in Σ.dofs]
    return gradscatpotential(P,Σ,G), scatpotential(P,Σ,G), x
end

""" Compute Lₙ and dLₙ as a row vector """

function Lₙ(P::Problem{2,1},x, R::Matrix{ComplexF64},dR::Matrix{ComplexF64}; n::Int, sgn) 
    # to integrate you must multiply by Lₙ
    αₙ,βₙ = seriesconstant(P,n)
    exp.(im*αₙ*x)' * (dR - sign(sgn) * im*βₙ* R) * P.L/length(x) /P.L
end
function Lₙ(P::Problem{3,2}, x, R::Matrix{ComplexF64},dR::Matrix{ComplexF64}; n::Tuple{Int64, Int64}, sgn)
    _,_,βₙ = seriesconstant(P,(n1,n2))
    x1  = [q[1] for q in x]
    x2  = [q[2] for q in x]
    exp.(im*αₙ₁*x1+im*αₙ₂*x2)' * (dR - sign(sgn)*im*βₙ*R) * P.L[1]*P.L[2]/prod(size(x)) /P.L[1]/P.L[2]
end

function dLₙ(P::Problem{2,1},x,R::Matrix{ComplexF64}; n::Int, sgn)
    αₙ,_ = seriesconstant(P,n)
    exp.(im*αₙ*x)' * (-im*sign(sgn)*R) * P.L/length(x) /P.L
end

""" Compute Ψₙ and dΨₙ as a column vector """

function Ψₙ(P::Problem{2,1},G::Vector,r1,n1; n::Int, sgn)
    αₙ,βₙ = seriesconstant(P,n)
    uₙ = [exp(im*(αₙ*r[1] - sign(sgn)*βₙ*r[2])) for r in r1] # uₙ|Γ₁
    duₙ = im*[(αₙ*ν[1] - sign(sgn)*βₙ*ν[2]) for ν in n1] .* uₙ # duₙ|Γ₁
    [uₙ; duₙ; zeros(2*length(G[2].dofs))]
end

function Ψₙ(P::Problem{3,2},G::Vector,r1,n1; n::Tuple{Int,Int}, sgn)
    αₙ₁,αₙ₂,βₙ = seriesconstant(P,(n1,n2))
    uₙ = [exp(im*(αₙ₁*r[1]+αₙ₂*r[2] - sign(sgn)*βₙ*r[3])) for r in r1] # uₙ|Γ₁
    duₙ = im*[(αₙ₁*ν[1]+αₙ₂*ν[2] - sign(sgn)*βₙ*ν[3]) for ν in n1] .* uₙ # duₙ|Γ₁
    [uₙ; duₙ; zeros(2*length(G[2].dofs)+2*length(G[4].dofs))]
end

function dΨₙ(P::Problem{2,1},G::Vector,r1,n1; n::Int, sgn)
    αₙ,_ = seriesconstant(P,n)
    uₙ = [-sign(sgn) * im*r[2] * exp(im*αₙ*r[1]) for r in r1] # uₙ|Γ₁
    vₙ = [-sign(sgn) * (-r1[j][2]*αₙ*n1[j][1] + im*n1[j][2]) * exp(im*αₙ*r1[j][1]) for j=1:length(r1)] # vₙ|Γ₁ #RE-DO
    [uₙ; vₙ; zeros(2*length(G[2].dofs))]
end



