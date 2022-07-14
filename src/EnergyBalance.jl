"""
    function for energy balance assesement
"""

function energytest(P::Problem{2,1},G::Vector,w::Window,σ::Vector{ComplexF64}; FRO = true)
    @assert G[1] isa Dict
    H = w.A*w.c
    Trap = WavePropBase.TrapezoidalOpen
    upp = StraightLine((0.5*P.L,+H),(-0.5*P.L,+H); M = 100, dimorder = 5, qrule = Trap)
    low = StraightLine((0.5*P.L,-H),(-0.5*P.L,-H); M = 100, dimorder = 5, qrule = Trap)
    xi  = [q.coords[1] for q in upp.dofs]

    W = wgfmatrix(G,w)
    u1 = scatpotential(P,upp,G)*W*σ
    u2 = scatpotential(P,low,G)*W*σ

    if FRO 
        # u1 += scatfiniterankoperator()
        # u2 += scatfiniterankoperator()
    else
        @info "Non corrected energy balance test"
    end

    energytest(P,H,u1,u2,xi)
end

function energytest(P::Problem{2,1},H::Float64,u1::Vector{ComplexF64},u2::Vector{ComplexF64}, x::Vector{Float64})
    nC = 20
    S = zero(ComplexF64)
    α,β = P.pde[1].k*P.dir
    for n=-nC:nC
        αₙ = α + n*2π/P.L
        if P.pde[1].k ≥ abs(αₙ)
            βₙ = sqrt(complex(P.pde[1].k^2-αₙ^2))
            B⁺ = exp(-im*βₙ*H) * P.L/length(x) * dot(exp.(im*αₙ*x), u1)/P.L
            B⁻ = exp(-im*βₙ*H) * P.L/length(x) * dot(exp.(im*αₙ*x), u2)/P.L
            if n==0
                global B₀ = B⁻
            end
            S = S + βₙ/β * (abs(B⁺)^2+abs(B⁻)^2)
        end
    end
    abs(2*real(B₀)+S)
end

function energytest(P::Problem{3,2},G::Vector,w::Window,σ::Vector{ComplexF64}; FRO = true)
    H = w.A*w.c
    Trap = WavePropBase.TrapezoidalOpen
    upp = HorizontalStraightPlane((-0.5*P.L[1],-0.5*P.L[2],+H),(0.5*P.L[1],0.5*P.L[2],+H); M = (20,20), dimorder = 3, qrule = Trap)
    low = HorizontalStraightPlane((-0.5*P.L[1],-0.5*P.L[2],-H),(0.5*P.L[1],0.5*P.L[2],-H); M = (20,20), dimorder = 3, qrule = Trap)
    xi  = [q.coords[1] for q in upp.dofs]
    yi  = [q.coords[2] for q in upp.dofs]

    W = wgfmatrix(G,w)
    u1 = scatpotential(P,upp,G)*W*σ
    u2 = scatpotential(P,low,G)*W*σ

    energytest(P,H,u1,u2,xi,yi)
end

function energytest(P::Problem{3,2},H::Float64,u1::Vector{ComplexF64},u2::Vector{ComplexF64},xi::Vector{Float64},yi::Vector{Float64})
    S = 0.0+im*0
    nC= 10
    α₁,α₂,β = P.pde[1].k*P.dir
    for n1 = -nC:nC, n2 = -nC:nC
        αₘ₁ = α₁ + n1*2π/P.L[1]
        αₘ₂ = α₂ + n2*2π/P.L[2]
        if P.pde[1].k^2 ≥ abs(αₘ₁)^2 + abs(αₘ₂)^2
            βₙ = sqrt((P.pde[1].k^2 - αₘ₁^2 - αₘ₂^2 ))
            @show βₙ, (n1,n2)
            B⁺ = exp(-im*βₙ*H)/P.L[1]/P.L[2] * sum( u1.*exp.(-im*αₘ₁*xi-im*αₘ₂*yi) ) *P.L[1]*P.L[2]/length(u1)
            B⁻ = exp(-im*βₙ*H)/P.L[1]/P.L[2] * sum( u2.*exp.(-im*αₘ₁*xi-im*αₘ₂*yi) ) *P.L[1]*P.L[2]/length(u1)
            if n1 == 0 && n2 == 0
                global B₀ = B⁻
            end
            S += βₙ/β * (abs(B⁺)^2+abs(B⁻)^2)
        end
    end
    abs(2*real(B₀)+S)
end