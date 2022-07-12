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
        u1 += scatfiniterankoperator()
        u2 += scatfiniterankoperator()
    else
        @info "Non corrected energy balance test"
    end

    energytest(P,H,u1,u2,xi)
end

function energytest(P::Problem{2,1},H::Float64,u1::Vector{ComplexF64},u2::Vector{ComplexF64}, x::Vector{Float64})
    nC = 20
    S = zero(ComplexF64)
    α,β = P.pde[1].k.*P.dir
    for n=-nC:nC
        αₙ = α + n*2π/P.L
        if P.pde[1].k ≥ abs(αₙ)
            βₙ = sqrt(complex(P.pde[1].k^2-αₙ^2))
            a⁺ = P.L/length(x) * dot(exp.(im*αₙ*x), u1)/P.L; B⁺ = exp(-im*βₙ*H) * a⁺
            a⁻ = P.L/length(x) * dot(exp.(im*αₙ*x), u2)/P.L; B⁻ = exp(-im*βₙ*H) * a⁻
            if n==0
                global B₀ = B⁻
            end
            S = S + βₙ/β * (abs(B⁺)^2+abs(B⁻)^2)
        end
    end
    return abs(2*real(B₀)+S) #Energy
end