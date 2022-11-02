"""
    function for energy balance assesement
"""

""" 
    Build horizontal curves/surfaces
"""

function testboundary(P::Problem{2,1}; h::Float64)
    # @assert G[1] isa Dict # @assert G isa extendedcell
    Trap = WavePropBase.TrapezoidalOpen
    upp = StraightLine((+0.5*P.L,+h),(-0.5*P.L,+h); M = 10, dimorder = 5, qrule = Trap)
    low = StraightLine((+0.5*P.L,-h),(-0.5*P.L,-h); M = 10, dimorder = 5, qrule = Trap)
    xi  = [q.coords[1] for q in upp.dofs]
    return upp,low, xi
end
function testboundary(P::Problem{3,2}; h::Float64)
    # Trap = WavePropBase.TrapezoidalOpen
    upp = HorizontalStraightPlane((-0.5*P.L[1],-0.5*P.L[2],+h),(0.5*P.L[1],0.5*P.L[2],+h); M = (20,20), dimorder = 4)
    low = HorizontalStraightPlane((-0.5*P.L[1],-0.5*P.L[2],-h),(0.5*P.L[1],0.5*P.L[2],-h); M = (20,20), dimorder = 4)
    xi  = [q.coords[1] for q in upp.dofs]
    yi  = [q.coords[2] for q in upp.dofs]
    return upp,low,(xi,yi)
end
function testboundary(P::Problem{3,1}; h::Float64)
    @warn "not implemented"
end

"""
    Evaluate the scattered field
    FRO = Finite Rank Operator
"""

function energytest(P::Problem{N,NP},G::Vector,w::Window,σ::Vector{ComplexF64}; FRO::Bool, h::Float64) where {N,NP}
    upp,low,ri = testboundary(P; h = h)

    σw = lmul!(wgfmatrix(G,w),σ)
    u1 = scatpotential(P,upp,G)*σw
    u2 = scatpotential(P,low,G)*σw

    if FRO 
        u1 += scatcorrection(P,G,[q.coords for q in upp.dofs],σw; h=+h, δ = 0.75*P.pde[1].k)
        u2 += scatcorrection(P,G,[q.coords for q in low.dofs],σw; h=+h, δ = 0.75*P.pde[1].k)
    # else
        # @info "Non corrected energy balance test"
    end

    energytest(P,h,u1,u2,ri)
end

""" computes coefficients and assets energy conservation """
function energytest(P::Problem{2,1},h::Float64,u1::Vector{ComplexF64},u2::Vector{ComplexF64}, x::Vector{Float64})
    nC = 30
    S = zero(ComplexF64)
    R = zero(ComplexF64)
    # T = zero(ComplexF64)
    α,β = P.pde[1].k*P.dir
    for n=-nC:nC
        αₙ = α + n*2π/P.L
        if P.pde[1].k ≥ abs(αₙ)
            βₙ = sqrt((P.pde[1].k^2-αₙ^2))
            B⁺ = exp(-im*βₙ*h) * P.L/length(x) * dot(exp.(im*αₙ*x), u1)/P.L
            B⁻ = exp(-im*βₙ*h) * P.L/length(x) * dot(exp.(im*αₙ*x), u2)/P.L
            if n==0
                global B₀⁻ = B⁻
            end
            S += βₙ/β * (abs(B⁺)^2+abs(B⁻)^2)
            R += βₙ/β * abs(B⁺)^2
            # T += βₙ/β * abs(B⁻)^2
        end
    end
    e = abs(2*real(B₀⁻)+S)
    T = 1 + 2*real(B₀⁻) + S - R
    # abs(2*real(B₀)+S)
    return e,R,T
end
function energytest(P::Problem{3,2},h::Float64,u1::Vector{ComplexF64},u2::Vector{ComplexF64},ri)
    xi,yi = ri
    nC= 10
    S = zero(ComplexF64)
    R = zero(ComplexF64)
    # T = zero(ComplexF64)
    α₁,α₂,β = P.pde[1].k*P.dir
    for n1 = -nC:nC, n2 = -nC:nC
        αₘ₁ = α₁ + n1*2π/P.L[1]
        αₘ₂ = α₂ + n2*2π/P.L[2]
        if P.pde[1].k^2 ≥ abs(αₘ₁)^2 + abs(αₘ₂)^2
            βₙ = sqrt(complex(P.pde[1].k^2 - αₘ₁^2 - αₘ₂^2 ))
            @show B⁺ = exp(-im*βₙ*h)/P.L[1]/P.L[2] * sum( u1.*exp.(-im*αₘ₁*xi-im*αₘ₂*yi) ) *P.L[1]*P.L[2]/length(u1)
            @show B⁻ = exp(-im*βₙ*h)/P.L[1]/P.L[2] * sum( u2.*exp.(-im*αₘ₁*xi-im*αₘ₂*yi) ) *P.L[1]*P.L[2]/length(u1)
            if n1 == 0 && n2 == 0
                global B₀⁻ = B⁻
            end
            S += βₙ/β * (abs(B⁺)^2+abs(B⁻)^2)
            R += βₙ/β * abs(B⁺)^2
            # T += βₙ/β * abs(B⁻)^2
        end
    end
    e = abs(2*real(B₀⁻)+S)
    T = 1 + 2*real(B₀⁻) + S - R
    # abs(2*real(B₀)+S)
    return e,R,T
end
function energytest(P::Problem{3,1},ρ::Float64,u::Vector{ComplexF64},ri)
    r,θ = ri
    α,β₁,β₂ = P.pde[1].k * P.dir
    θtilda = atan(β₂,β₁)
    nC = 20
    S = zero(ComplexF64)
    for n in -nC:nC
        αₙ = α + n*2π/P.L
        if P.pde[1].k ≥ abs(αₙ)
            βₙ = sqrt((P.pde[1].k^2-αₙ^2))
            for l in -nC:nC
                uₘˡ =  1/(2π*P.L)/hankelh1(l,βₙ*ρ) * sum( u.* exp.(-im*l*θ-im*αₙ*r) ) * 2π*P.L/length(u1)
                S += abs(uₘˡ)^2
                if n == 0
                    S₀ += im*exp(im*l*θtilda-π/2.) * uₘˡ
                end
            end
        end
    end
    e = imag(S₀) + S
    return e
end