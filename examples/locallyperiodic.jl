using PeriodicMedia

function locallyperiodic(obstaclesizes)
    
    # set physical params
    k1 = 9.5 
    k = [k1,17.]
    θ = [π/4.,π/4.]
    L = [0.5, 0.5]
    P = Problem(k,θ,L; ambdim = 3, geodim = 2)
    WGF = Window(0.5,20*(2π/k[1]))

    # Mesh params
    ppw, dim = 2, 1
    Sphere = PeriodicMedia.ParametricSurfaces.Sphere

    amps = Vector{Float64}(undef, 0)
    for ratio in obstaclesizes
        #create geometry
        PeriodicMedia.clear_entities!()
        Fig = Obstacle(Sphere, ratio * L)
        Γs = unitcell(P,Fig, WGF; ppw = ppw, dimorder = dim)
        # @show smat = length(Γs[1].dofs)*2+length(Γs[2].dofs)*2+length(Γs[4].dofs)*2
        Γt = extendedcell(P,Fig, WGF; ppw = ppw, dimorder = dim)

        ϕ = solver(P,Γs,Γt,WGF; FRO = true)
        @show eb = energytest(P,Γt,WGF, ϕ; FRO = true, H = 1.0)

        # Compute transmitted field
        @show coef = transmitted(P,Γt,WGF, ϕ; FRO = true, H = 1.0)
        push!(amps, coef)
    end
end

function transmitted(P::Problem{3,2},G::Vector,w::Window,σ::Vector{ComplexF64}; FRO = true, H::Float64)
    Trap = WavePropBase.TrapezoidalOpen
    upp = HorizontalStraightPlane((-0.5*P.L[1],-0.5*P.L[2],+H),(0.5*P.L[1],0.5*P.L[2],+H); M = (20,20), dimorder = 3, qrule = Trap)
    low = HorizontalStraightPlane((-0.5*P.L[1],-0.5*P.L[2],-H),(0.5*P.L[1],0.5*P.L[2],-H); M = (20,20), dimorder = 3, qrule = Trap)
    xi  = [q.coords[1] for q in upp.dofs]
    yi  = [q.coords[2] for q in upp.dofs]

    σw = lmul!(wgfmatrix(G,w),σ)
    u1 = scatpotential(P,upp,G)*σw
    u2 = scatpotential(P,low,G)*σw

    if FRO
        u1 += scatcorrection(P,G,[q.coords for q in upp.dofs],σw; H=+H, δ = 0.75*P.pde[1].k)
        u2 += scatcorrection(P,G,[q.coords for q in low.dofs],σw; H=+H, δ = 0.75*P.pde[1].k)
    else
        @info "Non corrected energy balance test"
    end

    transmitted(P,H,u1,u2,xi,yi)
end

function transmitted(P::Problem{3,2},H::Float64,u1::Vector{ComplexF64},u2::Vector{ComplexF64},xi::Vector{Float64},yi::Vector{Float64})
    S = 0.0+im*0
    nC= 10
    α₁,α₂,β = P.pde[1].k*P.dir
    for n1 = -nC:nC, n2 = -nC:nC
        αₘ₁ = α₁ + n1*2π/P.L[1]
        αₘ₂ = α₂ + n2*2π/P.L[2]
        if P.pde[1].k^2 ≥ abs(αₘ₁)^2 + abs(αₘ₂)^2
            βₙ = sqrt((P.pde[1].k^2 - αₘ₁^2 - αₘ₂^2 ))
            # B⁺ = exp(-im*βₙ*H)/P.L[1]/P.L[2] * sum( u1.*exp.(-im*αₘ₁*xi-im*αₘ₂*yi) ) *P.L[1]*P.L[2]/length(u1)
            B⁻ = exp(-im*βₙ*H)/P.L[1]/P.L[2] * sum( u2.*exp.(-im*αₘ₁*xi-im*αₘ₂*yi) ) *P.L[1]*P.L[2]/length(u1)
            if n1 == 0 && n2 == 0
                global B₀ = B⁻
            end
            S += βₙ/β * abs(B⁻)^2 # + βₙ/β * abs(B⁺)^2+
        end
    end
    # abs(2*real(B₀)+S)
    1 + 2*real(B₀) + S 
end



obstaclesizes = collect(0.1:0.05:0.25) * L
pts = locallyperiodic(obstaclesizes)

import Plots
Plots.scatter(obstaclesizes/L, pts)