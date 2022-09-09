function rw_poly(n; d,L)
    a = 1-d[1]^2-d[2]^2
    b = -4π*(n[1]*d[1]/L[1]+n[2]*d[2]/L[2])
    c = - (2π)^2 * (n[1]^2/L[1]^2+n[2]^2/L[2]^2)
    Δ = b^2 - 4*a*c
    @show x1 = (-b - sqrt(Δ))/(2a)
    @show x2 = (-b + sqrt(Δ))/(2a)
    return
end

function βₙ(k,n,d,L)
    a = 1-d[1]^2-d[2]^2
    b = -4π*(n[1]*d[1]/L[1]+n[2]*d[2]/L[2])
    c = - (2π)^2 * (n[1]^2/L[1]^2+n[2]^2/L[2]^2)
    sqrt(complex(a*k^2+b*k+c))
end

global θ = [π/4, π/4]
global L = [0.5,0.5]
global d = [sin(θ[1])*cos(θ[2]), sin(θ[1])*sin(θ[2])]

# Shows first RW anomaly
rw_poly([0,-1]; d=d,L=L)

n = collect(-3:3)
k = 9.199221756451442
k = 10.0
δ = 5/3*k
Cδ = Vector{Tuple{Int,Int}}(undef,0)
for i in n, j in n
    β = βₙ(k,[n[i+4],n[j+4]],d,L)
    abs(β) < δ ? push!(Cδ,(i,j)) : nothing
    abs(β) < 1e-14  ? println("RW: (i,j)=",(i,j)) : 
    real(β) > 0.0  ? println("Propagative: (i,j)=",(i,j)) : continue
end
@show length(Cδ)


println("\nVerification")
α₁,α₂ = k*d
for n1 = -3:3, n2 = -3:3
    αₘ₁ = α₁ + n1*2π/L[1]
    αₘ₂ = α₂ + n2*2π/L[2]
    βₘ = sqrt(complex(k^2 - αₘ₁^2 - αₘ₂^2 ))
    abs(βₘ) < 1e-6  ? println("RW: (i,j)=",(n1,n2)) : 
    real(βₘ) > 0.0  ? println("Propagative: (i,j)=",(n1,n2)) : continue
end