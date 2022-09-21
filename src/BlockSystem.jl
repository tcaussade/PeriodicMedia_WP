"""
    All functions related to integral block operator assembling and solving are defined here
"""

""" solver Core function """

# function solver(P::Problem{N,1}, G::Vector{NystromMesh{N,T,M,NM}}, Gt::Vector, w::Window; FRO = true) where {N,T,M,NM}
  
function matrixcreator(P::Problem{N,1},G,Wparams::Window) where {N}

    MB = integralblockassembler(P,G.Γs,G.Γt)
    E  = diagonal(P,G.Γs)
    W  = wgfmatrix(G.Γs,Wparams)
    
    return MB,W,E 
end
function matrixcreator(P::Problem{3,2}, G, Wparams::Window)
    MB = integralblockassembler(P,G.Γs,G.Γt)
    E  = diagonal(P,G.Γs)
    # b  = rightside(P,G.Γs)
    W  = wgfmatrix(G.Γs,Wparams)

    return MB,W,E 

    # if FRO
    #     @assert Gt[1] isa Dict{Tuple{Int,Int}}
    #     MB += finiterankoperator(P,G,Gt; δ = 0.75*P.pde[1].k, H = w.c*w.A)
    # else
    #     @info "Solving without corrections"
    # end

    # return gmres(E+MB*W,b; restart = size(MB,2), verbose = true, reltol = 1e-8)
end

""" diagonal accounts the identity term that ensures second-kind Fredholmness """

function diagonal(P::Problem{N,1},G::Vector{NystromMesh{N,T,M,NM}}) where {N,T,M,NM}
    N1 = length(G[1].dofs)
    N2 = length(G[2].dofs)
    Diagonal([ones(N1); 0.5*(1+P.η)*ones(N1); P.γ*ones(2*N2)])
end
function diagonal(P::Problem{3,2},G::Vector{NystromMesh{3,T,M,NM}}) where {T,M,NM}
    N1 = length(G[1].dofs)
    N2 = length(G[2].dofs)
    N3 = length(G[4].dofs)
    Diagonal([ones(N1); 0.5*(1+P.η)*ones(N1); P.γ[1]^3*ones(2*N2); P.γ[2]^3*ones(2*N3)])
end

""" rightside  accounts the incident planewave """

function rightside(P::Problem{N,1}, G::Vector{NystromMesh{N,T,M,NM}}) where {N,T,M,NM}
    N2 = length(G[2].dofs)
    ui = uinc(P,[q.coords for q in G[1].dofs])
    dui= ∂uinc(P,[q.coords for q in G[1].dofs],[WavePropBase.normal(q) for q in G[1].dofs])
    [ui;dui; zeros(2*N2)]
end
function rightside(P::Problem{3,2}, G::Vector{NystromMesh{3,T,M,NM}}) where {T,M,NM}
    N2 = length(G[2].dofs)
    N3 = length(G[4].dofs)
    ui = uinc(P,[q.coords for q in G[1].dofs])
    dui= ∂uinc(P,[q.coords for q in G[1].dofs],[WavePropBase.normal(q) for q in G[1].dofs])
    [ui;dui; zeros(2*N2); zeros(2*N3)]
end

""" integralblockassembler discretizes the integral operators """ 

function integralblockassembler(P::Problem{N,1},G::Vector{NystromMesh{N,T,M,NM}},Gt::Vector) where {N,T,M,NM}
    [transmissionequation(P,G); quasiperiodicequation(P,G)]
end
function integralblockassembler(P::Problem{3,2},G::Vector{NystromMesh{3,T,M,NM}},Gt::Vector) where {T,M,NM}
    [transmissionequation(P,G); xquasiperiodicequation(P,Gt); yquasiperiodicequation(P,Gt)]
end

""" transmissionequation computes the two equations related to imposing transmission conditions for dirichler and neumann trace
    methods: 
        - join blocks
        - discretize for 2D1D periodic arrays
        - discretize for 3D2D periodic arrays 
"""

function transmissionequation(P::Problem{N,NP},G::Vector{NystromMesh{N,T,M,NM}}) where {N,T,M,NM,NP}
    eq1 = transmissionequation(P,G, DoubleLayerOperator, SingleLayerOperator)
    eq2 = transmissionequation(P,G, HyperSingularOperator, AdjointDoubleLayerOperator)
    [eq1;eq2]
end
function transmissionequation(P::Problem{N,1}, G::Vector{NystromMesh{N,T,M,NM}}, opD, opN) where {N,T,M,NM}
    D1,N1 = Matrix( Nystrom.assemble_dim(opD(P.pde[1], G[1]))|> Nystrom.materialize ), Matrix( Nystrom.assemble_dim(opN(P.pde[1], G[1]))|> Nystrom.materialize )
    D2,N2 = Matrix( Nystrom.assemble_dim(opD(P.pde[2], G[1]))|> Nystrom.materialize ), Matrix( Nystrom.assemble_dim(opN(P.pde[2], G[1]))|> Nystrom.materialize )
    Mq1 = axpy!(-1.,D1,D2)
    Mq2 = axpby!(P.η,N1,-1.,N2)

    D12,N12 = Matrix( opD(P.pde[1],G[1],G[2]) ), Matrix( opN(P.pde[1],G[1],G[2]) )
    D13,N13 = Matrix( opD(P.pde[1],G[1],G[3]) ), Matrix( opN(P.pde[1],G[1],G[3]) )
    Mq3 = axpby!(P.γ,D13,-1.,D12)
    Mq4 = axpy!(-P.γ, N13,N12)

    [Mq1 Mq2 Mq3 Mq4]
end
function transmissionequation(P::Problem{3,2}, G::Vector{NystromMesh{3,T,M,NM}}, opD, opN) where {T,M,NM}
    D1,N1 = Matrix( Nystrom.assemble_dim(opD(P.pde[1], G[1]))|> Nystrom.materialize ), Matrix( Nystrom.assemble_dim(opN(P.pde[1], G[1]))|> Nystrom.materialize )
    D2,N2 = Matrix( Nystrom.assemble_dim(opD(P.pde[2], G[1]))|> Nystrom.materialize ), Matrix( Nystrom.assemble_dim(opN(P.pde[2], G[1]))|> Nystrom.materialize )
    Mp1 = axpy!(-1.,D1,D2)
    Mp2 = axpby!(P.η,N1,-1.,N2)

    D12,N12 = Matrix( opD(P.pde[1],G[1],G[2]) ), Matrix( opN(P.pde[1],G[1],G[2]) )
    D13,N13 = Matrix( opD(P.pde[1],G[1],G[3]) ), Matrix( opN(P.pde[1],G[1],G[3]) )
    Mp3 = axpby!(P.γ[1],D13,-1.,D12)
    Mp4 = axpy!(-P.γ[1],N13,N12)

    D14,N14 = Matrix( opD(P.pde[1],G[1],G[4]) ), Matrix( opN(P.pde[1],G[1],G[4]) )
    D15,N15 = Matrix( opD(P.pde[1],G[1],G[5]) ), Matrix( opN(P.pde[1],G[1],G[5]) )
    Mp5 = axpby!(P.γ[2],D15,-1.,D14)
    Mp6 = axpy!(-P.γ[2],N15,N14)

    [Mp1 Mp2 Mp3 Mp4 Mp5 Mp6]
end

""" imposes quasiperiodicity for 2D1D arrays """

function quasiperiodicequation(P::Problem{N,1}, G::Vector{NystromMesh{N,T,M,NM}}) where {N,T,M,NM}
    eq3 = quasiperiodicequation(P,G, DoubleLayerOperator, SingleLayerOperator)
    eq4 = quasiperiodicequation(P,G, HyperSingularOperator, AdjointDoubleLayerOperator)
    [eq3;eq4]
end
function quasiperiodicequation(P::Problem{N,1}, G::Vector{NystromMesh{N,T,M,NM}}, opD, opN) where {N,T,M,NM}
    D21,N21 = Matrix( opD(P.pde[1],G[2],G[1]) ), Matrix( opN(P.pde[1],G[2],G[1]) )
    D31,N31 = Matrix( opD(P.pde[1],G[3],G[1]) ), Matrix( opN(P.pde[1],G[3],G[1]) )
    Mq1 = axpby!(-P.γ,D21,-1,D31)
    Mq2 = axpby!(P.η*P.γ,N21,P.η,N31)

    D23,N23 = Matrix( opD(P.pde[1],G[2],G[3]) ), Matrix( opN(P.pde[1],G[2],G[3]) )
    D32,N32 = Matrix( opD(P.pde[1],G[3],G[2]) ), Matrix( opN(P.pde[1],G[3],G[2]) )
    Mq3 = axpby!(P.γ^2,D23,-1.,D32)
    Mq4 = axpy!(-P.γ^2,N23,N32)

    [Mq1 Mq2 Mq3 Mq4]
end

""" imposes quasiperiodicity for 3D2D arrays, i.e. in x-axis and y-axis """

function xquasiperiodicequation(P::Problem{3,2}, G::Vector) 
    eq3 = xquasiperiodicequation(P,G, DoubleLayerOperator, SingleLayerOperator)
    eq4 = xquasiperiodicequation(P,G, HyperSingularOperator, AdjointDoubleLayerOperator)
    [eq3;eq4]
end
function xquasiperiodicequation(P::Problem{3,2}, G::Vector, opD, opN)
    Γ₂,Γ₃ = [get(G[i],0,"") for i=2:3]
    # D21,N21 = superoperator(P.pde[1],Γ₂,G[1],P.γ[1], opD; axis = "x"), superoperator(P.pde[1],Γ₂,G[1],P.γ[1], opN; axis = "x")
    # D31,N31 = superoperator(P.pde[1],Γ₃,G[1],P.γ[1], opD; axis = "x"), superoperator(P.pde[1],Γ₃,G[1],P.γ[1], opN; axis = "x")
    D21,N21 = superoperator(P.pde[1],Γ₂,G[1],P.γ, opD), superoperator(P.pde[1],Γ₂,G[1],P.γ, opN)
    D31,N31 = superoperator(P.pde[1],Γ₃,G[1],P.γ, opD), superoperator(P.pde[1],Γ₃,G[1],P.γ, opN)
    Mp1 = lmul!(-P.γ[1],axpy!(P.γ[1]^3,D21,D31))
    Mp2 = lmul!(P.η*P.γ[1],axpy!(P.γ[1]^3,N21,N31))

    # D23,N23 = Matrix(opD(P.pde[1],Γ₂,Γ₃)), Matrix(opN(P.pde[1],Γ₂,Γ₃))
    # D32,N32 = Matrix(opD(P.pde[1],Γ₃,Γ₂)), Matrix(opN(P.pde[1],Γ₃,Γ₂))
    D23,N23 = superoperator(P.pde[1],Γ₂,G[3],P.γ[2],opD), superoperator(P.pde[1],Γ₂,G[3],P.γ[2],opN)
    D32,N32 = superoperator(P.pde[1],Γ₃,G[2],P.γ[2],opD), superoperator(P.pde[1],Γ₃,G[2],P.γ[2],opN)
    Mp3 = axpby!(P.γ[1]^6,D23,-1.,D32)
    Mp4 = axpy!(-P.γ[1]^6,N23,N32)

    D24,N24 = superoperator(P.pde[1],Γ₂,G[4],P.γ[1],opD), superoperator(P.pde[1],Γ₂,G[4],P.γ[1],opN)
    D34,N34 = superoperator(P.pde[1],Γ₃,G[4],P.γ[1],opD), superoperator(P.pde[1],Γ₃,G[4],P.γ[1],opN)
    D25,N25 = superoperator(P.pde[1],Γ₂,G[5],P.γ[1],opD), superoperator(P.pde[1],Γ₂,G[5],P.γ[1],opN)
    D35,N35 = superoperator(P.pde[1],Γ₃,G[5],P.γ[1],opD), superoperator(P.pde[1],Γ₃,G[5],P.γ[1],opN)
    td1 = axpy!(P.γ[1]^3,D24,D34)
    td2 = axpy!(P.γ[1]^3,D25,D35)
    # Mp5 = lmul!(-P.γ[1]/P.γ[2],axpy!(-P.γ[2]^3,td2,td1))
    Mp5 = axpby!(-P.γ[1]/P.γ[2],td1, +P.γ[1]*P.γ[2]^2,td2)
    tn1 = axpy!(P.γ[1]^3,N24,N34)
    tn2 = axpy!(P.γ[1]^3,N25,N35)
    # Mp6 = lmul!(+P.γ[1]/P.γ[2],axpy!(-P.γ[2]^3,tn2,tn1))
    Mp6 = axpby!(+P.γ[1]/P.γ[2],tn1, -P.γ[1]*P.γ[2]^2,tn2)

    [Mp1 Mp2 Mp3 Mp4 Mp5 Mp6]
end

function yquasiperiodicequation(P::Problem{3,2}, G::Vector) 
    eq5 = yquasiperiodicequation(P,G, DoubleLayerOperator, SingleLayerOperator)
    eq6 = yquasiperiodicequation(P,G, HyperSingularOperator, AdjointDoubleLayerOperator)
    [eq5;eq6]
end
function yquasiperiodicequation(P::Problem{3,2},G::Vector,opD,opN)
    Γ₄,Γ₅ = [get(G[i],0,"") for i=4:5]
    # D41,N41 = superoperator(P.pde[1],Γ₄,G[1],P.γ[1], opD; axis = "y"), superoperator(P.pde[1],Γ₄,G[1],P.γ[1], opN; axis = "y")
    # D51,N51 = superoperator(P.pde[1],Γ₅,G[1],P.γ[1], opD; axis = "y"), superoperator(P.pde[1],Γ₅,G[1],P.γ[1], opN; axis = "y")
    D41,N41 = superoperator(P.pde[1],Γ₄,G[1],P.γ, opD), superoperator(P.pde[1],Γ₄,G[1],P.γ, opN)
    D51,N51 = superoperator(P.pde[1],Γ₅,G[1],P.γ, opD), superoperator(P.pde[1],Γ₅,G[1],P.γ, opN)
    Mp1 = lmul!(-P.γ[2],axpy!(P.γ[2]^3,D41,D51))
    Mp2 = lmul!(P.η*P.γ[2],axpy!(P.γ[2]^3,N41,N51))

    D42,N42 = superoperator(P.pde[1],Γ₄,G[2],P.γ[2],opD), superoperator(P.pde[1],Γ₄,G[2],P.γ[2],opN)
    D52,N52 = superoperator(P.pde[1],Γ₅,G[2],P.γ[2],opD), superoperator(P.pde[1],Γ₅,G[2],P.γ[2],opN)
    D43,N43 = superoperator(P.pde[1],Γ₄,G[3],P.γ[2],opD), superoperator(P.pde[1],Γ₄,G[3],P.γ[2],opN)
    D53,N53 = superoperator(P.pde[1],Γ₅,G[3],P.γ[2],opD), superoperator(P.pde[1],Γ₅,G[3],P.γ[2],opN)
    td1 = axpy!(P.γ[2]^3,D42,D52)
    td2 = axpy!(P.γ[2]^3,D43,D53)
    # Mp3 = lmul!(-P.γ[2]/P.γ[1],axpy!(-P.γ[1]^3,td2,td1))
    Mp3 = axpby!(-P.γ[2]/P.γ[1],td1, +P.γ[2]*P.γ[1]^2,td2)
    tn1 = axpy!(P.γ[2]^3,N42,N52)
    tn2 = axpy!(P.γ[2]^3,N43,N53)
    # Mp4 = lmul!(-P.γ[2]/P.γ[1],axpy!(-P.γ[1]^3,tn2,tn1))
    Mp4 = axpby!(+P.γ[2]/P.γ[1],tn1, -P.γ[2]*P.γ[1]^2,tn2)

    # D45,N45 = Matrix(opD(P.pde[1],Γ₄,Γ₅)), Matrix(opN(P.pde[1],Γ₄,Γ₅))
    # D54,N54 = Matrix(opD(P.pde[1],Γ₅,Γ₄)), Matrix(opN(P.pde[1],Γ₅,Γ₄))
    D45,N45 = superoperator(P.pde[1],Γ₄,G[5],P.γ[1],opD), superoperator(P.pde[1],Γ₄,G[5],P.γ[1],opN)
    D54,N54 = superoperator(P.pde[1],Γ₅,G[4],P.γ[1],opD), superoperator(P.pde[1],Γ₅,G[4],P.γ[1],opN)
    Mp5 = axpby!(P.γ[2]^6,D45,-1.,D54)
    Mp6 = axpy!(-P.γ[2]^6,N45,N54)

    [Mp1 Mp2 Mp3 Mp4 Mp5 Mp6]
end


""" superoperators are used in 3D2D method to avoid corner singularities 
    if source isa Dict{Int}           : displace the density over the vertical boundary accordingly
    if source isa Dict{Tuple{Int,Int}}: displace densities over obstacles
"""

function superoperator(pde,target::NystromMesh,source::Dict{Int}, γ::ComplexF64,op)
    # sop = zeros(ComplexF64,length(target.dofs),length(get(source,0,"").dofs))
    # for k = -1:1
    #     sop += lmul!(γ^k, Matrix(op(pde,target,get(source,k,""))))
    # end
    # return sop
    sum( [lmul!(γ^k, Matrix(op(pde,target,get(source,k,"")))) for k=-1:1] )
end
function superoperator(pde,target::NystromMesh,source::Dict{Tuple{Int,Int}},γ::Vector{ComplexF64},op)
    # sop = zeros(ComplexF64,length(target.dofs),length(get(source,(0,0),"").dofs))
    # if axis == "x"
    #     for k = -1:1
    #         sop += lmul!(γ^k, Matrix( op(pde,target,get(source,(k,0),""))) )
    #     end
    # elseif axis == "y"
    #     for k = -1:1
    #         sop += lmul!(γ^k, Matrix( op(pde,target,get(source,(0,k),""))) )
    #     end
    # end
    # for i = -1:1, j = -1:1
    #     sop += lmul!(γ[1]^i*γ[2]^j, Matrix( op(pde,target,get(source,(i,j),""))) )
    # end
    # return sop
    sum( [lmul!(γ[1]^i*γ[2]^j, Matrix( op(pde,target,get(source,(i,j),""))) ) for i=-1:1, j=-1:1] )
end

