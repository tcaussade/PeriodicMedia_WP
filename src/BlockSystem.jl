"""
    All functions related to integral block operator assembling and solving are defined here
"""

function solver(P::Problem{N,1}, G::Vector{NystromMesh{N,T,M,NM}}, Gt::Vector, w::Window; FRO = true) where {N,T,M,NM}
    MB = integralblockassembler(P,G)
    E  = diagonal(P,G)
    b  = rightside(P,G)
    W  = wgfmatrix(G,w)
    
    if FRO
        @assert Gt[1] isa Dict
        MB += finiterankoperator(P,G,Gt; δ = 0.75*P.pde[1].k, H = w.c*w.A)
    else
        @info "Solving without corrections"
    end

    return (E+MB*W)\b
end

function diagonal(P::Problem{N,1},G::Vector{NystromMesh{N,T,M,NM}}) where {N,T,M,NM}
    N1 = length(G[1].dofs)
    N2 = length(G[2].dofs)
    Diagonal([ones(N1); 0.5*(1+P.η)*ones(N1); P.γ*ones(2*N2)])
end

function rightside(P::Problem{N,1}, G::Vector{NystromMesh{N,T,M,NM}}) where {N,T,M,NM}
    N2 = length(G[2].dofs)
    ui = uinc(P,[q.coords for q in G[1].dofs])
    dui= ∂uinc(P,[q.coords for q in G[1].dofs],[WavePropBase.normal(q) for q in G[1].dofs])
    [ui;dui; zeros(2*N2)]
end

function integralblockassembler(P::Problem{N,1},G::Vector{NystromMesh{N,T,M,NM}}) where {N,T,M,NM}
    [transmissionequation(P,G); quasiperiodicequation(P,G)]
end

function transmissionequation(P::Problem{N,1},G::Vector{NystromMesh{N,T,M,NM}}) where {N,T,M,NM}
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

