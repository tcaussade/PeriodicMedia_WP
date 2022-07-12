module devPeriodicMedia

using LinearAlgebra
using WavePropBase
using Nystrom
using Plots

include("Base.jl")
include("Configuration.jl")
include("BlockSystem.jl")
include("Evaluation.jl")
include("EnergyBalance.jl")
include("FiniteRankOperator.jl")

export 
Problem,
Window,
unitcell,
solver,
scatpotential

end