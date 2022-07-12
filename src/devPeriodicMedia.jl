module devPeriodicMedia

using LinearAlgebra
using StaticArrays
using WavePropBase
using Nystrom
import Plots
using IterativeSolvers

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