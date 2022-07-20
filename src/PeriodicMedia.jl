module PeriodicMedia

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
Obstacle,
Window,
unitcell,
extendedcell,
solver,
energytest,
cellsolution,
viewsolution,
XYviewsolution,
YZviewsolution,
XZviewsolution

end