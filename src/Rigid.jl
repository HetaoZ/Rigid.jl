module Rigid

# -------------------------------------------------
using LinearAlgebra, Statistics
using MathKits, PointInPoly
# -------------------------------------------------
const TEMPERAL_DIFF_SCHEME = "forward" # "central"

include("base.jl")
export RigidStructure
include("io_pre.jl")
export create_rigid, set_fixed_u!, set_fixed_omega!
include("solver.jl")
export advance!
include("utils.jl")

###
end
