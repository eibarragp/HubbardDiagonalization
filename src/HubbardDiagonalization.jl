module HubbardDiagonalization

# Plots must be imported before CUDA to prevent a dependency conflict with the GR backend.
import Plots
import CUDA

# Include our submodules
# Base modules (no dependencies)
include("CSVUtil.jl")
include("Graphs.jl")
include("StateEnumeration.jl")
include("SymmetricMatrices.jl")
include("Utils.jl")

# Higher-level modules (depend on other modules) (order matters here)
include("ExactDiagonalization.jl")
include("DataHelpers.jl")

end
