module HubbardDiagonalization

# Include our submodules
# Base modules (no dependencies)
include("CSVUtil.jl")
include("Graphs.jl")
include("StateEnumeration.jl")
include("SymmetricMatrices.jl")

# Higher-level modules (depend on other modules) (order matters here)
include("ExactDiagonalization.jl")
include("DataHelpers.jl")

end
