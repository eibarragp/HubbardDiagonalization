module HubbardDiagonalization

# Include our submodules (order matters here)
# Base modules (no dependencies)
include("CSVUtil.jl")
include("Graphs.jl")
include("StateEnumeration.jl")
include("SymmetricMatrices.jl")

# Higher-level modules (depend on the base modules)
include("DataHelpers.jl")
include("ExactDiagonalization.jl")

end
