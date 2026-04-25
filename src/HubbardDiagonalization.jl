module HubbardDiagonalization

# Use MKL because it plays more nicely with Julia's multi-threading
# This must be the first package loaded!
using MKL

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
