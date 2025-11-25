Here is a basic overview of the project layout.
```
- ./
  |- Makefile               # Build instructions
  |- Manifest.toml          # Package lock file
  |- Project.toml           # Project metadata and dependencies
  |- README.md              # Project overview documentation
  |- RepositoryLayout.md    # This file!
  |- SimulationConfig.toml  # Simulation parameters
  |- src/                   # Main project source files
  | |- CSVUtil.jl                 # Utility functions for working with CSV (and ZIP) files.
  | |- Graphs.jl                  # Defines a basic undirected graph structure and some common graph types.
  | |- HubbardDiagonalization.jl  # Entry point and main logic
  | |- StateEnumeration.jl        # Helper functions for enumerating states under various constraints
  | \- SymmetricMatrices.jl       # Defines a SymmetricMatrix type that stores only the lower-triangular part of the matrix.
  \- tests/                 # Test files
    \- grids                  # Test output against sample 2x2 grid results
      |- *.zip                  # Test datasets
      \- TestGrids.jl           # Test script
```
