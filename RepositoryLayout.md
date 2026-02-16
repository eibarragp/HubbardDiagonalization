Here is a basic overview of the project layout.
```
- ./
  |- Makefile               # Build instructions
  |- Manifest.toml          # Package lock file
  |- Project.toml           # Project metadata and dependencies
  |- README.md              # Project overview documentation
  |- RepositoryLayout.md    # This file!
  |- SimulationConfig.toml  # Simulation parameters
  |- slurm/                 # SLURM job scripts
  |  |- generate_slurm_script.py  # Helper script to generate SLURM job scripts for multiple clusters with varying requirements
  |  \- template.sh               # Template SLURM job script used by generate_slurm_script
  |- src/                   # Main project source files
  | |- CSVUtil.jl                 # Utility functions for working with CSV (and ZIP) files.
  | |- DataHelpers.jl             # Helper functions to handle data importing/exporting
  | |- ExactDiagonalization.jl    # Exact diagonalization logic and observable definitions
  | |- Graphs.jl                  # Defines a basic undirected graph structure and some common graph types.
  | |- HubbardDiagonalization.jl  # Main module definition (ensures that all files are loaded in the right order)
  | |- Main.jl                    # Main entry point; Handles command line parsing and high-level logic
  | |- StateEnumeration.jl        # Helper functions for enumerating states under various constraints
  | \- SymmetricMatrices.jl       # Defines a SymmetricMatrix type that stores only the lower-triangular part of the matrix.
  \- tests/                 # Test files
    \- grids                  # Test output against sample 2x2 grid results
      |- *.zip                  # Test datasets
      \- TestGrids.jl           # Test script
```
