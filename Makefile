.PHONY: clean debug default format repl run setup test_grids

default: run

clean:
	@rm -rf output

debug:
	@mkdir -p output
	@julia --project=. -e 'include("src/HubbardDiagonalization.jl"); using .HubbardDiagonalization; main("--debug")' 2>&1 | tee output/debug.log

format:
	@julia --project=. -e 'import JuliaFormatter; JuliaFormatter.format("."; always_for_in = true)'

repl:
	@julia --project=.

run:
	@julia --project=. -e 'include("src/HubbardDiagonalization.jl"); using

setup:
	@julia --project=. -e 'import Pkg; Pkg.instantiate()'

test_grids:
	@julia --project=. -e 'include("tests/grids/TestGrids.jl"); using .TestGrids'
