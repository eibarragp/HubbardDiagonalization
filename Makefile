.PHONY: clean debug default format repl run setup test_grids

default: run

JULIA=@julia --project=.

clean:
	@rm -rf output

debug:
	@mkdir -p output
	$(JULIA) src/Main.jl --debug 2>&1 | tee output/debug.log

format:
	$(JULIA) -e 'import JuliaFormatter; JuliaFormatter.format("."; always_for_in = true)'

repl:
	$(JULIA)

run:
	$(JULIA) src/Main.jl

setup:
	$(JULIA) -e 'import Pkg; Pkg.instantiate()'

test_grids:
	$(JULIA) tests/grids/TestGrids.jl
