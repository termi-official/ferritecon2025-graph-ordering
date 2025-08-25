using BenchmarkTools, Random

include("heat_assembly.jl")

for N in 1:10
    normal_grid = generate_grid(Quadrilateral, (2^N, 2^N));
    @info N
    begin
        cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(normal_grid)
        @btime assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
    end
    rng = MersenneTwister(1337)
    random_grid = Grid(normal_grid.cells[randperm(rng, getncells(normal_grid))], normal_grid.nodes);
    begin
        cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(random_grid)
        @btime assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
    end
end
