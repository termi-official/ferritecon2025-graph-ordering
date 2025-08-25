using BenchmarkTools, Random

include("heat_assembly.jl")

normal_times = []
random_times = []

for N in 1:10
    normal_grid = generate_grid(Quadrilateral, (2^N, 2^N));
    @info N
    begin
        cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(normal_grid)
        normal_time = @benchmark assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
        push!(normal_times, normal_time)
    end
    rng = MersenneTwister(1337)
    random_grid = Grid(normal_grid.cells[randperm(rng, getncells(normal_grid))], normal_grid.nodes);
    begin
        cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(random_grid)
        random_time = @benchmark assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
        push!(random_times, random_time)
    end
end

using GLMakie

# Number of elements vs time in ms
lines(4 .^collect(1:10), map(x->x.time, median.(normal_times)) ./ 1e6)
lines!(4 .^collect(1:10), map(x->x.time, median.(random_times))  ./ 1e6)
