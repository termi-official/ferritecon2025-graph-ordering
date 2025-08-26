using BenchmarkTools, Random

include("heat_assembly.jl")

normal_times = []
random_times = []

for N in 1:10
    normal_grid = generate_grid(Quadrilateral, (2^N, 2^N));
    @info N
    begin
        cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(normal_grid, RefQuadrilateral)
        normal_time = @benchmark assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
        push!(normal_times, normal_time)
    end
    rng = MersenneTwister(1337)
    random_grid = Grid(normal_grid.cells[randperm(rng, getncells(normal_grid))], normal_grid.nodes);
    begin
        cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(random_grid, RefQuadrilateral)
        random_time = @benchmark assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
        push!(random_times, random_time)
    end
end

using CairoMakie

# Number of elements vs time in ms
with_theme(theme_ggplot2()) do
    f = Figure(fontsize=22)
    ax = Axis(f[1,1], xlabel="Number of elements", ylabel="Assembly time in seconds", xscale=log2, yscale=log2)
    scatterlines!(ax, 4 .^collect(1:10), map(x->x.time, median.(normal_times)) ./ 1e9, label="Standard ordering")
    scatterlines!(ax, 4 .^collect(1:10), map(x->x.time, median.(random_times)) ./ 1e9, label="Random ordering")
    axislegend(ax, position=:lt)

    random_max = last(median.(random_times)).time / 1e9
    annotation!(ax, 2^11, random_max/2^5, 2^20, random_max,
        text = "$random_max s",
        path = Ann.Paths.Arc(0.3),
        style = Ann.Styles.LineArrow(),
        labelspace = :data
    )

    normal_max = last(median.(normal_times)).time / 1e9
    annotation!(ax, 2^15, normal_max/2^10, 2^20, normal_max,
        text = "$normal_max s",
        path = Ann.Paths.Arc(-0.3),
        style = Ann.Styles.LineArrow(),
        labelspace = :data
    )

    save("random-assembly-peformance.svg", f)
    f
end
