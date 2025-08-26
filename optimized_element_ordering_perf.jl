using BenchmarkTools, JLD2, Metis

include("heat_assembly.jl")
include("optimization_utils.jl")

meshpath = joinpath("data", "meshes")

function reverse_triangles!(grid)
    for cellid in 1:getncells(grid)
        nodes = grid.cells[cellid].nodes
        grid.cells[cellid] = Triangle((nodes[3], nodes[2], nodes[1]))
    end
end

# gen_times = []
# gen_num_elements = []

opt1_times = []
opt2_times = []
opt3_times = []
opt_num_elements = []


for fname in readdir(meshpath)
    # 
    endswith(fname, ".jld2") || continue

    @info fname
    if startswith(fname, "heart2d-initial-")
        grid = load(joinpath(meshpath, fname), "grid")
        reverse_triangles!(grid)
        cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(grid, RefTriangle)
        gen_time = @benchmark assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
        @info gen_time
        push!(gen_times, gen_time)
        push!(gen_num_elements, getncells(grid))
    elseif startswith(fname, "heart2d-optimized-")
        grid = load(joinpath(meshpath, fname), "grid")
        reverse_triangles!(grid)
        push!(opt_num_elements, getncells(grid))

        begin
            cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(grid, RefTriangle)
            opt_time = @benchmark assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
            @info opt_time
            push!(opt1_times, opt_time)
        end

        begin
            optimize_nodes!(grid)
            cellvalues, assembler, Ke, fe, dh = setup_heat_assembly(grid, RefTriangle)
            opt_time = @benchmark assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
            @info opt_time
            push!(opt2_times, opt_time)

            renumber!(dh, DofOrder.Ext{Metis}())
            ## Create a new compatible assembler
            assembler = start_assemble(allocate_matrix(dh), zeros(ndofs(dh)))
            opt_time = @benchmark assemble_heat_global($cellvalues, $assembler, $Ke, $fe, $dh)
            @info opt_time
            push!(opt3_times, opt_time)
        end
    end
end


gen_order = sortperm(gen_num_elements)
opt_order = sortperm(opt_num_elements)

using CairoMakie

# Number of elements vs time in ms
with_theme(theme_ggplot2()) do
    f = Figure(fontsize=22)
    ax = Axis(f[1,1], xlabel="Number of elements", ylabel="Assembly time in seconds", xscale=log2, yscale=log2)
    scatterlines!(gen_num_elements[gen_order], map(x->x.time, median.(gen_times))[gen_order] ./ 1e9, label="Standard ordering")
    scatterlines!(opt_num_elements[opt_order], map(x->x.time, median.(opt1_times))[opt_order] ./ 1e9, label="Optimized E ordering")
    scatterlines!(opt_num_elements[opt_order], map(x->x.time, median.(opt2_times))[opt_order] ./ 1e9, label="Optimized E+N ordering")
    scatterlines!(opt_num_elements[opt_order], map(x->x.time, median.(opt3_times))[opt_order] ./ 1e9, label="Optimized E+N+D ordering")
    # scatterlines!(optr_num_elements[opt_order], map(x->x.time, median.(optr_times))[opt_order] ./ 1e9, label="Rev. optimized element ordering")
    axislegend(ax, position=:lt)

    save("optimized-assembly-peformance.svg", f)
    f
end
