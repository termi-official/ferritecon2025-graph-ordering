using JLD2, Gecko, Ferrite, SparseArrays


struct ProgressCallbacks <: Gecko.AbstractGeckoLogger
end

function Gecko.begin_order(data::ProgressCallbacks, graph::GeckoGraph, cost::Real)
    @info "Begin Order ($cost)"
    return nothing
end
function Gecko.end_order(data::ProgressCallbacks, graph::GeckoGraph, cost::Real)
    @info "End Order ($cost)"
    return nothing
end
function Gecko.begin_iter(data::ProgressCallbacks, graph::GeckoGraph, iter::Unsigned, maxiter::Unsigned, window::Unsigned)
    @info "Begin Iter ($iter/$maxiter | $window)"
    return nothing
end
function Gecko.end_iter(data::ProgressCallbacks, graph::GeckoGraph, mincost::Real, cost::Real)
    @info "End Iter ($mincost | $cost)"
    return nothing
end
function Gecko.begin_phase(data::ProgressCallbacks, graph::GeckoGraph, name::Gecko.CxxWrap.StdLib.StdString)
    # @info "begin_phase $name"
    return nothing
end
function Gecko.end_phase(data::ProgressCallbacks, graph::GeckoGraph, show::Bool)
    # @info "end_phase $show"
    return nothing
end
function Gecko.quit(data::ProgressCallbacks)
    return false
end


meshpath = joinpath("data", "meshes")

for fname in filter(fname->startswith(fname, "heart2d-initial-") && endswith(fname, ".jld2"), readdir(meshpath))
    @info joinpath(meshpath, fname)
    sz = parse(Float64, replace(replace(fname, "heart2d-initial-" => ""), ".jld2" => ""))
    @info sz

    grid = load(joinpath(meshpath, fname), "grid")
    neighbormatrix = Ferrite.create_incidence_matrix(grid)

    gg = GeckoGraph(getncells(grid))
    # Add connectivity
    for cellid in 1:getncells(grid)
        for ni in nzrange(neighbormatrix, cellid)
            add_directed_edge!(gg, cellid, neighbormatrix.rowval[ni])
        end
    end

    # Optimize
    logger     = ProgressCallbacks()

    @info "Ordering..."
    order!(gg, GraphOrderingParameters(;iterations=4, window=4), logger)

    # Store
    cells = [grid.cells[new_index(gg, cellid)] for cellid in 1:getncells(grid)]
    newgrid = Grid(cells, grid.nodes)
    jldsave(joinpath(meshpath, "heart2d-optimized-$sz.jld2"); grid=newgrid)
end
