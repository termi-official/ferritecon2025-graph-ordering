using Ferrite, SparseArrays
# using MatrixBandwidth
using CuthillMcKee

function reorder_nodes_on_cell(cell::CellType, ordering) where CellType <: Ferrite.AbstractCell
    CellType(map(n->ordering[n], cell.nodes))
    # CellType(map(n->findfirst(x->x==n, ordering), cell.nodes))
end

function reorder_nodes!(g, ordering)
    oldnodes = copy(g.nodes)
    for i in eachindex(oldnodes)
        g.nodes[ordering[i]] = oldnodes[i]
        # g.nodes[i] = oldnodes[ordering[i]]
    end
    g.cells = [reorder_nodes_on_cell(cell, ordering) for cell in g.cells]
    return g
end

function optimize_nodes!(g)
    I, J = Int[], Int[]
    for cellid in 1:getncells(g)
        cell = getcells(g, cellid)
        for v1 in Ferrite.get_node_ids(cell)
            for v2 in Ferrite.get_node_ids(cell)
                push!(I, v1)
                push!(J, v2)
            end
        end
    end

    incidence_matrix = SparseArrays.spzeros!(Bool, I, J, getnnodes(g), getnnodes(g))
    fill!(incidence_matrix.nzval, true)
    @info "Compute ordering..."
    # res = minimize_bandwidth(incidence_matrix, Minimization.CuthillMcKee())
    # @show res
    # ordering = res.ordering
    ordering = symrcm(incidence_matrix, false, false)
    # ordering = [getnnodes(g)-o+1 for o in ordering]

    ordering2 = copy(ordering)
    for i in eachindex(ordering)
        ordering2[ordering[i]] = i
    end
    @info "Reorder..."
    reorder_nodes!(g, ordering2)
end

grid = generate_grid(Quadrilateral, (20, 20))

o = shuffle(Vector(1:getnnodes(grid)))
reorder_nodes!(grid, o)

optimize_nodes!(grid)
optimize_nodes!(grid)

VTKGridFile("FerriteCon-Debug", grid) do vtk
end
