using Ferrite, SparseArrays, MatrixBandwidth

function reorder_nodes_on_cell(cell::CellType, ordering) where CellType <: Ferrite.AbstractCell
    # CellType(map(n->ordering[n], cell.nodes))
    CellType(map(n->findfirst(x->x==n, ordering), cell.nodes))
end

function reorder_nodes!(g, ordering)
    oldnodes = copy(g.nodes)
    for i in eachindex(oldnodes)
        # g.nodes[ordering[i]] = oldnodes[i]
        g.nodes[i] = oldnodes[ordering[i]]
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
    res = minimize_bandwidth(incidence_matrix, Minimization.CuthillMcKee())
    @show res
    ordering = res.ordering

    reorder_nodes!(g, ordering)
end

# grid = generate_grid(Quadrilateral, (20, 20))

# o = shuffle(Vector(1:getnnodes(grid)))
# reorder_nodes!(grid, o)

# optimize_nodes!(grid)
# optimize_nodes!(grid)
