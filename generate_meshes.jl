using FerriteGmsh

for sz ∈ [0.1, 0.05, 0.03, 0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001]
    gmsh.initialize()
    gmsh.model.add("heart2d")

    gmsh.model.occ.addPoint( 0.00, 0.000, 0.0, sz)
    gmsh.model.occ.addPoint(-1.00, 0.675, 0.0, sz)
    gmsh.model.occ.addPoint(-0.28, 1.300, 0.0, sz)
    gmsh.model.occ.addPoint( 0.00, 0.850, 0.0, sz)
    gmsh.model.occ.addPoint( 1.00, 0.675, 0.0, sz)
    gmsh.model.occ.addPoint( 0.28, 1.300, 0.0, sz)

    gmsh.model.occ.addBSpline([1, 2, 3, 4], 10)
    gmsh.model.occ.addBSpline([1, 5, 6, 4], 11)

    gmsh.model.occ.addCurveLoop([10, 11])
    gmsh.model.occ.addPlaneSurface([1])

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.generate(3)

    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()

    gmsh.write("data/meshes/heart2d-initial-$sz.msh")

    gmsh.finalize()
end
