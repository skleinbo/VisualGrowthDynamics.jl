## Creates a mesh object for the Wigner-Seitz cell of an 
## FCC lattice.

using LinearAlgebra: normalize
using FileIO
using GeometryBasics, GrowthDynamics, Polyhedra

lattice = Lattices.FCCLattice(5)
# reference point and neighbors
ip0 = Lattices.midpoint(lattice)
p0 = Lattices.midpointcoord(lattice)
ipn = Lattices.neighbors(lattice, ip0)
pn = Lattices.coord(lattice, ipn)

# plane normals
hp_ns = normalize.(pn.-p0)

# 
poly = mapreduce(∩, hp_ns) do n
    HalfSpace(Array(n), 1)
end |> polyhedron

msh = GeometryBasics.mesh(Polyhedra.Mesh(poly))

save("fccvertex.obj", msh)