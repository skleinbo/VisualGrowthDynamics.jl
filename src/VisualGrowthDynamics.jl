module VisualGrowthDynamics

export TumorInspector, tumorplot

import Base.Iterators: product
using Colors
using GeometryBasics
import GLMakie
import GrowthDynamics
import GrowthDynamics: TumorConfigurations, Lattices
import GrowthDynamics.Lattices: midpoint, midpointcoord
import Graphs: neighborhood
import .Lattices: HexagonalLattice, CubicLattice, FCCLattice, RealLattice, Plane
import .Lattices: coord, dimension, neighbors, midpoint, spacings
import LinearAlgebra: dot, normalize
import Makie: plot!, default_theme, Attributes, to_value, SceneLike, lift, @lift, meshscatter!, arrows!, @recipe
using Makie
import StaticArrays: SVector
import .TumorConfigurations: TumorConfiguration

include("ColorFunctions.jl")
include("utility.jl")
include("recipes.jl")
include("inspector.jl")
import .ColorFunctions: color_depth, color_lineages

end # module
