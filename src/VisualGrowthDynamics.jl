module VisualGrowthDynamics

import Base.Iterators: product
import StaticArrays: SVector

import GrowthDynamics
import GrowthDynamics: TumorConfigurations, Lattices
import .TumorConfigurations: TumorConfiguration
import .Lattices: HexagonalLattice, CubicLattice, FCCLattice, RealLattice, dimension, neighbors, midpoint, spacings

import Graphs: neighborhood

import Makie: plot!, default_theme, Attributes, to_value, SceneLike, lift, @lift, meshscatter!, arrows!, @recipe
using Makie
using GeometryBasics
using Colors


include("ColorFunctions.jl")


DEFAULT_ATTRIBUTES = Attributes(
    markersize=1,
    size=(800,800),
    colorfunc=nothing,
    colormap=:watermelon
)

@recipe(TumorPlot) do scene
    merge(DEFAULT_ATTRIBUTES,
        Attributes(
            plane = nothing,
            genotype = nothing,
            radius = 0.0,
            shell = nothing
        )
    )
end
function Makie.plot!(p::TumorPlot{Tuple{<:TumorConfiguration{A}}}) where A<:RealLattice{Int64}
    state = p[1]
    sz = lift(s->size(s.lattice), state)

    geno_filter = lift(state, p[:genotype]) do state, genotype
        if isnothing(genotype) || iszero(genotype)
            state.meta.genotypes
        else
            state.meta.genotypes[neighborhood(state.phylogeny, genotype, typemax(Int64); dir=:in)]
        end
    end
    @debug geno_filter[]

    indices = @lift product(map(i->1:i, $sz)...)

    _BP = lift(p[:plane], indices) do P, indices
        if dimension(state[].lattice)==3 && isa(P, Lattices.Plane)
            BitArray(Lattices.euclidean_dist(SVector{3}(I), P) <= 1/2 for I in indices)
        else
            fill(true, sz[])
        end
    end

    mp = midpoint(state[].lattice)
    a = spacings(state[].lattice)[1]
    _BR = lift(p[:radius], indices) do r, indices
        if !iszero(r)
            BitArray(r-a/2 <= Lattices.dist(state[].lattice, CartesianIndex(I), mp) < r+a/2 for I in indices)
        else
            fill(true, sz[])
        end
    end

    flt = @lift begin
        @debug "Updating filter"
        BitArray( x!=0 && x in $geno_filter for x in $state.lattice.data) .& $_BP .& $_BR
    end

    if all(x->x==false, flt[])
        @warn "No cells to draw."
    end

    positions = @lift [ Lattices.coord($state.lattice, x) for x in CartesianIndices($flt) if $flt[x] ]

    colorfunc = p[:colorfunc]
    color = haskey(p, :color) ? p[:color] : Observable(nothing)
    if typeof(colorfunc[]) <: Function
        color = lift(colorfunc, state, flt) do colorfunc, state, flt
            begin
                c = colorfunc(state)[flt]
                reshape(c, length(c))
            end
        end
    elseif isnothing(colorfunc[]) && isnothing(color[])
        color = lift(state, flt) do state, flt
            reshape(state.lattice.data[flt], :)
        end
    else
        color = lift(p[:color], flt) do color, flt
            reshape(color[flt], :)
        end
    end
    p[:color] = color

    _mesh = meshvertex(typeof(state[].lattice))
    
    # meshscatter!(p, positions, marker=_mesh, markersize=1,
    # color=RGBA(0.0,0.0,0.0,1.0),
    # shading=true, raw=false,strokewidth=5,strokecolor=:black
    # )
    meshscatter!(p, positions, marker=_mesh, markersize=1.0,
    color=color, colormap=p[:colormap],
    shading=true, raw=false,strokewidth=5,strokecolor=:black
    )

    p
end

## Cubic lattice

function default_theme(scene, ::TumorPlot{Tuple{TumorConfiguration{CubicLattice{Int64, A}}}}) where A
end

meshvertex(::Type{HexagonalLattice{A,B}}) where {A,B} = GeometryBasics.mesh(map([
    # Vertices of a hex with flat edge down
    Point2f0(1,0),
    Point2f0(cos(pi/3),-sin(pi/3)),
    Point2f0(-cos(pi/3),-sin(pi/3)),
    Point2f0(-1,0),
    Point2f0(-cos(pi/3),sin(pi/3)),
    Point2f0(cos(pi/3),sin(pi/3))
    ]/2/Float32(cos(pi/6))) do x
        # Rotation by π/6 --> pointy edge down
        [Point2f0(cos(pi/6),sin(pi/6)) Point2f0(-sin(pi/6),cos(pi/6))]*x
    end
    )

meshvertex(::Type{CubicLattice{A,B}}) where {A,B} = GeometryBasics.mesh(Rect(Vec3f0(0), Vec3f0(1)))
meshvertex(::Type{FCCLattice{A,B}}) where {A,B} = GeometryBasics.Sphere(Point3f(0), 1/4)




end # module
