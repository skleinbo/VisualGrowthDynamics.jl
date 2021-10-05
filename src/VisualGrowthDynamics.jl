module VisualGrowthDynamics

import Base.Iterators: product
import StaticArrays: SVector

import GrowthDynamics
import GrowthDynamics: TumorConfigurations, Lattices
import .TumorConfigurations: TumorConfiguration
import .Lattices: HexagonalLattice, CubicLattice, neighbors, midpoint

import LightGraphs: neighborhood

import Makie: plot!, default_theme, Attributes, to_value, SceneLike, lift, @lift, meshscatter!, arrows!, @recipe
using Makie
using GeometryBasics
using Colors


include("ColorFunctions.jl")


DEFAULT_ATTRIBUTES = Attributes(
    markersize=1,
    size=(800,800),
    color=ColorFunctions.color_phylo
)

@recipe(TumorPlot) do scene
    merge(DEFAULT_ATTRIBUTES,
        Attributes(
            plane = nothing,
            genotype = nothing,
            radius = 0.0
        )
    )
end

function Makie.plot!(p::TumorPlot{Tuple{<:TumorConfiguration{HexagonalLattice{Int64}}}})
    state = to_value(p[1])
    L = state.lattice.Na

    hex_vertices = map([
            # Vertices of a hex with flat edge down
            Point2f0(1,0),
            Point2f0(cos(pi/3),-sin(pi/3)),
            Point2f0(-cos(pi/3),-sin(pi/3)),
            Point2f0(-1,0),
            Point2f0(-cos(pi/3),sin(pi/3)),
            Point2f0(cos(pi/3),sin(pi/3))
            ]/L/Float32(cos(pi/6))) do x
                # Rotation by π/6 --> pointy edge down
                [Point2f0(cos(pi/6),sin(pi/6)) Point2f0(-sin(pi/6),cos(pi/6))]*x
            end

    hex_mesh=GeometryBasics.mesh(hex_vertices)
    off = (sqrt(5)-sqrt(3))/2

    hex_positions = reshape(Point2f0[Point2f0(
        2*(q+ifelse(r%2==0,0.,1/2))/L-1.0,
        2*(r+1/2+off/2*ifelse(r==0,0.,-r))/L-1.0)
         for (r,q) in Iterators.product(0:L-1,0:L-1)], L^2)

    meshscatter!(p, hex_positions, marker=hex_mesh, markersize=1,
         color=lift(to_value(p[:color]), Node(state)),
        shading=false, raw=false)

    p
end

## Cubic lattice

function default_theme(scene, ::TumorPlot{Tuple{<:TumorConfiguration{CubicLattice{Int64, A}}}}) where A
    @info "Test"

end


function plot!(p::TumorPlot{Tuple{<:TumorConfiguration{CubicLattice{Int64, A}}}}) where A
    state = p[1]
    sz = lift(s->size(s.lattice), state)

    # if !haskey(p.attributes, :genotype)
    #     p.attributes[:genotype] = nothing
    # end
    # if !haskey(p.attributes, :radius)
    #     p.attributes[:radius] = 0.0
    # end

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
        if isa(P, Lattices.Plane)
            # arrow_pos = lift(P->fill(Point3f0(P.p), 3), P)
            # arrow_dir = lift(P->10.0*[P.u, P.v, P.w], P)
            # arrows!(p, arrow_pos, arrow_dir, arrowsize=0.1, arrowlength=10, linewidth=3)

            BitArray(Lattices.euclidean_dist(SVector{3}(I), P) <= 1/2 for I in indices)
            #@show size(B)

        else
            fill(true, sz[])
        end
    end

    mp = midpoint(state[].lattice)
    _BR = lift(p[:radius], indices) do r, indices
        if !iszero(r)
            BitArray(Lattices.dist(state[].lattice, CartesianIndex(I), mp) <= r for I in indices)
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

    positions = @lift Point3f0.(Tuple.(findall(x->x==true, $flt)))

    # colors = lift( (color, state, flt)->color(state)[reshape(flt, length(flt))],
    #  	p[:color], state, flt)
    color_func = p[:color]
    colors = lift(color_func, state, flt) do color_func, state, flt
        begin
            c = color_func(state)[flt]
            reshape(c, length(c))
        end
    end

    # cubemeshes = @lift begin
    # 	map(enumerate($positions[$flt])) do (j,p)
    # 		m = normal_mesh(Rect(Vec3f0(p), Vec3f0(1)))
    # 		pointmeta(m; color=fill(RGBAf0($colors[j]), length(coordinates(m))))
    # 	end
    # end
    # mesh_to_draw = @lift merge($cubemeshes)
    # mesh!(p, mesh_to_draw)
    cube = Rect(Vec3f0(0), Vec3f0(1))
    meshscatter!(p, positions, marker=cube, markersize=1,
          color=colors, raw=false)

    p
end






end # module
