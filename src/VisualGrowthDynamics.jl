module VisualGrowthDynamics

import Base.Iterators: product
import StaticArrays: SVector
import GrowthDynamics: TumorConfigurations, Lattices
import .TumorConfigurations: TumorConfiguration
import .Lattices: HexagonalLattice, CubicLattice, neighbors
import LightGraphs: neighborhood

using GeometryBasics
using Colors

import AbstractPlotting: Plot, plot!, default_theme, to_value
using Makie

include("ColorFunctions.jl")


function default_theme(scene, ::Type{<: Plot(TumorConfiguration{HexagonalLattice{Int64}})})
	Attributes(
		markersize=1,
		size=(800,800)
	)
end

function plot!(p::Plot(TumorConfiguration{HexagonalLattice{Int64}}))
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
function default_theme(scene::SceneLike, ::Type{<: Plot(TumorConfiguration{CubicLattice{Int64}})})
	Attributes(
		markersize=1,
		resolution=(800,800),
		genotype=1
	)
end

function AbstractPlotting.plot!(p::Plot(TumorConfiguration{CubicLattice{Int64}}))
	global myp = p
	@show typeof(p)
	state = to_value(p[1])
	L = state.lattice.Na

	geno_filter = state.meta.genotypes[
		neighborhood(state.Phylogeny, to_value(p[:genotype]), typemax(Int64); dir=:in)
		]

	flt = BitArray( x in geno_filter for x in state.lattice.data)

	if haskey(p.attributes, :plane)
		P = to_value(p[:plane])
		#B = Lattices.intersect_lattice_with_plane(state.lattice, P)
		indices = product(1:L, 1:L, 1:L)
		B = BitArray(Lattices.euclidean_dist(SVector{3}(I), P) <= 1/2 for I in indices)
		#@show size(B)
		flt = flt .& B

		arrows!(fill(Point3f0(P.p), 3), 10.0*[P.u, P.v, P.w], arrowsize=0.1, arrowlength=10, linewidth=3)
	end

	if all(x->x==false, flt)
		@info "No cells to draw."
	end

	positions = Point3f0[Point3f0(l,m,n)
				 for (l,m,n) in Iterators.product(0:L-1,0:L-1,0:L-1)][flt]
    positions = reshape(positions, length(positions))

	color = (to_value(p[:color])(state))[reshape(flt, length(flt))]

    cube = GeometryBasics.HyperRectangle(0,0,0, 1,1,1)
	meshscatter!(p, positions, marker=cube, markersize=1,
	 	 color=color, raw=false)

	p
end






end # module
