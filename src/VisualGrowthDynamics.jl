module VisualGrowthDynamics


# export plot!, plot

import GrowthDynamics: TumorConfigurations, Lattices
import .TumorConfigurations: TumorConfiguration
import .Lattices: HexagonalLattice, neighbors
import AbstractPlotting: Plot, plot!, default_theme, to_value

using GeometryTypes
using Colors
using Makie, GLMakie

include("ColorFunctions.jl")

function default_theme(scene::SceneLike, ::Type{<: Plot(TumorConfiguration{HexagonalLattice{Int64}})})
	Theme(
		markersize=1,
		size=(800,800)
	)
end

function plot!(p::Plot(TumorConfiguration{HexagonalLattice{Int64}}))
	state = to_value(p[1])
	L = state.lattice.Na
	@show to_value(p[:color])

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
	hex_mesh=GLNormalMesh(hex_vertices)
	off = (sqrt(5)-sqrt(3))/2

	hex_positions = reshape(Point2f0[Point2f0(
		2*(q+ifelse(r%2==0,0.,1/2))/L-1.0,
		2*(r+1/2+off/2*ifelse(r==0,0.,-r))/L-1.0)
		 for (r,q) in Iterators.product(0:L-1,0:L-1)], L^2)

	meshscatter!(p, hex_positions, marker=hex_mesh, markersize=1,
	 	color=lift(to_value(p[:color]), Node(state)),
		shading=false)
	p
end



end # module
