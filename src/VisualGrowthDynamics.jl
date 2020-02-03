module VisualGrowthDynamics

export plot!, plot

import GrowthDynamics: TumorConfigurations, Lattices
import .TumorConfigurations: TumorConfiguration
import .Lattices: HexagonalLattice
import AbstractPlotting: Plot, plot!, default_theme, to_value

using GeometryTypes
using Makie, GLMakie

function default_theme(scene::SceneLike, ::Type{<: Plot(TumorConfiguration{HexagonalLattice{Int64}})})
	Theme(
		markersize=1
	)
end

function plot!(p::Plot(TumorConfiguration{HexagonalLattice{Int64}}); color)
	state = to_value(p[1])
	L = state.lattice.Na
	@show color

    hex_vertices = map(x->[Point2f0(cos(pi/6),sin(pi/6)) Point2f0(-sin(pi/6),cos(pi/6))]*x, [Point2f0(1,0),
		    Point2f0(cos(pi/3),-sin(pi/3)),
		    Point2f0(-cos(pi/3),-sin(pi/3)),
		    Point2f0(-1,0),
		    Point2f0(-cos(pi/3),sin(pi/3)),
		    Point2f0(cos(pi/3),sin(pi/3))]/L/Float32(cos(pi/6)))
	hex_mesh=GLNormalMesh(hex_vertices)
	off = (sqrt(5)-sqrt(3))/2

	hex_positions = reshape(Point2f0[Point2f0(2*(xi+ifelse(yi%2==0,0.,1/2))/L-1.,
		2*(yi+1/2+off/2*ifelse(yi==0,0.,-yi))/L-1.)
		 for (xi,yi) in Iterators.product(0:L-1,0:L-1)], L^2)

	meshscatter!(p, hex_positions, marker=hex_mesh, markersize=1, color=lift(color, state))
	p
end



end # module
