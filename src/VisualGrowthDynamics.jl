module VisualGrowthDynamics

import Base.Iterators: product
import StaticArrays: SVector

import GrowthDynamics
import GrowthDynamics: TumorConfigurations, Lattices
import GrowthDynamics.Lattices: midpoint, midpointcoord
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

@recipe(TumorPlot, state) do scene
    merge(DEFAULT_ATTRIBUTES,
        Attributes(
            plane = nothing,
            genotype = nothing,
            filter_radius = false,
            midpoint = nothing,
            radius = 0.0,
            shell = nothing
        )
    )
end
function Makie.plot!(p::TumorPlot{Tuple{<:TumorConfiguration{A}}}) where A<:RealLattice{Int}
    # @show typeof(p)
    state = p[:state]
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

    p[:midpoint][] = midpointcoord(state[].lattice)

    a = spacings(state[].lattice)[1]
    _BR = lift(p[:filter_radius], p[:radius], p[:midpoint], indices) do bflt_r, r, mp, indices
        if bflt_r && !iszero(r)
            idx_mp = Tuple(Lattices.index(state[].lattice, mp))
            BitArray(r-a/2 <= Lattices.dist(state[].lattice, CartesianIndex(I), idx_mp) < r+a/2 for I in indices)
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

function validate_xyz(s::AbstractString)
    v = split(s, ',')
    if length(v)==3 && all(x->!isnothing(tryparse(Float64, x)), v)
        return true
    end
    return false
end

xyz_string_to_point3(s::AbstractString) = Point3(tryparse.(Float64, split(s, ',')))

function TumorInspector(state::TumorConfiguration{A}, args...; kwargs...) where A<:RealLattice{Int}
    fig = Figure()
    grid_plots = fig[1,1] = GridLayout()
    grid_controls = fig[1,2] = GridLayout(;tellheight=false)
    tumor_ax, tumor_plot = tumorplot(grid_plots[1,1], Observable(state), args...; kwargs...)
    
    ## Radius filter controls
    Label(grid_controls[1,1], text="Filter radius", tellheight=false)
    toggle_radius = Toggle(grid_controls[1,2])
    Label(grid_controls[2,1], text="Radius")
    text_radius = Textbox(grid_controls[2,2], placeholder="99", validator=Float64)
    btn_radius_plus = Button(grid_controls[2,3], label="+")
    btn_radius_minus = Button(grid_controls[2,4], label="-")
    on(btn_radius_plus.clicks) do _
        r = (tumor_plot[:radius][] += 1)
        text_radius.stored_string[] = string(r)
        text_radius.displayed_string[] = string(r)
    end
    on(btn_radius_minus.clicks) do _
        r = max(0, (tumor_plot[:radius][] - 1))
        tumor_plot[:radius][] = r
        text_radius.stored_string[] = string(r)
        text_radius.displayed_string[] = string(r)
    end
    Label(grid_controls[3,1], text="Midpoint")
    text_mp = Textbox(grid_controls[3,2], placeholder="x,y,z", validator=validate_xyz)
    btn_reset_mp = Button(grid_controls[3,3], label="Reset")
    connect!(tumor_plot[:filter_radius], toggle_radius.active)
    on(text_radius.stored_string) do s
        @info "radius changed to $s"
        tumor_plot[:radius][] = parse(Float64, s)
    end
    on(text_mp.stored_string) do s
        @info "MP changed to $(xyz_string_to_point3(s))"
        @show typeof(xyz_string_to_point3(s))
        tumor_plot[:midpoint][] = xyz_string_to_point3(s)
    end
    on(btn_reset_mp.clicks) do _
        @info "Reset btn clicked"
        v = Lattices.midpointcoord(tumor_plot[:state][].lattice)
        #tumor_plot[:midpoint][] = v
        text_mp.stored_string[] = join(v,',')
        text_mp.displayed_string[] = text_mp.stored_string[]
    end

    ## Plane filter controls
    Label(grid_controls[4,1], text="Plane direction")
    btn_plane_x = Button(grid_controls[4,2], label="x")
    btn_plane_y = Button(grid_controls[4,3], label="y")
    btn_plane_z = Button(grid_controls[4,4], label="z")
    # Label(grid_controls[5,1], text="Origin")
    Label(grid_controls[5,1], text="Offset")
    text_plane_offset = Textbox(grid_controls[5,2], placeholder="0", validator=Float64)
    

    return fig
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
