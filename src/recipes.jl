DEFAULT_ATTRIBUTES = Attributes(
    markersize=1,
    size=(800,800),
    colorfunc=ColorFunctions.color_depth
)

@recipe(TumorPlot, state) do scene
    merge(DEFAULT_ATTRIBUTES,
        Attributes(
            plane = nothing,
            genotype = nothing,
            filter_radius = false,
            midpoint = nothing,
            radius = 0.0,
            radius_eq = false,
            shell = nothing
        )
    )
end

function Makie.plot!(p::TumorPlot{Tuple{<:TumorConfiguration{G, A}}}) where {G, A<:RealLattice{Int}}
    state = p[:state]
    sz = lift(s->size(s.lattice), state)

    geno_filter = lift(state, p[:genotype]) do state, genotype
        if isnothing(genotype) || iszero(genotype)
            state.meta[:, :genotype]
        else
            state.meta[neighborhood(state.phylogeny, genotype, typemax(Int64); dir=:in), :genotype]
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
    _BR = lift(p[:filter_radius], p[:radius], p[:radius_eq], p[:midpoint], indices) do bflt_r, r, radius_eq, mp, indices
        if bflt_r && !iszero(r)
            idx_mp = Tuple(Lattices.index(state[].lattice, mp))
            rmin = radius_eq ? r-a/2 : 0
            BitArray(rmin <= Lattices.dist(state[].lattice, CartesianIndex(I), idx_mp) < r+a/2 for I in indices)
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

    positions = @lift begin 
        v = [ Lattices.coord($state.lattice, x) for x in CartesianIndices($flt) if $flt[x] ]
        if !isnothing(p[:plane][])
            v = map(x->project(p[:plane][], x), v)
        end
        v
    end

    colorfunc = p[:colorfunc]
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
    
    # Determine marker size
    markersize = norm(coord(state[].lattice, neighbors(state[].lattice, (1,1,1))[1]))/2

    meshscatter!(p, positions, marker=_mesh, markersize=markersize,
    color=color
    )

    p
end

meshvertex(::Type{HexagonalLattice}) = GeometryBasics.mesh(map([
    # Vertices of a hex with flat edge down
    Point2f(1,0),
    Point2f(cos(pi/3),-sin(pi/3)),
    Point2f(-cos(pi/3),-sin(pi/3)),
    Point2f(-1,0),
    Point2f(-cos(pi/3),sin(pi/3)),
    Point2f(cos(pi/3),sin(pi/3))
    ]/2/Float32(cos(pi/6))) do x
        # Rotation by π/6 --> pointy edge down
        [Point2f(cos(pi/6),sin(pi/6)) Point2f(-sin(pi/6),cos(pi/6))]*x
    end
    )

meshvertex(::Type{CubicLattice{A,B}}) where {A,B} = GeometryBasics.mesh(Rect(Vec3f(-1), Vec3f(2)))
# TODO: Vertex for FCC lattice
# meshvertex(::Type{FCCLattice{A,B}}) where {A,B} = GeometryBasics.Sphere(Point3f(0), 1/4)
meshvertex(::Type{FCCLattice{A,B}}) where {A,B} = fccvertex