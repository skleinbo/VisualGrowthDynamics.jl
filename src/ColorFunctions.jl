module ColorFunctions
	using Colors
	import ..Lattices: intersect_lattice_with_plane, dimension, density
	using Graphs
	import GrowthDynamics.TumorConfigurations: gindex, TumorConfiguration

	is_leaf(g::SimpleDiGraph, s) = !mapreduce(|, vertices(g)) do in; has_edge(g, in, s) end

	mutable struct ColorMapping{T}
		inner::Vector{T}
		full::Dict{T, <:Colorant}
	end
	ColorMapping(seed::Dict{T, <:Colorant}) where T = ColorMapping(collect(keys(seed)), copy(seed))
	min_colors(palette) = ColorMapping(Dict([0=>colorant"transparent"]))

	default_palette = distinguishable_colors(256, Colors.JULIA_LOGO_COLORS[:blue])

	"""
		lighten(C, p)

		Lighten a color `C` by `p` percent. Internally converts to HSV colorspace,
		reduces saturation, increases value (brightness) and converts back to
		the orginal colorspace.
	"""
	function lighten(C,p)
	    T = typeof(C)
	    hsv = HSV(C)
	    convert(T, HSV(hsv.h, min(hsv.s * (1-0*p), 1.0), min(hsv.v * (1+p), 1.0) ))
	end


	"""
		color_lineages!(C, state; roots=[1])
	
	Assign to each genotype in `roots` a distinct color from `palette`;
	assign all children the same color as their root.
	Previous assignments are kept.

	Expects a `ColorMapping` as its first argument.

	# Additional keyword arguments
	* `palette`: vector of colors for the root nodes.
	"""
	function color_lineages!(C::ColorMapping{T}, state;
		 roots=[one(T)],
		 palette=default_palette,
		 default=colorant"transparent",
		 kwargs...) where T

		P = state.phylogeny

		# Indices of root nodes
		iroots = gindex.(Ref(state.meta), roots)
		# no. of root colors in use
		j = length(C.inner)
		# root colors
		inner_colors = filter(x->x[1] in C.inner, pairs(C.full))
		C.full = Dict(inner_colors)

		# Assign colors to new genotypes within the inner ring.
		for (root, iroot) in zip(roots, iroots)
			if !in(root, C.inner)
				push!(C.inner, root)
				C.full[root] = current_color = palette[(j-1)%length(palette)+1]
			else
				current_color = C.full[root]
			end
			# Fill in descendants with the same color
			for child in neighborhood(P, iroot, nv(P), dir=:in)
				g_child = state.meta[child, :genotype]
				C.full[g_child] = current_color
			end
			j += 1
		end

		map(state.lattice.data) do s
			get(C.full, s, default)
		end
	end
	color_lineages(state; palette=default_palette, kwargs...) = color_lineages!(min_colors(palette), state; kwargs...)
	
	function color_depth!(C::ColorMapping, state; depth=2, kwargs...)
		inner = neighborhood(state.phylogeny, 1, depth, dir=:in)
		ginner = state.meta[inner, :genotypes]
		color_lineages!(C, state; roots=ginner, kwargs...)
	end
	color_depth(state; depth=2, palette=default_palette, kwargs...) = color_depth!(min_colors(palette), state; depth, palette, kwargs...)

	function color_phylo_by_genomic_distance!(C::ColorMapping, state; depth=1,
		 palette=reverse(colormap("Blues", 20))
		)
		P = state.phylogeny

		inner_colors = filter(x->x[1] in C.inner, pairs(C.full))
		C.full = Dict(inner_colors)

		for j in 2:length(state.meta)
			l = length(state.meta[j, :snps])
			# push!(C.full, state.meta.genotypes[j] => palette[l%length(palette)+1])
			push!(C.full, state.meta[j, :genotype] => palette[min(l, 20)] )
		end

		map(state.lattice.data) do s
			get(C.full, s , default)
		end
	end
	color_phylo_by_genomic_distance(state; depth=1, palette=default_palette, kwargs...) = color_phylo_by_genomic_distance!(min_colors(palette), state; depth, palette, kwargs...)

	function color_generic(state;palette=rand(RGB,maximum(state.lattice.data)))
		v0 = state.lattice.data
		map(v0) do x
			if x==0
				RGBA{Float32}(1.0,1.0,1.0,0.0)
			else
				RGBA{Float32}(palette[x])
			end
		end
	end

	function color_genotype(state, g; color=colorant"blue")
		v0 = state.lattice.data
		color = RGBA(color)
		map(v0) do x
			if x!=g
				RGBA{Float32}(1.0,1.0,1.0,0.0)
			else
				color
			end
		end
	end

	function color_neighbors(state, I::CartesianIndex)
		nn = neighbors(state.lattice, I)
		C = similar(state.lattice.data, RGBA{Float32})
		fill!(C, colorant"gray")

		C[nn] .= colorant"cornflowerblue"
		return reshape(C, length(C))
	end

	function color_column(state, q)
		C = similar(state.lattice.data, RGBA{Float32})
		fill!(C, colorant"gray")

		C[:,q] .= colorant"cornflowerblue"
		return reshape(C, length(C))
	end

	function color_row(state, r)
		C = similar(state.lattice.data, RGBA{Float32})
		fill!(C, colorant"gray")

		C[r,:] .= colorant"cornflowerblue"
		return reshape(C, length(C))
	end

	function color_rand(state)
		C = similar(state.lattice.data, RGBA{Float32})
		reshape(rand(RGB{Float32}, size(state.lattice.data)), length(state.lattice.data))
	end

	function hide!(v, A::BitArray)
		for i in eachindex(v)
			if !A[i]
				v[i] = RGBA{Float32}(0., 0., 0., 0.)
			end
		end
		v
	end


#= End module =#
end
