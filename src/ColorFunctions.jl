module ColorFunctions
	using Colors
	import ..Lattices: intersect_lattice_with_plane, dimension, density
	using Graphs

	is_leaf(g::SimpleDiGraph, s) = !mapreduce(|, vertices(g)) do in; has_edge(g, in, s) end

	mutable struct ColorMapping{T}
		inner::Vector{T}
		full::Dict{T, <:Colorant}
	end
	ColorMapping(seed::Dict{T, <:Colorant}) where T = ColorMapping(collect(keys(seed)), copy(seed))
	min_colors(palette) = ColorMapping(Dict([0=>colorant"transparent", 1 => palette[1]]))

	default_palette = distinguishable_colors(256, colorant"midnightblue")

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

	function color_phylo!(C::ColorMapping{T}, state; depth=2,
		 palette=default_palette) where T

		P = state.phylogeny

		## Find all vertices which are at most `depth` removed from the root.
		## Each of which is going to be assigned its own color.
		## They form the 'inner' part.
		if depth == -1
			depth = nv(P)
		end
		ring = neighborhood(P, 1, depth, dir=:in)
		# How many colors are currently in use?
		j = length(C.inner)

		### Prune inner colors ###
		inner_colors = filter(x->x[1] in C.inner, pairs(C.full))
		C.full = Dict(inner_colors)

		# Assign colors to new genotypes within the inner ring.
		for ringv in ring
			g = state.meta[ringv, :genotype]
			if !in(g, C.inner)
				push!(C.inner, g)
				C.full[g] = current_color = palette[j%length(palette)+1]
			else
				current_color = C.full[g]
			end
			# Fill in descendants with the same color
			for child in neighborhood(P, g, nv(P), dir=:in)
				g_child = state.meta[child, :genotype]
				C.full[g_child] = current_color
			end
			j += 1
		end

		map(state.lattice.data) do s
			C.full[s]
		end
	end
	color_phylo(state; depth=2, palette=default_palette, kwargs...) = color_phylo!(min_colors(palette), state; depth, palette, kwargs...)

	function color_phylo_by_genomic_distance!(C::ColorMapping, state; depth=1,
		 palette=reverse(colormap("Blues", 20))
		)
		P = state.phylogeny

		C.full = copy(C.inner)

		for j in 2:length(state.meta.genotypes)
			l = length(state.meta.snps[j])
			# push!(C.full, state.meta.genotypes[j] => palette[l%length(palette)+1])
			push!(C.full, state.meta.genotypes[j] => palette[min(l, 20)] )
		end

		map(state.lattice.data) do s
			C.full[s]
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
