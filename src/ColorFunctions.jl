module ColorFunctions
	using Colors
	import ..Lattices: intersect_lattice_with_plane, dimension
	using LightGraphs

	is_leaf(g::SimpleDiGraph, s) = !mapreduce(|, vertices(g)) do in; has_edge(g, in, s) end

	mutable struct ColorMapping
		inner::Dict{Any, Color}
		full::Dict{Any, Color}
	end
	ColorMapping(seed::Dict{<:Any, <:Color}) = ColorMapping(seed, copy(seed))
	min_colors(palette) = ColorMapping(Dict([0 =>colorant"white", 1 => palette[1]]))

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

	function color_phylo!(C::ColorMapping, state; depth=1,
		 palette=default_palette, transition=(x::Color, d)->x)

		P = state.phylogeny

		## Find all vertices which are at most `depth` removed from the root.
		## Each of which is going to be assigned its own color.
		## They form the 'inner' part.
		if depth == -1
			depth = nv(P)
		end
		ring = neighborhood(P, 1, depth, dir=:in)
		# How many colors are currently in use?
		j = length(unique(values(C.inner))) - 1

		### Prune inner colors ###
		# Assign genotypes to indices.
		inner_gs = map(x->state.meta.genotypes[x], ring)
		# If a genotype is no longer present, remove it from the color-map.
		for k in keys(C.inner)
			if k!=0 && !(k in inner_gs)
				delete!(C.inner, k)
				j -= 1
				# @info "$k pruned"
			end
		end
		# Assign colors to new genotypes within the inner ring.
		for ringv in ring
			g = state.meta.genotypes[ringv]
			if !haskey(C.inner ,g)
				push!(C.inner, g => palette[j%length(palette)+1])
			end
			j += 1
		end

		### Build a new _full_ color-map from the inner colors.
		C.full = copy(C.inner)
		# @info C.inner
		for ringv in ring
			nb = neighborhood(P, ringv, nv(P), dir=:in)
			g_parent = state.meta.genotypes[ringv]
			color_range = range(C.inner[g_parent], colorant"white", length=round(Int, 1.5*length(nb)))
			for (d,n) in enumerate(nb)
				g = state.meta.genotypes[n]
				newcolor = color_range[d]
				#oc = RGBA(C.inner[g_parent])
				#newcolor = RGBA(oc.r, oc.g, oc.b, clamp( oc.alpha*(1-0.3), 0., 1.))
				push!(C.full, g => newcolor)
			end
		end

		map(state.lattice.data) do s
			C.full[s]
		end
	end
	color_phylo(state; depth=1, palette=default_palette, kwargs...) = color_phylo!(min_colors(palette), state; depth=depth, kwargs...)

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
