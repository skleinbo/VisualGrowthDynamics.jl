function color_gen(state;palette=rand(RGB,maximum(state.lattice.data)))
	v0 = state.lattice.data
	map(v0) do x
		if x==0
			RGBA{Float32}(1.0,1.0,1.0,1.0)
		else
			RGBA{Float32}(palette[x])
		end
	end |> X->reshape(X, length(v0))
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
