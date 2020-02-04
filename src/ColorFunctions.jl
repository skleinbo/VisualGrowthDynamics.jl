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
