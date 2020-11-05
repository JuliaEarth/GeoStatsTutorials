### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ f2ea0ebb-8502-4c0b-99fb-dbc673e564f7
using Pkg; Pkg.instantiate(); Pkg.precompile()

# ╔═╡ 6585c350-1a90-11eb-174a-274f60f662d9
begin
	using GeoStats
	using PlutoUI
	using Plots
	gr(size=(700,400),
	   xlim=(0,1),
	   ylim=(0,1))
end;

# ╔═╡ 0e4ece36-520d-4a55-9057-0cf3c24754dd
using Random; Random.seed!(2000);

# ╔═╡ 6c9b41d8-043f-4626-8798-4a3ae9da922b
md"""
# Gaussian processes

Gaussian process *regression* and Simple Kriging are essentially [two names for the same concept](https://en.wikipedia.org/wiki/Kriging). The derivation of Kriging estimators, however; does **not** require distributional assumptions. It is a beautiful coincidence that simple Kriging estimates are equal to the conditional mean of a Gaussian process with the same mean and covariance.

Gaussian process regression is usually explained in 1D, and we will do the same here.
"""

# ╔═╡ b069e90a-5be6-4ef9-8063-d6271fa3b386
md"""
## The data

Let's start by drawing a set of points in the real line, each with a known value:
"""

# ╔═╡ ce127b5b-4cfd-4b00-8ee2-f73bf756de19
begin
	x = range(0.1, stop=.9, length=10)
	z = rand(10)
	
	scatter(x, z, label="data")
end

# ╔═╡ ba9c726e-1f70-11eb-3d63-ab3c21c143d3
md"""
And a set of unseen locations where to make predictions:
"""

# ╔═╡ 8dad2424-1f70-11eb-3dfe-d5719a07e959
begin
	n  = 200
	xs = range(0., stop=1, length=n)
end

# ╔═╡ d358e5bb-aaac-485c-a952-b29170861053
md"""
## Kriging estimators

GeoStats.jl provides various types of [Kriging estimators](https://juliaearth.github.io/GeoStats.jl/stable/kriging/estimators). We will demonstrate two of these types: `SimpleKriging` and `OrdinaryKriging`.

For `SimpleKriging`, we need to specify the covariance (or variogram) function and the mean of the stochastic process. For `OrdinaryKriging`, we only need to specify the covariance as the mean is derived from linear constraints.

Let's create a variogram model and vary the sill and range to understand the effect on the estimates:
"""

# ╔═╡ c2e80350-1a90-11eb-0a0a-ddf37ddac7d8
md"sill = $(@bind s Slider(range(.0, stop=.1, length=50); default=0.05, show_value=true))"

# ╔═╡ c3966bc0-1a90-11eb-05b3-0358a9025603
md"range = $(@bind r Slider(range(.1, stop=.2, length=50); default=0.15, show_value=true))"

# ╔═╡ 2b90a50a-1f74-11eb-3375-bf5b937a963a
γ = GaussianVariogram(sill=Float64(s), range=Float64(r))

# ╔═╡ 3f2b8950-1a91-11eb-2312-5f961fd3ac19
begin
	# Kriging estimators
	sk = SimpleKriging(x', z, γ, 0.5)
	ok = OrdinaryKriging(x', z, γ)
	
	#conditional mean and variance
	μ1, μ2 = Float64[], Float64[]
	σ1, σ2 = Float64[], Float64[]
	
	for xₒ in xs
		μ₁, σ₁² = predict(sk, [xₒ])
		μ₂, σ₂² = predict(ok, [xₒ])
		
		push!(μ1, μ₁)
		push!(μ2, μ₂)
		push!(σ1, √σ₁²)
		push!(σ2, √σ₂²)
	end
	
	# plot results
	p1 = scatter(x, z, label="data")
	plot!(xs, μ1, ribbon=σ1, fillalpha=.5, label="Simple Kriging")
	p2 = scatter(x, z, label="data")
	plot!(xs, μ2, ribbon=σ2, fillalpha=.5, label="Ordinary Kriging")
	
	plot(p1, p2, layout=(2,1))
end

# ╔═╡ dea78a79-903f-4157-a6a8-924c2a82aa31
md"""
# Other generalizations

- We can also have random fields with non-constant mean. In Universal Kriging, the mean is a polynomial on the coordinates of the problem up to a certain degree specified by the user. In External Drift Kriging, the mean can be any external variable (a.k.a. external drift).

- All methods in GeoStats.jl accept general distance functions. We can for example use the Haversine distance to compute covariances between latitude/longitude coordinates directly.
"""

# ╔═╡ Cell order:
# ╟─f2ea0ebb-8502-4c0b-99fb-dbc673e564f7
# ╠═6585c350-1a90-11eb-174a-274f60f662d9
# ╠═0e4ece36-520d-4a55-9057-0cf3c24754dd
# ╟─6c9b41d8-043f-4626-8798-4a3ae9da922b
# ╟─b069e90a-5be6-4ef9-8063-d6271fa3b386
# ╠═ce127b5b-4cfd-4b00-8ee2-f73bf756de19
# ╟─ba9c726e-1f70-11eb-3d63-ab3c21c143d3
# ╠═8dad2424-1f70-11eb-3dfe-d5719a07e959
# ╟─d358e5bb-aaac-485c-a952-b29170861053
# ╟─c2e80350-1a90-11eb-0a0a-ddf37ddac7d8
# ╟─c3966bc0-1a90-11eb-05b3-0358a9025603
# ╠═2b90a50a-1f74-11eb-3375-bf5b937a963a
# ╠═3f2b8950-1a91-11eb-2312-5f961fd3ac19
# ╟─dea78a79-903f-4157-a6a8-924c2a82aa31
