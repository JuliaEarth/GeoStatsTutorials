### A Pluto.jl notebook ###
# v0.12.4

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

# ╔═╡ 0e4ece36-520d-4a55-9057-0cf3c24754dd
using Random; Random.seed!(2000);

# ╔═╡ 6585c350-1a90-11eb-174a-274f60f662d9
begin
	using GeoStats
	using Plots; gr(size=(900,400), xlim=(0,1), ylim=(0,1))
	using PlutoUI
end

# ╔═╡ 6c9b41d8-043f-4626-8798-4a3ae9da922b
md"""
# Gaussian processes

In this tutorial, we explain how to obtain the functionality of other packages such as GaussianProcesses.jl with GeoStats.jl. Gaussian process regression and Simple Kriging are essentially [two names for the same concept](https://en.wikipedia.org/wiki/Kriging). The derivation of Kriging estimators, however; does **not** require distributional assumptions. [Matheron](https://en.wikipedia.org/wiki/Georges_Matheron) and other important geostatisticians have generalized Gaussian processes to random fields with locally-varying mean and for situations where the mean is unknown. GeoStats.jl includes Gaussian process regression as a special case as well as other more practical Kriging variants.

In machine learning, Gaussian processes are usually explained in 1D, and we will do the same here.
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

# ╔═╡ d358e5bb-aaac-485c-a952-b29170861053
md"""
We will skip an important step known in geostatistics as *variogram modeling*, and will proceed to the estimation step directly in order to avoid having too many concepts in this tutorial. These concepts are likely new to practioners of machine learning and can be very useful in practice. For more information on variogram modeling, please check the other tutorials in this folder.

## Simple Kriging

In simple Kriging, we need to specify a mean and a covariance-type function in order to interpolate the stochastic process at unseen locations. We choose a covariance (or variogram) model in GeoStats.jl and interact with the parameters:
"""

# ╔═╡ 980ffed0-1a90-11eb-1907-e506cd2dd222
md"points where to make predictions"

# ╔═╡ b8b64d60-1a90-11eb-3578-479b2690911f
begin
	n  = 200
	xs = range(0., stop=1, length=n)
end

# ╔═╡ c2e80350-1a90-11eb-0a0a-ddf37ddac7d8
md"s = $(@bind s Slider(range(.0, stop=.1, length=50); default=0.05, show_value=true))"

# ╔═╡ c3966bc0-1a90-11eb-05b3-0358a9025603
md"r = $(@bind r Slider(range(.1, stop=.2, length=50); default=0.15, show_value=true))"

# ╔═╡ 30a8efd0-1a91-11eb-3873-a7d4675c0add
md"create Kriging estimator"

# ╔═╡ 3f2b8950-1a91-11eb-2312-5f961fd3ac19
begin
	μ = 0.5
	γ = GaussianVariogram(sill=Float64(s), range=Float64(r))
	sk = SimpleKriging(x', z, γ, μ)
	
	#conditional mean and variance
	μs = zeros(n)
	σ² = zeros(n)
	ps = [predict(sk, [xₒ]) for xₒ in xs]
	μs = first.(ps)
	σs = sqrt.(last.(ps))
end

# ╔═╡ 69d2b9d0-1a91-11eb-3d57-f57a351205d7
begin
	scatter(x, z, label="data")
	plot!(xs, μs, ribbon=σs, fillalpha=.5, label="Simple Kriging")
end

# ╔═╡ f04ae1cb-4d7f-4617-a90a-ac4f74c38930
md"""
## Ordinary Kriging

An immediate generalization of Simple Kriging is Ordinary Kriging. In this case, we still have the assumption that the mean of the field is constant, but this time we treat it as an unknown. Below is the estimation with Ordinary Kriging where we only need to specify the covariance-type function:
"""

# ╔═╡ 2f827bc0-1a92-11eb-01fc-69c5cd3b335b
md"s1 = $(@bind s1 Slider(range(.0, stop=.1, length=50); default=0.05, show_value=true))"

# ╔═╡ 2d5d51d2-1a92-11eb-0c77-238f5cfe48ce
md"r1 = $(@bind r1 Slider(range(.1, stop=.2, length=50); default=0.15, show_value=true))"

# ╔═╡ 4aa706be-93d3-4e3f-86e3-a73de75cfd6e
begin
	# create Kriging estimator
	γ1 = GaussianVariogram(sill=s1, range=r1)
	ok = OrdinaryKriging(x', z, γ1)
	
	# conditional mean and variance
	μs1 = zeros(n)
	σ²1 = zeros(n)
	ps1 = [predict(ok, [xₒ]) for xₒ in xs]
	μs1 = first.(ps1)
	σs1 = sqrt.(last.(ps1))
end

# ╔═╡ 49dc4860-1a93-11eb-3981-55dd70ff76c8
begin
	scatter(x, z, label="data")
	plot!(xs, μs1, ribbon=σs1, fillalpha=.5, label="Ordinary Kriging")
end

# ╔═╡ dea78a79-903f-4157-a6a8-924c2a82aa31
md"""
# Other generalizations

- We can also have random fields with non-constant mean. In Universal Kriging, the mean is a polynomial on the coordinates of the problem up to a certain degree specified by the user. In External Drift Kriging, the mean can be any external variable (a.k.a. external drift).

- All methods in GeoStats.jl accept general distance functions. We can for example use the Haversine distance to compute covariances between latitude/longitude coordinates directly.

- Other methods for co-estimation and co-simulation will be available in future releases.
"""

# ╔═╡ Cell order:
# ╠═f2ea0ebb-8502-4c0b-99fb-dbc673e564f7
# ╠═0e4ece36-520d-4a55-9057-0cf3c24754dd
# ╠═6585c350-1a90-11eb-174a-274f60f662d9
# ╟─6c9b41d8-043f-4626-8798-4a3ae9da922b
# ╟─b069e90a-5be6-4ef9-8063-d6271fa3b386
# ╠═ce127b5b-4cfd-4b00-8ee2-f73bf756de19
# ╟─d358e5bb-aaac-485c-a952-b29170861053
# ╟─980ffed0-1a90-11eb-1907-e506cd2dd222
# ╠═b8b64d60-1a90-11eb-3578-479b2690911f
# ╟─c2e80350-1a90-11eb-0a0a-ddf37ddac7d8
# ╟─c3966bc0-1a90-11eb-05b3-0358a9025603
# ╟─30a8efd0-1a91-11eb-3873-a7d4675c0add
# ╠═3f2b8950-1a91-11eb-2312-5f961fd3ac19
# ╠═69d2b9d0-1a91-11eb-3d57-f57a351205d7
# ╟─f04ae1cb-4d7f-4617-a90a-ac4f74c38930
# ╟─2f827bc0-1a92-11eb-01fc-69c5cd3b335b
# ╟─2d5d51d2-1a92-11eb-0c77-238f5cfe48ce
# ╠═4aa706be-93d3-4e3f-86e3-a73de75cfd6e
# ╠═49dc4860-1a93-11eb-3981-55dd70ff76c8
# ╟─dea78a79-903f-4157-a6a8-924c2a82aa31
