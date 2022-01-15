### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ a47c3e6d-9687-4bb9-b71f-3ee242ede575
begin
	# instantiate environment
	using Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
end

# ╔═╡ fc7b9a21-34dd-40ba-9d57-9f785904e307
begin
	# packages used in this notebook
	using GeoStats
	
	# default plot settings
	using Plots; gr(size=(700,400))
	
	# make sure that results are reproducible
	using Random; Random.seed!(2021)
end;

# ╔═╡ e72d525c-5d3b-4938-9664-e6ea9f055d4b
md"""
# Directional variograms

In this tutorial, we demonstrate how empirical variograms can be computed along a specific direction (a.k.a. *directional variograms*). This computation is extremely efficient in GeoStats.jl and can be performed in less than a second with many thousands of geospatial data points arranged in arbitrary locations.
"""

# ╔═╡ b3fa8609-b27e-4ad2-9836-b7841d6a8db0
md"""
## Synthetic data

In order to illustrate the functionality, we create synthetic data with predefined anisotropy. Later on in this tutorial, we will pretend that we don't know the anisotropy ratio and will then try to discover it with directional variograms. Here, we generate Gaussian realizations as follows.

First, we consider an anisotropic variogram model with ranges ``30`` and ``10`` aligned with the horizontal and vertical directions respectively:
"""

# ╔═╡ 342bc591-d16e-4e73-809c-6fbebdc90f0d
γ = GaussianVariogram(MetricBall((30.,10.)))

# ╔═╡ c12699ae-6b89-470e-82f7-10db271b7d2e
md"""
With this anisotropic model, we generate $3$ realizations using Gaussian simulation:
"""

# ╔═╡ 1c49fba0-1e11-11eb-100f-b5319133da0c
ensemble = let
	problem = SimulationProblem(CartesianGrid(100,100), :Z=>Float64, 3)
	
	solver  = LUGS(:Z=>(variogram=γ,))
	
	solve(problem, solver)
end

# ╔═╡ 588a7c36-82a2-494d-86aa-f0e322aa88e0
md"""
We observe that the "blobs" in the realizations are indeed stretched horizontally, and that they approximately fit into bounding boxes of size $30\times 10$:
"""

# ╔═╡ c63a6768-e680-451b-89aa-8f433f0afc82
plot(ensemble, size=(700,250))

# ╔═╡ 217cfa52-3b2a-41b2-b7ab-9f8edbf539e8
md"""
We will now use one of these realizations as our geospatial data, and will pretend that we don't know the anisotropy ratio of $30 / 10 = 3$:
"""

# ╔═╡ e6c873d1-784a-4741-aec2-1954f49bef64
𝒮 = ensemble[1]

# ╔═╡ ea6cf614-c5b9-4b87-a9e3-3055d73f5266
md"""
## Variogram calculations

We estimate the horizontal and vertical variograms from the data:
"""

# ╔═╡ 6454bed0-1a0b-11eb-0cf9-e1009c77d4d8
γhor = DirectionalVariogram((1.,0.), 𝒮, :Z, maxlag=50.)

# ╔═╡ 689fd600-1a0b-11eb-06d6-3533ef872c88
γver = DirectionalVariogram((0.,1.), 𝒮, :Z, maxlag=50.)

# ╔═╡ 72583d42-1a0b-11eb-2f2f-89224d7a7f22
begin
	plot(γhor, showbins=false, label="horizontal")
	plot!(γver, showbins=false, label="vertical")
end

# ╔═╡ 271a2412-5932-41da-a492-7351e638eee3
md"""
The plot shows clearly that the horizontal and vertical ranges are approximately $30$ and $10$, which is a satisfactory result (we know the correct answer). Finally, we can fit theoretical variogram models and estimate the anisotropy ratio:
"""

# ╔═╡ cf1e0110-1a09-11eb-38c1-a76c2f8d44e4
γₕ = fit(GaussianVariogram, γhor)

# ╔═╡ ccd0b5b0-1a09-11eb-0e34-59f93b46d9f9
γᵥ = fit(GaussianVariogram, γver)

# ╔═╡ edc4ba62-8b57-4855-9a9d-58ca1d91c393
begin
	plot(γₕ, label="horizontal")
	plot!(γᵥ, label="vertical")
end

# ╔═╡ 73def5f0-1279-4d9c-9f64-cedee35d5148
ratio = range(γₕ) / range(γᵥ)

# ╔═╡ a573757a-8f6c-4d3e-a953-e8afa93cd55c
md"""
## Varioplanes

In the previous section, we focused our attention to two major directions $(1,0)$ and $(0,1)$, and noticed that the  horizontal direction $(1,0)$ had a larger range when compared to the vertical direction $(0,1)$. Alternatively, we can investigate all possible directions in a single plot known as the *varioplane*.

In this plot, we compute the empirical variogram for all angles $\theta \in [0,2\pi]$, and optionally estimate the ranges:
"""

# ╔═╡ bec5961e-1a09-11eb-2d0c-6f76d78a7a16
γₚ = EmpiricalVarioplane(𝒮, :Z, maxlag=50.)

# ╔═╡ 5d623878-053d-4001-84ee-4ac91af22202
plot(γₚ, size=(500,500))

# ╔═╡ f344409c-343e-4569-9536-b7611e628730
md"""
As can be seen from the plot, the major direction of correlation is horizontal $range(0^\circ) \approx 30$ and the minor direction of correlation is vertical $range(90^\circ) \approx 10$. All other directions have ranges in between these two extremes $10 \le range(\theta) \le 30$.
"""

# ╔═╡ b20de14c-39cf-480a-abd6-7107518d492b
md"""
## Remarks

- Directional variograms can be computed very efficiently in GeoStats.jl with any geospatial data (e.g. point set data, grid data)

- They are useful to estimate anisotropy, particularly when a clear image is not available showing "blobs", but only sparse samples

- Variogram plane plots give a good overview of the ranges as a function of direction, for all directions in a plane
"""

# ╔═╡ Cell order:
# ╟─a47c3e6d-9687-4bb9-b71f-3ee242ede575
# ╠═fc7b9a21-34dd-40ba-9d57-9f785904e307
# ╟─e72d525c-5d3b-4938-9664-e6ea9f055d4b
# ╟─b3fa8609-b27e-4ad2-9836-b7841d6a8db0
# ╠═342bc591-d16e-4e73-809c-6fbebdc90f0d
# ╟─c12699ae-6b89-470e-82f7-10db271b7d2e
# ╠═1c49fba0-1e11-11eb-100f-b5319133da0c
# ╟─588a7c36-82a2-494d-86aa-f0e322aa88e0
# ╠═c63a6768-e680-451b-89aa-8f433f0afc82
# ╟─217cfa52-3b2a-41b2-b7ab-9f8edbf539e8
# ╠═e6c873d1-784a-4741-aec2-1954f49bef64
# ╟─ea6cf614-c5b9-4b87-a9e3-3055d73f5266
# ╠═6454bed0-1a0b-11eb-0cf9-e1009c77d4d8
# ╠═689fd600-1a0b-11eb-06d6-3533ef872c88
# ╠═72583d42-1a0b-11eb-2f2f-89224d7a7f22
# ╟─271a2412-5932-41da-a492-7351e638eee3
# ╠═cf1e0110-1a09-11eb-38c1-a76c2f8d44e4
# ╠═ccd0b5b0-1a09-11eb-0e34-59f93b46d9f9
# ╠═edc4ba62-8b57-4855-9a9d-58ca1d91c393
# ╠═73def5f0-1279-4d9c-9f64-cedee35d5148
# ╟─a573757a-8f6c-4d3e-a953-e8afa93cd55c
# ╠═bec5961e-1a09-11eb-2d0c-6f76d78a7a16
# ╠═5d623878-053d-4001-84ee-4ac91af22202
# ╟─f344409c-343e-4569-9536-b7611e628730
# ╟─b20de14c-39cf-480a-abd6-7107518d492b
