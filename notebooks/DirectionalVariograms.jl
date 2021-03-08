### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 2b9ebaaa-1f9c-11eb-39ff-a5ee1ecb90ee
begin
	using Distributed
	pids = [myid()]
	
	md"""
	Running on processes: $pids
	
	Use `pids = addprocs(n)` to run the notebook with `n` parallel processes.
	"""
end

# ╔═╡ a47c3e6d-9687-4bb9-b71f-3ee242ede575
@everywhere pids begin
	using Pkg; Pkg.activate(@__DIR__)
	Pkg.instantiate(); Pkg.precompile()
end

# ╔═╡ fc7b9a21-34dd-40ba-9d57-9f785904e307
@everywhere pids begin
	# packages used in this notebook
	using GeoStats
	
	# default plot settings
	using Plots; gr(size=(700,400))
	
	# make sure that results are reproducible
	using Random; Random.seed!(2021)
end

# ╔═╡ e72d525c-5d3b-4938-9664-e6ea9f055d4b
md"""
# Directional variograms

In this tutorial, we demonstrate how empirical variograms can be computed along a specific direction (a.k.a. *directional variograms*). This computation is extremely efficient in GeoStats.jl and can be performed in less than a second with many thousands of spatial data points arranged in arbitrary locations.
"""

# ╔═╡ b3fa8609-b27e-4ad2-9836-b7841d6a8db0
md"""
## Synthetic data

In order to illustrate the functionality, we create synthetic data with predefined anisotropy. Later on in this tutorial, we will pretend that we don't know the anisotropy ratio and will then try to discover it with directional variograms. Here, we generate Gaussian realizations as follows.

First, we consider a base (isotropic) variogram model:
"""

# ╔═╡ fe880642-1a08-11eb-1a6d-89cd44643007
GaussianVariogram(range=10.)

# ╔═╡ f33cd7c5-7e57-4962-8e30-2e7155c01484
md"""
To convert this isotropic model into an anisotropic model, we use an ellipsoid distance. Given that the range of the variogram is $10$, we will stretch the $x$ axis by a factor of $3$ to produce an effective horizontal range of $30$, and the $y$ axis by a factor of $1$, which will leave the vertical range untouched and equal to $10$. We set the angle of the ellipsoid to $0$ so that the anisotropy is aligned with the coordinate system:
"""

# ╔═╡ 342bc591-d16e-4e73-809c-6fbebdc90f0d
γ = GaussianVariogram(range=10., distance=aniso2distance([3.,1.], [0.]))

# ╔═╡ c12699ae-6b89-470e-82f7-10db271b7d2e
md"""
With this anisotropic model, we generate $3$ realizations using direct Gaussian simulation:
"""

# ╔═╡ 1c49fba0-1e11-11eb-100f-b5319133da0c
problem = SimulationProblem(CartesianGrid(100,100), :Z=>Float64, 3)

# ╔═╡ 19fcfe60-1e11-11eb-366f-b745e9a95425
solver  = LUGS(:Z=>(variogram=γ,))

# ╔═╡ 24a088f0-1e11-11eb-35e6-f17264f4dcba
solution = solve(problem, solver)

# ╔═╡ 588a7c36-82a2-494d-86aa-f0e322aa88e0
md"""
We observe that the "blobs" in the realizations are indeed stretched horizontally, and that they approximately fit into bounding boxes of size $30\times 10$:
"""

# ╔═╡ c63a6768-e680-451b-89aa-8f433f0afc82
plot(solution, size=(700,250))

# ╔═╡ 217cfa52-3b2a-41b2-b7ab-9f8edbf539e8
md"""
We will now use one of these realizations as our spatial data, and will pretend that we don't know the anisotropy ratio of $30 / 10 = 3$:
"""

# ╔═╡ e6c873d1-784a-4741-aec2-1954f49bef64
𝒮 = solution[1]

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

- Directional variograms can be computed very efficiently in GeoStats.jl with any spatial data (e.g. point set data, regular grid data)

- They are useful to estimate anisotropy, particularly when a clear image is not available showing "blobs", but only sparse samples

- Variogram plane plots give a good overview of the ranges as a function of direction, for all directions in a plane
"""

# ╔═╡ Cell order:
# ╟─2b9ebaaa-1f9c-11eb-39ff-a5ee1ecb90ee
# ╟─a47c3e6d-9687-4bb9-b71f-3ee242ede575
# ╠═fc7b9a21-34dd-40ba-9d57-9f785904e307
# ╟─e72d525c-5d3b-4938-9664-e6ea9f055d4b
# ╟─b3fa8609-b27e-4ad2-9836-b7841d6a8db0
# ╠═fe880642-1a08-11eb-1a6d-89cd44643007
# ╟─f33cd7c5-7e57-4962-8e30-2e7155c01484
# ╠═342bc591-d16e-4e73-809c-6fbebdc90f0d
# ╟─c12699ae-6b89-470e-82f7-10db271b7d2e
# ╠═1c49fba0-1e11-11eb-100f-b5319133da0c
# ╠═19fcfe60-1e11-11eb-366f-b745e9a95425
# ╠═24a088f0-1e11-11eb-35e6-f17264f4dcba
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
