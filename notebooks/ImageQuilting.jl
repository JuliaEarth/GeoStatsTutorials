### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 1ed84657-f480-4920-8719-d9123ca0d7d4
begin
	using Distributed
	pids = [myid()]
	
	md"""
	Running on processes: $pids
	
	Use `pids = addprocs(n)` to run the notebook with `n` parallel processes.
	"""
end

# ╔═╡ 39aebddc-1fa3-11eb-2f85-813c7473ffab
@everywhere pids begin
	using Pkg; Pkg.activate(@__DIR__)
	Pkg.instantiate(); Pkg.precompile()
end

# ╔═╡ 418e9d24-1fa3-11eb-3ef9-51b170d57b9c
@everywhere pids begin
	# packages used in this notebook
	using GeoStats
	using ImageQuilting
	using GeoStatsImages
	
	# default plot settings
	using Plots; gr(size=(700,400), c=:cividis)
	
	# make sure that results are reproducible
	using Random; Random.seed!(2000)
end

# ╔═╡ 2a154f98-20d8-4277-a582-02dd4a6824a5
md"""
# Image quilting

In this tutorial, we demonstrate the use of the image quilting solver, which is the fastest multiple-point simulation solver in the literature. For more information about the algorithm and its performance, please watch:
"""

# ╔═╡ 91dca1b8-1fa3-11eb-3896-7fd43f394a02
html"""
<iframe width="560" height="315" src="https://www.youtube.com/embed/YJs7jl_Y9yM" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>"""

# ╔═╡ 3389c9f1-441d-46da-9348-53a14923619a
md"""
## Problem definition

### Conditional simulation
"""

# ╔═╡ 80e60650-1e69-11eb-376a-c1eb22aa1d0a
begin
	coords = [(50.,50.),(190.,50.),(150.,70.),(150.,190.)]
	facies = [1,0,1,1]
	
	𝒮 = georef((facies=facies,), coords)
	
	𝒟 = CartesianGrid(250, 250)
	
	𝒫₁ = SimulationProblem(𝒮, 𝒟, :facies, 3)
end

# ╔═╡ c39ecc8a-f1f4-4fe7-8a0f-a05807f8aca4
md"""
### Unconditional simulation
"""

# ╔═╡ ac17a4ac-5cab-403f-8533-3f294c905a11
𝒫₂ = SimulationProblem(𝒟, :facies => Int, 3)

# ╔═╡ 6330d732-f002-4552-b6b8-2a03e1e85ac2
md"""
## Solving the problem

Like most other multiple-point simulation solvers, image quilting is parametrized with a training image. The [GeoStatsImages.jl](https://github.com/juliohm/GeoStatsImages.jl) package provides various training images from the literature for fast experimentation in Julia. **Please give credit to the data sources if you use these images in your research**.

We load a famous image from the geostatistics literature:
"""

# ╔═╡ 149fe30a-364f-4f95-9775-21137150af45
ℐ = geostatsimage("Strebelle")

# ╔═╡ 4e13c550-1e69-11eb-1dd9-7d438c2fad76
plot(ℐ)

# ╔═╡ f56dcc7b-1d9d-42d3-9fc5-30a09c823b2f
md"""
and define our solver:
"""

# ╔═╡ 4a85e404-23b9-40b5-ad38-210139e1100a
solver = IQ(:facies => (trainimg=ℐ, tilesize=(30,30)))

# ╔═╡ d859e4a4-2e24-4336-9b3e-e6830b026faf
md"""
The solver can be used for conditional simulation:
"""

# ╔═╡ 20088240-1e69-11eb-3350-73ececd3a537
sol₁ = solve(𝒫₁, solver)

# ╔═╡ ed40a963-1570-4f8c-807f-2c3fd83e2a61
plot(sol₁, size=(700,200))

# ╔═╡ 6e308e55-8665-401a-973a-c9bb58a8c835
md"""
as well as unconditional simulation:
"""

# ╔═╡ 293a6bdb-ecb8-463b-88f7-2bfaf0918ab4
sol₂ = solve(𝒫₂, solver)

# ╔═╡ 0acd2cf0-1e69-11eb-0b4c-c1e407eb77bb
plot(sol₂, size=(700,200))

# ╔═╡ 6213b523-8c93-424d-83c5-4fe32204accf
md"""
## Remarks

- Image quilting is an efficient solver that is particularly useful when hard data is scarce. It is capable of reproducing complex texture present in training images and is very easy to tune.

- For more details, please check:
  - *Hoffimann et al. 2017.* [Stochastic simulation by image quilting of process-based geological models](http://www.sciencedirect.com/science/article/pii/S0098300417301139).
"""

# ╔═╡ Cell order:
# ╟─1ed84657-f480-4920-8719-d9123ca0d7d4
# ╟─39aebddc-1fa3-11eb-2f85-813c7473ffab
# ╠═418e9d24-1fa3-11eb-3ef9-51b170d57b9c
# ╟─2a154f98-20d8-4277-a582-02dd4a6824a5
# ╟─91dca1b8-1fa3-11eb-3896-7fd43f394a02
# ╟─3389c9f1-441d-46da-9348-53a14923619a
# ╠═80e60650-1e69-11eb-376a-c1eb22aa1d0a
# ╟─c39ecc8a-f1f4-4fe7-8a0f-a05807f8aca4
# ╠═ac17a4ac-5cab-403f-8533-3f294c905a11
# ╟─6330d732-f002-4552-b6b8-2a03e1e85ac2
# ╠═149fe30a-364f-4f95-9775-21137150af45
# ╠═4e13c550-1e69-11eb-1dd9-7d438c2fad76
# ╟─f56dcc7b-1d9d-42d3-9fc5-30a09c823b2f
# ╠═4a85e404-23b9-40b5-ad38-210139e1100a
# ╟─d859e4a4-2e24-4336-9b3e-e6830b026faf
# ╠═20088240-1e69-11eb-3350-73ececd3a537
# ╠═ed40a963-1570-4f8c-807f-2c3fd83e2a61
# ╟─6e308e55-8665-401a-973a-c9bb58a8c835
# ╠═293a6bdb-ecb8-463b-88f7-2bfaf0918ab4
# ╠═0acd2cf0-1e69-11eb-0b4c-c1e407eb77bb
# ╟─6213b523-8c93-424d-83c5-4fe32204accf
