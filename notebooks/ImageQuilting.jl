### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# ╔═╡ 1ed84657-f480-4920-8719-d9123ca0d7d4
using Pkg; Pkg.instantiate(); Pkg.precompile()

# ╔═╡ 6d3a3573-a307-44a1-a4ab-b975d115cbcc
using Random; Random.seed!(2000);

# ╔═╡ 293b3ec0-1e69-11eb-3326-9d65d2fe99e3
begin
	using GeoStats
	using ImageQuilting
	using Plots; gr(c=:cividis)
	using GeoStatsImages
end

# ╔═╡ 2a154f98-20d8-4277-a582-02dd4a6824a5
md"""
# Image quilting

In this tutorial, we demonstrate the use of the image quilting solver, which is the fastest multiple-point simulation solver in the literature. For more information about the algorithm and its performance, please watch [this video](https://www.youtube.com/watch?v=YJs7jl_Y9yM).
"""

# ╔═╡ 3389c9f1-441d-46da-9348-53a14923619a
md"""
## Problem definition

### Conditional simulation
"""

# ╔═╡ 82554dc0-1e69-11eb-16f4-6f30549fc0fe
md"create some artificial data"

# ╔═╡ 80e60650-1e69-11eb-376a-c1eb22aa1d0a
X = [50. 190. 150. 150.
	     50. 50.  70.  190.]

# ╔═╡ 7c090f62-1e69-11eb-156a-c7b8a7852569
z = [1,0,1,1]

# ╔═╡ 79b2c352-1e69-11eb-0d05-4b479a1ee254
𝒮 = georef((facies=z,), X)

# ╔═╡ 742abf50-1e69-11eb-001f-11f8ac49ddac
𝒟 = RegularGrid(250, 250)

# ╔═╡ 71494a40-1e69-11eb-0f76-bff1a4f3a5fa
𝒫₁ = SimulationProblem(𝒮, 𝒟, :facies, 3)

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

# ╔═╡ 4e13c550-1e69-11eb-1dd9-7d438c2fad76
ℐ = geostatsimage("Strebelle")

# ╔═╡ 502455d0-1e69-11eb-10d8-cf92a45ae646
plot(ℐ)

# ╔═╡ f56dcc7b-1d9d-42d3-9fc5-30a09c823b2f
md"""
and define our solver:
"""

# ╔═╡ 4a85e404-23b9-40b5-ad38-210139e1100a
begin
	# convert spatial variable to Julia array
	TI = reshape(ℐ[:facies], size(domain(ℐ)))
	
	solver = ImgQuilt(:facies => (TI=TI, tilesize=(30,30)))
end

# ╔═╡ d859e4a4-2e24-4336-9b3e-e6830b026faf
md"""
The solver can be used for conditional simulation:
"""

# ╔═╡ 20088240-1e69-11eb-3350-73ececd3a537
sol₁ = solve(𝒫₁, solver)

# ╔═╡ ed40a963-1570-4f8c-807f-2c3fd83e2a61
plot(sol₁, size=(900,300))

# ╔═╡ 6e308e55-8665-401a-973a-c9bb58a8c835
md"""
as well as unconditional simulation:
"""

# ╔═╡ 293a6bdb-ecb8-463b-88f7-2bfaf0918ab4
sol₂ = solve(𝒫₂, solver)

# ╔═╡ 0acd2cf0-1e69-11eb-0b4c-c1e407eb77bb
plot(sol₂, size=(900,300))

# ╔═╡ 6213b523-8c93-424d-83c5-4fe32204accf
md"""
## Conclusions

Image quilting is an efficient solver that is particularly useful when hard data is scarce. It is capable of reproducing complex texture present in training images and is very easy to tune. For more details, please refer to [Hoffimann 2017](http://www.sciencedirect.com/science/article/pii/S0098300417301139).
"""

# ╔═╡ Cell order:
# ╠═1ed84657-f480-4920-8719-d9123ca0d7d4
# ╠═6d3a3573-a307-44a1-a4ab-b975d115cbcc
# ╠═293b3ec0-1e69-11eb-3326-9d65d2fe99e3
# ╟─2a154f98-20d8-4277-a582-02dd4a6824a5
# ╟─3389c9f1-441d-46da-9348-53a14923619a
# ╟─82554dc0-1e69-11eb-16f4-6f30549fc0fe
# ╠═80e60650-1e69-11eb-376a-c1eb22aa1d0a
# ╠═7c090f62-1e69-11eb-156a-c7b8a7852569
# ╠═79b2c352-1e69-11eb-0d05-4b479a1ee254
# ╠═742abf50-1e69-11eb-001f-11f8ac49ddac
# ╠═71494a40-1e69-11eb-0f76-bff1a4f3a5fa
# ╟─c39ecc8a-f1f4-4fe7-8a0f-a05807f8aca4
# ╠═ac17a4ac-5cab-403f-8533-3f294c905a11
# ╟─6330d732-f002-4552-b6b8-2a03e1e85ac2
# ╠═4e13c550-1e69-11eb-1dd9-7d438c2fad76
# ╠═502455d0-1e69-11eb-10d8-cf92a45ae646
# ╟─f56dcc7b-1d9d-42d3-9fc5-30a09c823b2f
# ╠═4a85e404-23b9-40b5-ad38-210139e1100a
# ╟─d859e4a4-2e24-4336-9b3e-e6830b026faf
# ╠═20088240-1e69-11eb-3350-73ececd3a537
# ╠═ed40a963-1570-4f8c-807f-2c3fd83e2a61
# ╟─6e308e55-8665-401a-973a-c9bb58a8c835
# ╠═293a6bdb-ecb8-463b-88f7-2bfaf0918ab4
# ╠═0acd2cf0-1e69-11eb-0b4c-c1e407eb77bb
# ╟─6213b523-8c93-424d-83c5-4fe32204accf
