### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 57e916d0-183c-11eb-207a-91fa96032b25
begin
	using Distributed
	pids = [myid()]
	
	md"""
	Running on processes: $pids
	
	Use `pids = addprocs(n)` to run the notebook with `n` parallel processes.
	"""
end

# ╔═╡ 8a7dbf80-1fa1-11eb-104a-bf5115d9e6a6
@everywhere pids begin
	using Pkg; Pkg.activate(@__DIR__)
	Pkg.instantiate(); Pkg.precompile()
end

# ╔═╡ 94719c30-183c-11eb-19b7-8da27516dab1
@everywhere pids begin
	# packages used in this notebook
	using GeoStats
	
	# default plot settings
	using Plots; gr(size=(700,400))
	
	# make sure that results are reproducible
	using Random; Random.seed!(2020)
end

# ╔═╡ 873628b2-183c-11eb-2805-133d253e1e7a
md"""
# Estimation problems

In this tutorial, we explain one of the greatest features of GeoStats.jl: the ability to setup geostatistical problems indenpendently of the solution strategy.

If you are an experienced user of geostatistics or if you do research in the field, you know how hard it is to compare algorithms fairly. Often a new algorithm is proposed in the literature, and yet the task of comparing it with the state of the art is quite demanding. Even when a comparison is made by the author after a great amount of effort, it is inevitably biased.

Part of this issue is attributed to the fact that a general definition of the problem is missing. What is it that we call an "estimation problem" in geostatistics? The answer to this question is given below in the form of code.
"""

# ╔═╡ 904cd1b0-183c-11eb-350a-1f2488af5218
md"""## Problem definition

An estimation problem in geostatistics is a triplet:

1. Spatial data (i.e. data with coordinates)
2. Spatial domain (e.g. regular grid, unstructured grid)
3. Target variables (or variables to be estimated)

It doesn't involve variograms, training images, or any tuning parameter. These concepts belong to solvers. Let's make it a concrete example, and create some data. We will use the `readgeotable` utility function to read a CSV from disk and convert it into a spatial data set based on two columns `:x` and `:y` with coordinates:
"""

# ╔═╡ 99512b30-183c-11eb-1810-5b2f9c759229
𝒮 = readgeotable("data/precipitation.csv", coordnames=(:x,:y))

# ╔═╡ 9dbf33b0-183c-11eb-0c73-45f17d57e4cd
md"""Next, we define the domain in which the variables will be estimated. One of the many choices possible is the regular grid:"""

# ╔═╡ a30db120-183c-11eb-21ac-b7da79707d4e
𝒟 = CartesianGrid(100, 100)

# ╔═╡ a84dd6b0-183c-11eb-2da2-511cc3f7ec3a
md"""Notice that by default, a regular grid has zero origin and unit spacing. Alternatively, these parameters can be specified explicitly:"""

# ╔═╡ ae8617e0-183c-11eb-140c-f76c37b4d82c
CartesianGrid((100, 100), (0., 0.), (1., 1.))

# ╔═╡ b52c5d20-183c-11eb-1d37-47e6106139b0
md"""Regular grids are lightweight objects. They do not require any memory space other than the space used to save the input parameters (i.e. dimensions, origin and spacing):"""

# ╔═╡ b6ec83b0-183c-11eb-068a-0d4642430675
@allocated CartesianGrid(10^6, 10^6)

# ╔═╡ ba491c30-183c-11eb-2067-51f15b79cba6
md"""Looping over a regular grid or finding the coordinates of a given location is done analytically for maximum performance."""

# ╔═╡ c15b69b0-183c-11eb-20b4-fb650ad19121
md"""Finally, we define the estimation problem for the precipitation variable:"""

# ╔═╡ c7799332-183c-11eb-2ed2-95ee6c83a881
problem = EstimationProblem(𝒮, 𝒟, :precipitation)

# ╔═╡ cb2de8a0-183c-11eb-1941-f10a606fd760
md"""## Solving the problem

Now that the problem is unambiguously defined, we can solve it with various estimation solvers. In this tutorial, we will use the polyalgorithm Kriging solver distributed with GeoStats.jl. In Kriging, each variable of the problem is (optionally) parametrized by a mean and a variogram:
"""

# ╔═╡ cf6ddc3e-183c-11eb-1aec-89743897361c
solver = Kriging(
    :precipitation => (variogram=GaussianVariogram(range=35.),)
)

# ╔═╡ d6e7ec40-183c-11eb-0922-0305be3aad3a
md"""The line above translates to *"solve the precipitation variable using a Gaussian variogram"*. When only the variogram is specified, Ordinary Kriging is triggered. The user can specify the mean (e.g. `mean=.5`) for Simple Kriging, the polynomial degree (e.g. `degree=1`) for Universal Kriging, and the drift functions (e.g. `drifts=[x -> 1 + x[1], x -> 2x[2]]`) for External Drift Kriging. For more solver options, please consult the GeoStats.jl documentation.

The solution to the problem is easily obtained with:
"""

# ╔═╡ d8fbff30-183c-11eb-1cf3-2d2db481b2a4
solution = solve(problem, solver)

# ╔═╡ dc01ec30-183c-11eb-38d6-3b5b6b3361f1
md"""It is stored in an efficient format with all the necessary information to reconstruct the estimates spatially. Results for specific properties can be accessed:
"""

# ╔═╡ e05eb6a0-183c-11eb-3075-a3f261e0b8d5
μ, σ² = solution[:precipitation]

# ╔═╡ e7b0a530-183c-11eb-2a1e-5930ae3286ac
md"""However, very often we just want to visualize the results. In GeoStats.jl, solutions can be plotted directly in a standardized format for comparison between different solvers and parameter settings:
"""

# ╔═╡ 00691350-183d-11eb-0cec-911bccf5d86b
plot(solution)

# ╔═╡ 06b1a830-183d-11eb-19ce-ab6e9ae03bed
md"""Thanks to the integration with [Plots.jl](https://github.com/JuliaPlots/Plots.jl), many plot types are available:"""

# ╔═╡ 0c3c1d30-183d-11eb-171e-652252f6caa6
contour(solution, clabels=true)

# ╔═╡ 132510c0-183d-11eb-029b-71a0dcb0d2d1
contourf(solution, clabels=true)

# ╔═╡ 19d95fc0-183d-11eb-395c-87b08a38bbde
surface(solution, camera=(50,70))

# ╔═╡ 1b0583b0-183d-11eb-328d-359a7e550a0d
md"""
## Remarks

- The ability to work at the level of the problem definition is quite desirable. Users can switch between different solvers without having to learn new syntax. This approach also guarantees that the problem being solved is the same and that the comparison is fair.
"""

# ╔═╡ Cell order:
# ╟─57e916d0-183c-11eb-207a-91fa96032b25
# ╟─8a7dbf80-1fa1-11eb-104a-bf5115d9e6a6
# ╠═94719c30-183c-11eb-19b7-8da27516dab1
# ╟─873628b2-183c-11eb-2805-133d253e1e7a
# ╟─904cd1b0-183c-11eb-350a-1f2488af5218
# ╠═99512b30-183c-11eb-1810-5b2f9c759229
# ╟─9dbf33b0-183c-11eb-0c73-45f17d57e4cd
# ╠═a30db120-183c-11eb-21ac-b7da79707d4e
# ╟─a84dd6b0-183c-11eb-2da2-511cc3f7ec3a
# ╠═ae8617e0-183c-11eb-140c-f76c37b4d82c
# ╟─b52c5d20-183c-11eb-1d37-47e6106139b0
# ╠═b6ec83b0-183c-11eb-068a-0d4642430675
# ╟─ba491c30-183c-11eb-2067-51f15b79cba6
# ╟─c15b69b0-183c-11eb-20b4-fb650ad19121
# ╠═c7799332-183c-11eb-2ed2-95ee6c83a881
# ╟─cb2de8a0-183c-11eb-1941-f10a606fd760
# ╠═cf6ddc3e-183c-11eb-1aec-89743897361c
# ╟─d6e7ec40-183c-11eb-0922-0305be3aad3a
# ╠═d8fbff30-183c-11eb-1cf3-2d2db481b2a4
# ╟─dc01ec30-183c-11eb-38d6-3b5b6b3361f1
# ╠═e05eb6a0-183c-11eb-3075-a3f261e0b8d5
# ╟─e7b0a530-183c-11eb-2a1e-5930ae3286ac
# ╠═00691350-183d-11eb-0cec-911bccf5d86b
# ╟─06b1a830-183d-11eb-19ce-ab6e9ae03bed
# ╠═0c3c1d30-183d-11eb-171e-652252f6caa6
# ╠═132510c0-183d-11eb-029b-71a0dcb0d2d1
# ╠═19d95fc0-183d-11eb-395c-87b08a38bbde
# ╟─1b0583b0-183d-11eb-328d-359a7e550a0d
