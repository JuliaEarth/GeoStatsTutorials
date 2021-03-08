### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ a83fe188-1fa4-11eb-1483-2f11ed5794cb
begin
	using Distributed
	pids = [myid()]
	
	md"""
	Running on processes: $pids
	
	Use `pids = addprocs(n)` to run the notebook with `n` parallel processes.
	"""
end

# ╔═╡ 0daa81e3-14be-48b2-a154-470e931468d9
@everywhere pids begin
	using Pkg; Pkg.activate(@__DIR__)
	Pkg.instantiate(); Pkg.precompile()
end

# ╔═╡ b293b968-1fa4-11eb-2011-934f29cff5d8
@everywhere pids begin
	# packages used in this notebook
	using GeoStats
	using ImageQuilting
	using GeoStatsImages
	
	# default plot settings
	using Plots; gr(size=(700,400), c=:cividis)
	
	# make sure that results are reproducible
	using Random; Random.seed!(1234)
end

# ╔═╡ 3ca175d4-e123-4491-aca1-51b81a6813e5
md"""
# Cookie-cutter

In this tutorial, we illustrate a simplistic, yet useful procedure for generating geostatistical realizations of the subsurface taking into account various lithology. The procedure known in the literature as "cookie-cutter" [Begg 1992](https://www.onepetro.org/conference-paper/SPE-24698-MS) consists of first simulating a categorical variable (a.k.a. the facies) and then populating each simulated category with a separate simulation solver.

Regardless of the artifacts that the procedure may create, it is very useful for building complex geological models that show strong contrasts (e.g. highly permeable channels within less permeable media). These contrasts are relevant for flow simulation studies and other applications.
"""

# ╔═╡ f42e0618-68d7-405a-bcba-a5ce4772ee99
md"""
## Problem definition

We define an unconditional simulation problem for simplicity, but the same steps apply in the presence of data.
"""

# ╔═╡ 7900db30-1c76-11eb-3ce7-6f3528ef7446
𝒟 = CartesianGrid(100, 100)

# ╔═╡ 7c101700-1c76-11eb-2b3d-a31150a46ee9
vars = (:facies => Int, :porosity => Float64)

# ╔═╡ 80196a40-1c76-11eb-0bfb-a7eb71bc08c6
problem = SimulationProblem(𝒟, vars, 3)

# ╔═╡ 7f1a7b30-d5b4-4630-bdaa-cfd9929a1912
md"""
In this problem, we will simulate the `facies` variable and use its realizations to guide the simulation of `porosity`.
"""

# ╔═╡ fa5b5983-c6fc-4004-b720-10a3ed6bbfb6
md"""
## Solving the problem

We define the facies simulation solver based on a training image that has two categories:
"""

# ╔═╡ 29bab540-1c77-11eb-399a-4943758c72ac
begin
	ℐ = geostatsimage("Strebelle")
	plot(ℐ)
end

# ╔═╡ 930d162e-5cca-463f-bc00-e39740b4380e
md"""
Image quilting is a good default for training-image-based simulation:
"""

# ╔═╡ ae30e250-1c76-11eb-25e2-6149bdd11c4e
fsolver = IQ(:facies => (trainimg=ℐ, tilesize=(30,30)))

# ╔═╡ 727a67da-4424-439b-b86c-3ab1a0de9f99
md"""
Because there are two categorical values (0 and 1), we define two solvers:
"""

# ╔═╡ 99333240-1c76-11eb-0920-cb7761577958
psolver₀ = LUGS(
	:porosity => (variogram=SphericalVariogram(range=20., sill=.2),)
)

# ╔═╡ 9c917870-1c76-11eb-35ea-29dd220e13de
psolver₁ = LUGS(
	:porosity => (
		variogram=SphericalVariogram(
			range=20., distance=aniso2distance([10.,1.],[0.])
		),
	)
)

# ╔═╡ 7a8cc266-b367-49ff-93fd-aaa709ea850e
md"""
Finally, we create the cookie-cutter procedure by specifying the master (a.k.a. facies) solver, and a solver for each categorical value:
"""

# ╔═╡ c4df62a6-dad3-4b3e-8a73-a8af28cf0928
solver = CookieCutter(fsolver, Dict(0 => psolver₀, 1 => psolver₁))

# ╔═╡ 1ab667e1-e306-4046-af68-2c02e8ebc463
md"""
GeoStats.jl will generate the realizations in parallel as usual for each solver, and then it will combine the results into a single solution:
"""

# ╔═╡ 5943fd0b-eb3c-4242-816e-fa1f1173f8e5
solution = solve(problem, solver)

# ╔═╡ c5fcb38a-7826-48c1-9dfc-9fefaa5393c1
md"""
We can plot each variable to confirm that the procedure was effective:
"""

# ╔═╡ 60797bf5-6520-4194-b648-c4e9b62db5e6
plot(solution, size=(700,500))

# ╔═╡ f8f53d25-04ef-4606-9ba9-8254c240ad0b
md"""
## Remarks

- The cookie-cutter procedure in GeoStats.jl gives users the ability to create very complex combinations of patterns that are difficult to generate otherwise with a single simulation algorithm.

- Any simulation solver adhering to the interface proposed in the framework can be used as a building block for cookie-cutter. This feature opens the door to a whole new set of models, which cannot be generated in other software.
"""

# ╔═╡ Cell order:
# ╟─a83fe188-1fa4-11eb-1483-2f11ed5794cb
# ╟─0daa81e3-14be-48b2-a154-470e931468d9
# ╠═b293b968-1fa4-11eb-2011-934f29cff5d8
# ╟─3ca175d4-e123-4491-aca1-51b81a6813e5
# ╟─f42e0618-68d7-405a-bcba-a5ce4772ee99
# ╠═7900db30-1c76-11eb-3ce7-6f3528ef7446
# ╠═7c101700-1c76-11eb-2b3d-a31150a46ee9
# ╠═80196a40-1c76-11eb-0bfb-a7eb71bc08c6
# ╟─7f1a7b30-d5b4-4630-bdaa-cfd9929a1912
# ╟─fa5b5983-c6fc-4004-b720-10a3ed6bbfb6
# ╠═29bab540-1c77-11eb-399a-4943758c72ac
# ╟─930d162e-5cca-463f-bc00-e39740b4380e
# ╠═ae30e250-1c76-11eb-25e2-6149bdd11c4e
# ╟─727a67da-4424-439b-b86c-3ab1a0de9f99
# ╠═99333240-1c76-11eb-0920-cb7761577958
# ╠═9c917870-1c76-11eb-35ea-29dd220e13de
# ╟─7a8cc266-b367-49ff-93fd-aaa709ea850e
# ╠═c4df62a6-dad3-4b3e-8a73-a8af28cf0928
# ╟─1ab667e1-e306-4046-af68-2c02e8ebc463
# ╠═5943fd0b-eb3c-4242-816e-fa1f1173f8e5
# ╟─c5fcb38a-7826-48c1-9dfc-9fefaa5393c1
# ╠═60797bf5-6520-4194-b648-c4e9b62db5e6
# ╟─f8f53d25-04ef-4606-9ba9-8254c240ad0b
