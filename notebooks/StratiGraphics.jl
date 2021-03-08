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
	using StratiGraphics
	
	# default plot settings
	using Plots; gr(size=(700,400), format=:png)
	
	# make sure that results are reproducible
	using Random; Random.seed!(2021)
end

# ╔═╡ b3fa8609-b27e-4ad2-9836-b7841d6a8db0
md"""
## Creating surfaces

First, we demonstrate how Markov-Poisson sampling can be used to generate random surfaces of a 3D stratigraphic model.

We start by defining a set of 2D geostatistical processes:
"""

# ╔═╡ 61daf6c2-71a3-11eb-3810-ff67b544d72b
begin
	solver₁ = FFTGS(:land => (variogram=GaussianVariogram(range=100.,sill=3e-2),))
	solver₂ = FFTGS(:land => (variogram=GaussianVariogram(range=100.,sill=3e-2),))
	
	procs = [GeoStatsProcess(solver) for solver in [solver₁, solver₂]]
end;

# ╔═╡ 700cff42-71a3-11eb-0bd0-43662099751a
md"""
In this example, we define two Gaussian processes. This processes represent smooth "depositional" processes as indicated by the GaussianVariogram model. The range and sill parameters determine the "frequency" and "height" of oscillation in the horizons, respectively.

Having a set of processes defined to evolve a landscape, we need transition probabilities between these processes in order to build stratigraphy over time. Here we say that the two processes are balanced (i.e. 50% / 50% chance in the long run):
"""

# ╔═╡ 7cfdd530-71a3-11eb-212f-5f05c6e9a32b
P = [0.5 0.5
     0.5 0.5]

# ╔═╡ 8195379e-71a3-11eb-007f-556d235edf2e
md"""
Finally, we specify a duration process that determines the (random) amount of time that each process takes before a new transition:
"""

# ╔═╡ 87c0f5b0-71a3-11eb-2bd9-5b17609aa37e
ΔT = ExponentialDuration(1.0)

# ╔═╡ 8caa48b0-71a3-11eb-13b6-ad167093b2b5
md"""
We create a geological environment with the components defined above:
"""

# ╔═╡ 943ce9c0-71a3-11eb-1f36-c1575df9d64d
env = Environment(procs, P, ΔT);

# ╔═╡ a39ee6c0-71a3-11eb-37a6-6b10b691816c
md"""
Given an initial state (i.e. flat land), and a number of epochs, we can simulate the environment to produce a geological record:
"""

# ╔═╡ 98ad8a52-71a3-11eb-11f7-697f1f42f72e
begin
	nepochs = 10 # number of surfaces to simulate
	
	init = LandState(zeros(500,500)) # each surface is a 500x500 image
	
	record = simulate(env, init, nepochs)
end;

# ╔═╡ ac62b4d0-71a3-11eb-3f18-45756be15fba
md"""
From the record, we can extract the surfaces that make the stratigraphic model. Two options are available for stacking the surfaces, they are :erosional (default) in which case the surfaces are eroded backward in time, and :depositional in which case the surfaces are deposited forward in time:
"""

# ╔═╡ b43481c0-71a3-11eb-358e-6314ea693e5e
begin
	strata = Strata(record)
	
	plot(strata, size=(600,600))
end

# ╔═╡ c1813350-71a3-11eb-1d48-ab432f7a2f4b
md"""
We can convert the stratigraphic model into a 3D voxel model by specifying a vertical resolution:
"""

# ╔═╡ ca98edc0-71a3-11eb-25d8-bde416e1b51a
begin
	model = voxelize(strata, 250) # 500x500x250 voxel model
	
	xslice = rotr90(model[25,:,:])
	yslice = rotr90(model[:,25,:])
	
	px = heatmap(xslice, title="xline")
	py = heatmap(yslice, title="yline")
	
	plot(px, py, size=(950,200), aspect_ratio=:equal, clim=(0,nepochs))
end

# ╔═╡ d06dd940-71a3-11eb-25ba-5d1cf54e15e2
md"""
## Putting it all together

We can create many such stratigraphic models by defining a simulation problem for stratigraphy:
"""

# ╔═╡ e00676a0-71a3-11eb-1355-4b949a460052
problem = SimulationProblem(CartesianGrid(500,500,250), :strata => Float64, 3)

# ╔═╡ e79e47d0-71a3-11eb-1fa7-6da8d2f1a167
md"""
The StratSim solver is compliant with the GeoStats.jl API:
"""

# ╔═╡ ee764b70-71a3-11eb-2d35-43056ba7f528
solver = StratSim(:strata => (environment=env,))

# ╔═╡ fad1e640-71a3-11eb-1518-0d2eb7b95875
solution = solve(problem, solver)

# ╔═╡ 05b8944e-71a4-11eb-269a-cb9412e8e223
md"""
The solution contains 3 realizations of stratrigraphy, which can be visualized with:
"""

# ╔═╡ fb8cea20-71a9-11eb-2a4c-4beb0a2f20d9
begin
	plts = []
	for (i,real) in enumerate(solution)
		# reshape flat realization to cube
		r = reshape(real[:strata], 500, 500, 250)
		
		# take vertical slices
	    xslice = rotr90(r[25,:,:])
	    yslice = rotr90(r[:,25,:])
	    
	    px = heatmap(xslice, title="realization $i (xline)", clim=(0,nepochs))
	    py = heatmap(yslice, title="realization $i (yline)", clim=(0,nepochs))
	    
	    plt = plot(px, py, aspect_ratio=:equal, size=(950,200))
	    
	    push!(plts, plt)
	end
	plot(plts..., layout=(3,1), size=(750, 700))
end

# ╔═╡ Cell order:
# ╟─2b9ebaaa-1f9c-11eb-39ff-a5ee1ecb90ee
# ╟─a47c3e6d-9687-4bb9-b71f-3ee242ede575
# ╠═fc7b9a21-34dd-40ba-9d57-9f785904e307
# ╟─b3fa8609-b27e-4ad2-9836-b7841d6a8db0
# ╠═61daf6c2-71a3-11eb-3810-ff67b544d72b
# ╟─700cff42-71a3-11eb-0bd0-43662099751a
# ╠═7cfdd530-71a3-11eb-212f-5f05c6e9a32b
# ╟─8195379e-71a3-11eb-007f-556d235edf2e
# ╠═87c0f5b0-71a3-11eb-2bd9-5b17609aa37e
# ╟─8caa48b0-71a3-11eb-13b6-ad167093b2b5
# ╠═943ce9c0-71a3-11eb-1f36-c1575df9d64d
# ╟─a39ee6c0-71a3-11eb-37a6-6b10b691816c
# ╠═98ad8a52-71a3-11eb-11f7-697f1f42f72e
# ╟─ac62b4d0-71a3-11eb-3f18-45756be15fba
# ╠═b43481c0-71a3-11eb-358e-6314ea693e5e
# ╟─c1813350-71a3-11eb-1d48-ab432f7a2f4b
# ╠═ca98edc0-71a3-11eb-25d8-bde416e1b51a
# ╟─d06dd940-71a3-11eb-25ba-5d1cf54e15e2
# ╠═e00676a0-71a3-11eb-1355-4b949a460052
# ╟─e79e47d0-71a3-11eb-1fa7-6da8d2f1a167
# ╠═ee764b70-71a3-11eb-2d35-43056ba7f528
# ╠═fad1e640-71a3-11eb-1518-0d2eb7b95875
# ╟─05b8944e-71a4-11eb-269a-cb9412e8e223
# ╠═fb8cea20-71a9-11eb-2a4c-4beb0a2f20d9
