### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 01c88650-1fa3-11eb-1156-bb632230292b
using Random; Random.seed!(2019); # make sure this tutorial is reproducible

# ╔═╡ 221618e0-6e73-11eb-047e-3b329eefa7d3
begin
	using GeoStats
	using StratiGraphics
	using SpectralGaussianSimulation
	
	solver₁ = FFTGS(:land => (variogram=GaussianVariogram(range=100.,sill=3e-2),))
	solver₂ = FFTGS(:land => (variogram=GaussianVariogram(range=100.,sill=3e-2),))
	
	procs = [GeoStatsProcess(solver) for solver in [solver₁, solver₂]];
end

# ╔═╡ cc7eef10-6e72-11eb-2cd4-abf33bd7f1ad
md"""
# Installation & Setup

Before we proceed, please install the following packages:
"""

# ╔═╡ 6c9b41d8-043f-4626-8798-4a3ae9da922b
md"""
# Creating surfaces

First, we demonstrate how Markov-Poisson sampling can be used to generate random surfaces of a 3D stratigraphic model.

We start by defining a set of 2D geostatistical processes:
"""

# ╔═╡ 34e7730e-6e73-11eb-14e4-9f7bfe169470
md"""
In this example, we define two Gaussian processes. This processes represent smooth "depositional" processes as indicated by the GaussianVariogram model. The range and sill parameters determine the "frequency" and "height" of oscillation in the horizons, respectively.

Having a set of processes defined to evolve a landscape, we need transition probabilities between these processes in order to build stratigraphy over time. Here we say that the two processes are balanced (i.e. 50% / 50% chance in the long run):
"""

# ╔═╡ 497c2820-6e73-11eb-3b09-4facec99f59a
P = [0.5 0.5
     0.5 0.5]

# ╔═╡ 519a6760-6e73-11eb-2987-fd6e15456f15
md"""
Finally, we specify a duration process that determines the (random) amount of time that each process takes before a new transition:
"""

# ╔═╡ 6029f930-6e73-11eb-3834-a355dd9550fa
ΔT = ExponentialDuration(1.0)

# ╔═╡ 6778b370-6e73-11eb-03fa-738c64b466a1
md"""
We create a geological environment with the components defined above:
"""

# ╔═╡ 7bc5b120-6e73-11eb-29c1-7de2d8ab37c4
env = Environment(procs, P, ΔT)

# ╔═╡ 8c6fbe80-6e73-11eb-388f-bb64eaf84d96
md"""
Given an initial state (i.e. flat land), and a number of epochs, we can simulate the environment to produce a geological record:
"""

# ╔═╡ ab57bb40-6e73-11eb-208e-b76a2961ad6c
begin
	nepochs = 10 # number of surfaces to simulate
	
	init = LandState(zeros(500,500)) # each surface is a 500x500 image
	
	record = simulate(env, init, nepochs);
end

# ╔═╡ bb05dc70-6e73-11eb-0e35-25fb65c02fd4
begin
	using Plots; gr(format=:png)
	
	strata = Strata(record)
	
	plot(strata, size=(600,600))
end

# ╔═╡ b419dba0-6e73-11eb-3dd7-1f5a9141066b
md"""
From the record, we can extract the surfaces that make the stratigraphic model. Two options are available for stacking the surfaces, they are :erosional (default) in which case the surfaces are eroded backward in time, and :depositional in which case the surfaces are deposited forward in time:
"""

# ╔═╡ ccec3a60-6e73-11eb-3bf6-ff8d31deaca2
md"""
We can convert the stratigraphic model into a 3D voxel model by specifying a vertical resolution:
"""

# ╔═╡ e648fed0-6e73-11eb-0907-1d7e194a38c6
begin
	model = voxelize(strata, 250) # 500x500x250 voxel model
	
	xslice = rotr90(model[25,:,:])
	yslice = rotr90(model[:,25,:])
	
	px = heatmap(xslice, title="xline")
	py = heatmap(yslice, title="yline")
	
	plot(px, py, size=(950,200), aspect_ratio=:equal, clim=(0,nepochs))
end

# ╔═╡ ef7b6d30-6e73-11eb-0427-f7d965c99f8f
md"""
# Putting it all together
We can create many such stratigraphic models by defining a simulation problem for stratigraphy:
"""

# ╔═╡ fe0288a0-6e75-11eb-1bcb-e12f7baa44cb
problem = SimulationProblem(RegularGrid{Float64}(500,500,250), :strata => Float64, 3)

# ╔═╡ 08c6ba90-6e76-11eb-059d-6f09ad18a136
solver = StratSim(:strata => (environment=env,))

# ╔═╡ 1a4221b0-6e76-11eb-26e1-6da7b3ddffe3
begin
	Random.seed!(2000)
	
	solution = solve(problem, solver)
end

# ╔═╡ 21c77c50-6e76-11eb-0621-b1a96184c856
for (i,real) in enumerate(solution[:strata])
    xslice = rotr90(real[25,:,:])
    yslice = rotr90(real[:,25,:])
    
    px = heatmap(xslice, title="realization $i (xline)", clim=(0,nepochs))
    py = heatmap(yslice, title="realization $i (yline)", clim=(0,nepochs))
    
    p = plot(px, py, aspect_ratio=:equal, size=(950,200))
    
    display(p)
end

# ╔═╡ Cell order:
# ╟─cc7eef10-6e72-11eb-2cd4-abf33bd7f1ad
# ╠═01c88650-1fa3-11eb-1156-bb632230292b
# ╟─6c9b41d8-043f-4626-8798-4a3ae9da922b
# ╠═221618e0-6e73-11eb-047e-3b329eefa7d3
# ╟─34e7730e-6e73-11eb-14e4-9f7bfe169470
# ╠═497c2820-6e73-11eb-3b09-4facec99f59a
# ╟─519a6760-6e73-11eb-2987-fd6e15456f15
# ╠═6029f930-6e73-11eb-3834-a355dd9550fa
# ╟─6778b370-6e73-11eb-03fa-738c64b466a1
# ╠═7bc5b120-6e73-11eb-29c1-7de2d8ab37c4
# ╟─8c6fbe80-6e73-11eb-388f-bb64eaf84d96
# ╠═ab57bb40-6e73-11eb-208e-b76a2961ad6c
# ╟─b419dba0-6e73-11eb-3dd7-1f5a9141066b
# ╠═bb05dc70-6e73-11eb-0e35-25fb65c02fd4
# ╟─ccec3a60-6e73-11eb-3bf6-ff8d31deaca2
# ╠═e648fed0-6e73-11eb-0907-1d7e194a38c6
# ╟─ef7b6d30-6e73-11eb-0427-f7d965c99f8f
# ╠═fe0288a0-6e75-11eb-1bcb-e12f7baa44cb
# ╠═08c6ba90-6e76-11eb-059d-6f09ad18a136
# ╠═1a4221b0-6e76-11eb-26e1-6da7b3ddffe3
# ╠═21c77c50-6e76-11eb-0621-b1a96184c856
