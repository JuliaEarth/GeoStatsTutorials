### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# ╔═╡ 19a01138-bf5f-4e81-b795-66d6fdf8d7a8
begin
	using Distributed
	
	nprocs()
end

# ╔═╡ 4e30334d-6953-4624-98e6-3d6f11e8f091
begin
    using Pkg; Pkg.instantiate()
    using Random; Random.seed!(myid())
    
    using GeoStats
    using GeoStatsImages
    using ImageQuilting
    using DirectGaussianSimulation
	using Plots
end

# ╔═╡ 6a1481e3-431d-4851-b766-c1935567a5a2
md"""
# Parallel stochastic simulation

Stochastic simulation is an [embarrassingly parallel problem](https://en.wikipedia.org/wiki/Embarrassingly_parallel) in which realizations are only a function of the random seed, and can be generated indenpendently one from another. Although this is a well-known fact, there has not been reasonable effort in geostatistical software to exploit modern hardware such as HPC clusters and the cloud (e.g. AWS, MS Azure).

In GeoStats.jl, *all* stochastic simulation algorithms generate realizations in parallel by default. The package exploits Julia's built-in support for parallel execution, and works seamlessly on personal laptops with multiple cores as well as on high-performance computer clusters with multiple nodes.

In this tutorial, we demonstrate how to generate realizations with sequential Gaussian simulation in parallel. The same script can be run on a computer cluster where thousands of processes are available.
"""

# ╔═╡ 6dad5816-1a34-4a18-b1ae-0bcb49a49026
md"""
## Parallel setup

When you start Julia, it starts with a single process:
"""

# ╔═╡ aed5eddc-f6fc-442b-9566-74c85dbc9f07
md"""
In order to run simulations in parallel, the first thing we need to do is increase the number of processes in the pool *before* loading the environment. The command `addprocs` adds a given number of processes for parallel execution:
"""

# ╔═╡ c7c00ea0-1db9-11eb-0fb9-f749bbb8b8a0
if myid() == 1 && nworkers() == 1
	addprocs(topology=:master_worker, exeflags="--project=$(Base.active_project())")
end

# ╔═╡ eb6fa22c-9946-4a3b-bf52-dc71cc830e60
md"""
Notice that `addprocs()` when called without an argument, adds the number of *logical* cores available in the machine. The `topology` option specifies the type of communication allowed, in this case only the master process can communicate with workers.

> **WARNING:** If you are on Windows 8 or an older version of the operating system, you will likely experience a slow down. Please add the number of *physical* cores instead that can be found externally in computer settings or via packages such as [CpuId.jl](https://github.com/m-j-w/CpuId.jl) and [Hwloc.jl](https://github.com/JuliaParallel/Hwloc.jl).

On a HPC cluster, computing resources are generally requested via a resource manager (e.g. SLURM, PBS). In this case, the package [ClusterManagers.jl](https://github.com/JuliaParallel/ClusterManagers.jl) provides variants of the built-in `addprocs()` for adding processes to the pool effortlessly. For example, we can use `addprocs_slurm(1000)` to request 1000 processes in a SLURM job.

Now that the processes are available, we can instantiate the environment in all processes:
"""

# ╔═╡ 797bf86b-5bdf-44f1-af17-68b158bf8a84
md"""
## Problem definition

We define a simulation problem with 30 realizations of `facies` and `porosity` over a 100x100 grid:
"""

# ╔═╡ 505d6d84-2334-4abb-bc4b-a5ea262f3ab6
𝒫 = SimulationProblem(RegularGrid(100,100), (:facies=>Int, :porosity=>Float64), 30)

# ╔═╡ 1da3e28f-c045-4e5d-b194-481dcbefa5e2
md"""
## Solving the problem

We simulate `facies` with image quilting, and then for each facies type (0 or 1), we specify a different variogram model for `porosity`:
"""

# ╔═╡ c5bd3cee-1dbd-11eb-0537-e5b3d0d5d037
md"model for facies"

# ╔═╡ cbef1580-1dbd-11eb-1d10-e5c2bff79a6b
begin
	ℐ  = geostatsimage("Ellipsoids")
	TI = reshape(ℐ[:Z], size(domain(ℐ)))
	f  = ImgQuilt(:facies => (TI=TI, tilesize=(30,30)))
end

# ╔═╡ db54bc00-1dbd-11eb-3cdf-c599ac7b006d
md"model for porosity within facies 0"

# ╔═╡ e6cde610-1dbd-11eb-242f-d1322b340c13
begin	
	γ₀ = SphericalVariogram(range=20., sill=.2)
	p₀ = DirectGaussSim(:porosity => (variogram=γ₀,))
end

# ╔═╡ ef779c70-1dbd-11eb-2d73-55eb79e02b2b
md"model for porosity within facies 1"

# ╔═╡ f8f01480-1dbd-11eb-1e63-a1e571998685
begin	
	γ₁ = SphericalVariogram(range=20., distance=Ellipsoidal([10.,1.],[0.]))
	p₁ = DirectGaussSim(:porosity => (variogram=γ₁,))
end

# ╔═╡ 0145de80-1dbe-11eb-19c4-e379b371e3b5
md"combined cookie-cutter model"

# ╔═╡ 093daa00-1dbe-11eb-2192-adbb2001347d
ℳ = CookieCutter(f, Dict(0 => p₀, 1 => p₁))

# ╔═╡ b791b275-bb23-4c4a-a5c5-f822c3da1a94
sol = solve(𝒫, ℳ)

# ╔═╡ 9ef8cb73-4c03-4c66-9776-63cd115dc2f1
plot(sol, c=:cividis, size=(900,600))

# ╔═╡ 55aaba80-7174-4035-9921-24271f993166
md"""
## Remarks

- Users can utilize hardware resources such as HPC clusters and the cloud to speed up simulation jobs. In theory, the speedup is upper bounded by the number of computing nodes/cores/processes available. In practice, this translates into generating multiple realizations within the time required for one.
- On a personal laptop or desktop, serial (non-parallel) simulation can often outperform parallel simulation.
"""

# ╔═╡ Cell order:
# ╟─6a1481e3-431d-4851-b766-c1935567a5a2
# ╟─6dad5816-1a34-4a18-b1ae-0bcb49a49026
# ╠═19a01138-bf5f-4e81-b795-66d6fdf8d7a8
# ╟─aed5eddc-f6fc-442b-9566-74c85dbc9f07
# ╠═c7c00ea0-1db9-11eb-0fb9-f749bbb8b8a0
# ╟─eb6fa22c-9946-4a3b-bf52-dc71cc830e60
# ╠═4e30334d-6953-4624-98e6-3d6f11e8f091
# ╟─797bf86b-5bdf-44f1-af17-68b158bf8a84
# ╠═505d6d84-2334-4abb-bc4b-a5ea262f3ab6
# ╟─1da3e28f-c045-4e5d-b194-481dcbefa5e2
# ╟─c5bd3cee-1dbd-11eb-0537-e5b3d0d5d037
# ╠═cbef1580-1dbd-11eb-1d10-e5c2bff79a6b
# ╟─db54bc00-1dbd-11eb-3cdf-c599ac7b006d
# ╠═e6cde610-1dbd-11eb-242f-d1322b340c13
# ╟─ef779c70-1dbd-11eb-2d73-55eb79e02b2b
# ╠═f8f01480-1dbd-11eb-1e63-a1e571998685
# ╟─0145de80-1dbe-11eb-19c4-e379b371e3b5
# ╠═093daa00-1dbe-11eb-2192-adbb2001347d
# ╠═b791b275-bb23-4c4a-a5c5-f822c3da1a94
# ╠═9ef8cb73-4c03-4c66-9776-63cd115dc2f1
# ╟─55aaba80-7174-4035-9921-24271f993166
