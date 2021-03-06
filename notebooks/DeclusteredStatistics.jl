### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ bc351d82-394f-45b6-a78e-88b2e0ff8eaf
begin
	using Distributed
	pids = [myid()]
	
	md"""
	Running on processes: $pids
	
	Use `pids = addprocs(n)` to run the notebook with `n` parallel processes.
	"""
end

# ╔═╡ 7f452f9e-1fa2-11eb-2456-85df17947152
@everywhere pids begin
	using Pkg; Pkg.activate(@__DIR__)
	Pkg.instantiate(); Pkg.precompile()
end

# ╔═╡ 8ced4438-1fa2-11eb-0f34-5d4b6ae2e1d3
@everywhere pids begin
	# packages used in this notebook
	using GeoStats
	using GeoStatsImages
	
	# default plot settings
	using Plots; gr(size=(700,400), c=:cividis)
	
	# make sure that results are reproducible
	using Random; Random.seed!(2020)
end

# ╔═╡ 58cbb430-74bc-45b0-b371-46d8db6ea2bf
md"""
# How much Gold there is left?

We will pretend that we are a mining company interested in extracting a mineral resource such as Gold (Au) from a mine $\mathcal{D}$. Based on the fact that we are trying to maximize profit, we end up collecting samples $\mathcal{S} = \{Au(x_1),Au(x_2),\ldots,Au(x_N)\}$ at nearby locations $x_1,x_2,\ldots,x_N$ where we believe we will find the mineral we are after.

Our data acquisition process suffers from *sampling bias*: most of our samples reflect high concentrations of Gold.

**Problem statement**:

Estimate the remaining amount of Gold in the mine $\mathcal{D}$ from the spatial samples $\mathcal{S}$.

*Can we do any better than multiplying the sample average by the volume of the mine?*

$$\underbrace{\left(\frac{1}{N}\sum_{i=1}^{N} Au(x_i)\right)}_{\mu_\mathcal{S}:\text{ sample average}} \times \mathcal{V}(\mathcal{D})$$
"""

# ╔═╡ fc60b540-1a0f-11eb-3bd0-23c1ea40c42c
begin
	# reference image representing a mine
	ℐ = geostatsimage("WalkerLakeTruth")
	
	# rename variable Z to Au
	𝒟 = georef((Au=ℐ[:Z],), domain(ℐ))
	
	# sample from spatial data
	𝒮 = sample(𝒟, 50, 𝒟[:Au], replace=false)
	
	plot(plot(𝒟), plot(𝒮))
end

# ╔═╡ a779b5df-19d5-4de1-90c1-499fef8e0a99
md"""
The mean value of Gold in the mine is:
"""

# ╔═╡ c5a1c13e-cd05-41bc-a20a-ce1d2912f3b4
μ𝒟 = mean(𝒟[:Au])

# ╔═╡ 06b89221-c543-45be-baf9-513ccc7a8762
md"""
whereas the sample average is much higher:
"""

# ╔═╡ 1f984003-b063-4400-bd26-9d5888bd923a
μ𝒮 = mean(𝒮[:Au])

# ╔═╡ c6d99824-249b-4517-982d-ea951a826dbb
md"""
## Spatial declustering

Notice that besides suffering from sampling bias, our sampling process leads to samples that are *clustered in space*. To quantify this observation, let's partition $\mathcal{S}$ into blocks $\mathcal{B}_1,\mathcal{B}_2,\ldots,\mathcal{B}_M$ of size $b^2$ and count the number of samples that share a block:
"""

# ╔═╡ 41cfd2ee-1a10-11eb-2509-911a08eb0e16
begin
	ℬ = partition(𝒮, BlockPartition(50.,50.))
	p₁ = plot(ℬ, colorbar=false, xlabel="x", ylabel="y")
	p₂ = bar(nelements.(ℬ), xlabel="block", ylabel="counts", legend=false)
	plot(p₁, p₂)
end

# ╔═╡ 0857c193-0866-4c79-be0d-74750e0cc5ef
md"""
Samples that are close to each other are redundant, and shouldn't receive the same "importance" in the mean estimate compared with isolated samples. We can use the block counts $|\mathcal{B}_j|$ to assign importance weights $w_b(x_i) = \frac{1}{|\mathcal{B}_j|}$ to the samples $Au(x_i)$ based on their locations $x_i \in \mathcal{B}_j$:
"""

# ╔═╡ 16a8cebb-0427-4cab-8697-be5b8274812e
begin
	𝒲 = weight(𝒮, BlockWeighting(50.,50.))
	plot(𝒲, c=:Oranges)
end

# ╔═╡ f7fd71ae-2995-4176-984b-9524cad8c199
md"""
These weights $w_b(x_i)$, which are a function of the block size $b^2$, can be used in a weighted average

$$\mu_\mathcal{B} = \frac{\sum_{i=1}^N w_b(x_i) Au(x_i)}{\sum_{i=1}^N w_b(x_i)}$$

that generalizes the sample average $\mu_\mathcal{S} = \lim_{b\to 0} \mu_\mathcal{B}$.

We can plot the weighted average for increasing block sizes to notice that the sample average (i.e. uniform weights) is recovered when the block size is too small (each sample is its own block), or when the block size is too large (all samples are in a single block):
"""

# ╔═╡ bb3b8830-1a0d-11eb-3578-3fb0b1f8ffe0
begin
	bs = range(1, stop=120, length=100)
	ms = [mean(𝒮, :Au, b) for b in bs]
	
	plot(xlabel="block size", ylabel="mean estimate", legend=:bottomright)
	plot!(bs, ms, c=:green, label="weighted average")
	hline!([μ𝒮], c=:red, ls=:dash, label="sample average")
	hline!([μ𝒟], c=:black, ls=:dash, label="true average")
end

# ╔═╡ 5d87af07-ffba-4e02-ac1b-b168b5d0b304
md"""
In case the block size is ommited, GeoStats.jl uses a heuristic to select a "reasonable" block size for the given spatial configuration:
"""

# ╔═╡ a1c3e8c0-1a0d-11eb-2163-659f72901110
μℬ = mean(𝒮, :Au)

# ╔═╡ a52f0030-1a0d-11eb-2412-c90fa06d295a
md"""
| Sample average | Weighted average | True average |
|:--------------:|:----------------:|:------------:|
| $μ𝒮            | $μℬ              | $μ𝒟          |
"""

# ╔═╡ a9bca08e-2a23-49b4-9f99-2267ff827e27
md"""
We can compare the difference, in volume of Gold, between the two statistics:
"""

# ╔═╡ 4644c09e-1a0d-11eb-3db4-afbc0bbaad23
𝒱 = measure(boundingbox(𝒮))

# ╔═╡ 096d397d-370c-4d6a-b01b-45cac8981adc
(μ𝒮 - μℬ) * 𝒱

# ╔═╡ d9daf7aa-271b-4833-a15f-a3b71dd9da4a
md"""
### Declustered statistics

The idea of assigning importance weights to samples via spatial declustering is general, and holds for any statistic of interest. Hence, the term *declustered statistics*. To give another example, we can obtain better estimates of any quantile of the Gold distribution by considering the coordinates of the samples:

#### Non-spatial quantile
"""

# ╔═╡ d10fd364-dfe2-4011-a96c-d986a3b1e8a1
quantile(𝒮[:Au], [0.25,0.50,0.75])

# ╔═╡ d52c29ee-535d-4b6f-a277-00e39c1a44cd
md"""
#### Spatial quantile
"""

# ╔═╡ e87cbc5e-b1f3-4ee0-90b2-6e83ff2ac862
quantile(𝒮, :Au, [0.25,0.50,0.75])

# ╔═╡ 7084342c-622b-464c-9c63-c9434bad3b2e
md"""
## Remarks

- Spatial samples can be weighted based on their coordinates to improve volumetric estimates
- Spatial declustering is particularly useful in the presence of sampling bias and correlation
- GeoStats.jl changes the semantics of statistics such as `mean`, `var` and `quantile`
"""

# ╔═╡ Cell order:
# ╟─bc351d82-394f-45b6-a78e-88b2e0ff8eaf
# ╟─7f452f9e-1fa2-11eb-2456-85df17947152
# ╠═8ced4438-1fa2-11eb-0f34-5d4b6ae2e1d3
# ╟─58cbb430-74bc-45b0-b371-46d8db6ea2bf
# ╠═fc60b540-1a0f-11eb-3bd0-23c1ea40c42c
# ╟─a779b5df-19d5-4de1-90c1-499fef8e0a99
# ╠═c5a1c13e-cd05-41bc-a20a-ce1d2912f3b4
# ╟─06b89221-c543-45be-baf9-513ccc7a8762
# ╠═1f984003-b063-4400-bd26-9d5888bd923a
# ╟─c6d99824-249b-4517-982d-ea951a826dbb
# ╠═41cfd2ee-1a10-11eb-2509-911a08eb0e16
# ╟─0857c193-0866-4c79-be0d-74750e0cc5ef
# ╠═16a8cebb-0427-4cab-8697-be5b8274812e
# ╟─f7fd71ae-2995-4176-984b-9524cad8c199
# ╠═bb3b8830-1a0d-11eb-3578-3fb0b1f8ffe0
# ╟─5d87af07-ffba-4e02-ac1b-b168b5d0b304
# ╠═a1c3e8c0-1a0d-11eb-2163-659f72901110
# ╟─a52f0030-1a0d-11eb-2412-c90fa06d295a
# ╟─a9bca08e-2a23-49b4-9f99-2267ff827e27
# ╠═4644c09e-1a0d-11eb-3db4-afbc0bbaad23
# ╠═096d397d-370c-4d6a-b01b-45cac8981adc
# ╟─d9daf7aa-271b-4833-a15f-a3b71dd9da4a
# ╠═d10fd364-dfe2-4011-a96c-d986a3b1e8a1
# ╟─d52c29ee-535d-4b6f-a277-00e39c1a44cd
# ╠═e87cbc5e-b1f3-4ee0-90b2-6e83ff2ac862
# ╟─7084342c-622b-464c-9c63-c9434bad3b2e
