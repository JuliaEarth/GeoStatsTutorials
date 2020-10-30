### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 07589519-3c07-4a44-a3a2-f312213687c1
using Pkg; Pkg.instantiate(); Pkg.precompile()

# ╔═╡ 4a66263a-7fba-4294-89d3-55a5cc82bd27
using Random; Random.seed!(2020);

# ╔═╡ 15ed2ae0-1a81-11eb-13c3-01781cc979f6
begin
	using GeoStats
	using Plots; gr(c=:cividis, ms=2)
	using GeoStatsImages
end

# ╔═╡ 10ac0071-273b-4c9e-9b97-c94071dc826f
md"""
# What is the average grain size?

We have a Micro-CT image of a rock sample such as the *Ketton* image provided by the [PERM Research Group at Imperial College London](http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/micro-ct-images-and-networks), and wish to estimate the average grain size within this sample. The average grain size can be an important parameter for subsurface characterization, but assigning a single value for the sample objectively is not always trivial.
"""

# ╔═╡ 644bcfa0-1a83-11eb-2a92-03074995e820
html"""
<img src="https://raw.githubusercontent.com/JuliaEarth/GeoStatsImages.jl/master/src/data/Ketton.png" width="30%">
<div style="text-align: center"> <b>Image credit:</b> <i>PERM at Imperial College London</i> </div>
"""

# ╔═╡ 54a3f000-1a83-11eb-2d1c-55fbe1506269
md"""
**Problem statement**:

Estimate the average grain size (e.g. average grain radius) in the rock sample.

*Can we do any better than visual estimation?*
"""

# ╔═╡ 49634762-156e-4f01-a4b9-a4b2e74597fa
𝒦 = geostatsimage("Ketton")

# ╔═╡ d5870dd6-a294-4cda-af18-88f51abbeaa5
md"""
## Two-point statistics

We will solve this problem by estimating the probability of two locations ``\newcommand{u}{\mathbf{u}}\newcommand{v}{\mathbf{v}}\newcommand{\R}{\mathbb{R}} u,v\in\R^3`` to be on a grain, and by understanding how this probability changes as the distance ``h=||u-v||`` between these locations is made larger.

Consider two Bernoulli random variables ``Z_u`` and ``Z_v`` that indicate *grain* with probability ``p`` and *non-grain* with probability ``1-p`` at locations ``u`` and ``v``, respectively. The variogram between these locations is defined as:

``
\newcommand{\E}{\mathbb{E}}
\gamma(h) = \frac{1}{2} \E\left[(Z_u - Z_v)^2\right] = \frac{1}{2}\left(\E\left[Z_u^2\right] + \E\left[Z_v^2\right] - 2\E\left[Z_u Z_v\right]\right)
``

We can derive the expectations:

- ``\E\left[Z_u^2\right] = \E[Z_u Z_u] = Pr(Z_u=1, Z_u=1) = Pr(Z_u=1)``
- ``\E[Z_u Z_v] = Pr(Z_u=1, Z_v=1) = Pr(Z_u=1)Pr(Z_v=1 | Z_u=1)``
- ``\E\left[Z_v^2\right] = Pr(Z_v=1)`` (exactly like above)

and use the fact that ``Pr(Z_u=1) = Pr(Z_v=1) = p`` to rewrite the variogram as:

``
\gamma(h) = p\cdot\big(1 - Pr(Z_v=1 | Z_u=1)\big) = p\cdot Pr(Z_v=0|Z_u=1)
``

or alternatively:

``
Pr(Z_v=0 | Z_u=1) = \frac{\gamma(h)}{p}
``

This last equation shows that we can calculate the probability that a location ``v`` is *non-grain* given that a location ``u`` is *grain* (left-hand-side) provided that we can efficiently estimate the variogram ``\gamma(h)`` and the marginal probability ``p`` from spatial data (right-hand-side).
"""

# ╔═╡ cf926398-2cf7-409f-b5c1-c0a25955098a
md"""
### Partition variograms

Efficient parallel algorithms for variogram estimation do exist [Hoffimann & Zadrozny 2019](https://www.sciencedirect.com/science/article/pii/S0098300419302936). Here, we consider 1M samples from the Ketton rock to demonstrate GeoStats.jl's computational performance:
"""

# ╔═╡ de1c4fec-dcbd-4e93-a939-9297afaf2d97
𝒮 = sample(𝒦, 1_000_000, replace=false)

# ╔═╡ 77965f06-ed4c-4b42-90cf-1c34122f2037
md"""
Due to compaction and other physical processes, grains are ellipsoids with major and minor axes. We would like to characterize these axes for different sections (or planes) of the Ketton rock. For example, we can consider planes that are normal to the vertical direction:
"""

# ╔═╡ d6c60d3e-51db-4d36-8b20-73e78e8a892b
begin
	𝒫 = partition(𝒮, PlanePartitioner((0.,0.,1.)))
	
	plot(𝒫[1])
	plot!(𝒫[2])
	plot!(𝒫[3])
end

# ╔═╡ 2ab80306-b3e1-431a-bdb1-1fe8e379b411
md"""
and efficiently estimate the variogram on these planes:
"""

# ╔═╡ bab9a493-2338-444b-9d45-6de34dcf4114
@time EmpiricalVariogram(𝒫, :grain)

# ╔═╡ 91753cb1-8529-4dff-ac1f-22b491fae8be
md"""
These two steps can be performed together using a `PlanarVariogram`, where the first argument is the normal direction to the planes of interest:
"""

# ╔═╡ 6eb63d85-322a-460c-8905-f1e82109558f
γₑ = PlanarVariogram((0.,0.,1.), 𝒮, :grain)

# ╔═╡ f6c8e110-1a83-11eb-3422-734f673f07b0
plot(γₑ)

# ╔═╡ 76f11c9d-c19a-4416-8f76-eb022711bbbc
md"""
In order to evaluate the variogram at any distance (or lag), we fit a theoretical variogram model:
"""

# ╔═╡ 02074c10-1a84-11eb-07f4-3b896aacd1ab
γₜ = fit(ExponentialVariogram, γₑ)

# ╔═╡ 5bb91508-ce5e-4e62-865f-b77a75b38b04
begin
	plot(γₑ, label="empirical")
	plot!(γₜ, 0., 30., c=:green, label="theoretical")
end

# ╔═╡ f3bfaadd-edda-4402-8e11-9476aa302179
md"""
Finally, we estimate the marginal probability $p$ using the proportion of *grain* in the rock, and scale our variogram accordingly:
"""

# ╔═╡ 1561362f-f577-4a12-99ce-281c067356f9
p = mean(𝒮[:grain])

# ╔═╡ 0d6218b0-1a84-11eb-21ac-29e7bcea5fcb
γₚ = (1/p) * γₜ

# ╔═╡ 01026525-2eee-4217-bb9d-6e07ee57b5c1
md"""
### Probability versus distance

The scaled variogram $\gamma_p(h)$ is the probability $Pr(Z_v=0|Z_u=1)$ as a function of the distance $h = ||u-v||$:
"""

# ╔═╡ 592dc540-9905-4372-be16-70fac9d28df7
begin
	plot(γₚ, label="Pr(Zᵥ = 0 | Zᵤ = 1)")
	hline!([sill(γₚ)], ls=:dash, label="1 - p (proportion of non-grain)")
	vline!([range(γₚ)], c=:purple, ls=:dash, label="average grain radius")
end

# ╔═╡ ae66baba-2a86-41f5-9f11-0bdf4548d5db
md"""
The probability levels out at some positive distance (known as the variogram range), after which the chance of encountering a *non-grain* at location $v$ given a *grain* at location $u$ is simply the marginal probability of *non-grain* $1-p$ (i.e. the variogram sill):
"""

# ╔═╡ a74e3730-4e10-490e-8b1a-f580d31a6ff7
1 - p

# ╔═╡ 0a7a54bd-bb1d-4510-9e4d-cb3dfdb5ded9
md"""
Additionally, we define the "average grain radius" on the planes of interest as the variogram range.
"""

# ╔═╡ 2fe75541-7860-4738-acdb-abb2ee01f88a
md"""
### Average grain shape

We can create a function to compute the average grain radius for any set of planes given a normal direction:
"""

# ╔═╡ a7de2c67-e03b-4dcf-8b63-4d2846dace89
function radius(normal)
    γₑ = PlanarVariogram(normal, 𝒮, :grain)
    γₜ = fit(ExponentialVariogram, γₑ)
    p = mean(𝒮[:grain])
    range((1/p) * γₜ)
end

# ╔═╡ 11b1039e-1a88-11eb-387c-f395bd21662d
begin
	normal = (1.,0.,0.)
	md"$normal → $(radius(normal))"
end

# ╔═╡ 64a5e112-1a84-11eb-1c53-d128549d1c62
begin
	normal1 = (0.,1.,0.)
	md"$normal1 → $(radius(normal1))"
end

# ╔═╡ 30fdc030-1a89-11eb-3f70-21d77c5188ea
begin
	normal2 = (0.,0.,1.)
	md"$normal2 → $(radius(normal2))"
end

# ╔═╡ 2e04841e-de35-4d9b-870f-636b1492150c
md"""
We notice that the average grain radius is larger on horizontal sections of the rock compared with vertical sections. To visualize this anisotropy, we choose a set of vertical planes and estimate the variogram for all directions on these planes:
"""

# ╔═╡ a8568e21-dcfe-412d-8b5f-e6d75d0a60eb
γθ = @time EmpiricalVarioplane(𝒮, :grain, maxlag=30., normal=(1.,0.,0.))

# ╔═╡ 1f24d134-6b35-417a-b21b-9ceabc40f9bb
md"""
The following plot depicts the "average grain shape" on vertical sections of the rock:
"""

# ╔═╡ d1342474-b475-4c82-9004-4cfa2c43432d
plot(γθ, model=ExponentialVariogram)

# ╔═╡ 084ffdad-0fb2-4df5-ab78-c0e1c8bb850c
md"""
## Remarks

- Variograms can be useful to efficiently estimate two-point spatial probabilities.
- We've shown that the variogram range represents the "average grain radius".
- GeoStats.jl provides variogram estimators that can handle more than 1M points.
- For more details, please check:
  - *Hoffimann & Zadrozny. 2019.* [Efficient variography with partition variograms](https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""

# ╔═╡ Cell order:
# ╠═07589519-3c07-4a44-a3a2-f312213687c1
# ╠═4a66263a-7fba-4294-89d3-55a5cc82bd27
# ╠═15ed2ae0-1a81-11eb-13c3-01781cc979f6
# ╟─10ac0071-273b-4c9e-9b97-c94071dc826f
# ╟─644bcfa0-1a83-11eb-2a92-03074995e820
# ╟─54a3f000-1a83-11eb-2d1c-55fbe1506269
# ╠═49634762-156e-4f01-a4b9-a4b2e74597fa
# ╟─d5870dd6-a294-4cda-af18-88f51abbeaa5
# ╟─cf926398-2cf7-409f-b5c1-c0a25955098a
# ╠═de1c4fec-dcbd-4e93-a939-9297afaf2d97
# ╟─77965f06-ed4c-4b42-90cf-1c34122f2037
# ╠═d6c60d3e-51db-4d36-8b20-73e78e8a892b
# ╟─2ab80306-b3e1-431a-bdb1-1fe8e379b411
# ╠═bab9a493-2338-444b-9d45-6de34dcf4114
# ╟─91753cb1-8529-4dff-ac1f-22b491fae8be
# ╠═6eb63d85-322a-460c-8905-f1e82109558f
# ╠═f6c8e110-1a83-11eb-3422-734f673f07b0
# ╟─76f11c9d-c19a-4416-8f76-eb022711bbbc
# ╠═02074c10-1a84-11eb-07f4-3b896aacd1ab
# ╠═5bb91508-ce5e-4e62-865f-b77a75b38b04
# ╟─f3bfaadd-edda-4402-8e11-9476aa302179
# ╠═1561362f-f577-4a12-99ce-281c067356f9
# ╠═0d6218b0-1a84-11eb-21ac-29e7bcea5fcb
# ╟─01026525-2eee-4217-bb9d-6e07ee57b5c1
# ╠═592dc540-9905-4372-be16-70fac9d28df7
# ╟─ae66baba-2a86-41f5-9f11-0bdf4548d5db
# ╠═a74e3730-4e10-490e-8b1a-f580d31a6ff7
# ╟─0a7a54bd-bb1d-4510-9e4d-cb3dfdb5ded9
# ╟─2fe75541-7860-4738-acdb-abb2ee01f88a
# ╠═a7de2c67-e03b-4dcf-8b63-4d2846dace89
# ╟─11b1039e-1a88-11eb-387c-f395bd21662d
# ╟─64a5e112-1a84-11eb-1c53-d128549d1c62
# ╟─30fdc030-1a89-11eb-3f70-21d77c5188ea
# ╟─2e04841e-de35-4d9b-870f-636b1492150c
# ╠═a8568e21-dcfe-412d-8b5f-e6d75d0a60eb
# ╟─1f24d134-6b35-417a-b21b-9ceabc40f9bb
# ╠═d1342474-b475-4c82-9004-4cfa2c43432d
# ╟─084ffdad-0fb2-4df5-ab78-c0e1c8bb850c
