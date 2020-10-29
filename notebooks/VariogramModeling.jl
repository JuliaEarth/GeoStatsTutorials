### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ a89ee010-1840-11eb-1516-e11f57e7213a
using Pkg; Pkg.instantiate(); Pkg.precompile()

# ╔═╡ b9a194c0-1840-11eb-250d-e5c7477e7899
using Random; Random.seed!(2017);

# ╔═╡ 11cd3600-184a-11eb-31ee-0b05bb76e7da
using GeoStats

# ╔═╡ 198427f0-184a-11eb-029f-995776ad5462
using GeoStatsImages

# ╔═╡ 1afce542-184a-11eb-048e-85fea03cbfca
using Plots

# ╔═╡ 2c733260-1841-11eb-19c7-c1934da53ce3
using PlutoUI

# ╔═╡ bd2c9040-1840-11eb-3b28-7f9d8681a235
md"""# Variogram modeling

In this tutorial, we illustrate one possible workflow for fitting variograms with GeoStats.jl."""

# ╔═╡ cb75ca40-1840-11eb-07b6-65ffdf8b1480
md"""## The data

Let's consider a simple 2D problem in which properties of a field are sampled at random locations. We will be using an image of the Walker Lake in Nevada available in [GeoStatsImages.jl](https://github.com/JuliaEarth/GeoStatsImages.jl) as our field:"""

# ╔═╡ 206af8a0-184a-11eb-2223-616291fcd27e
𝒟 = geostatsimage("WalkerLake")

# ╔═╡ 241f9c30-184a-11eb-04ef-2139b4670992
plot(𝒟)

# ╔═╡ f90ceec2-1840-11eb-02e8-777f5d33918f
md"""We sample a thousand points from the image and plot them:"""

# ╔═╡ 13459440-1841-11eb-3f1b-ef74d30c10b2
begin
	𝒮 = sample(𝒟, 1000)
	plot(𝒮)
end

# ╔═╡ 1d0fd5d0-1841-11eb-271f-ab7825a69fec
md"""## Empirical variogram

We consider a simple omnidirectional variogram. Variogram objects in GeoStats.jl are integrated with the Julia Plots.jl framework. We can plot them directly and interact with the result as we vary the number of lags used for estimation:"""

# ╔═╡ 2480c170-186a-11eb-021f-b7a99755a0d4
md"nlags = $(@bind nlags Slider(1:30; default=15, show_value=true))"

# ╔═╡ c0c8ad50-1841-11eb-18f7-37e327fd0ae5
plot(EmpiricalVariogram(𝒮, :Z, nlags=nlags, maxlag=200.))

# ╔═╡ c74aefd0-1841-11eb-1865-dfaae510218d
md"""Besides the variogram itself, GeoStats.jl presents the bin counts (scaled) as a measure of confidence about the estimated points. This frequency plot can be deactived by passing the option `showbins=false` to the plot command. We encourage users to keep the bin counts option activated as it has zero cost.

This empirical variogram was constructed using the Euclidean distance between data locations. We can also specify a custom distance to estimate the variogram when points are embedded on different coordinate systems. Please consult the documentation for more distance functions.

After interacting with the plot, we select a number of lags and proceed to fitting a theoretical model:"""

# ╔═╡ d35b50d0-1841-11eb-0814-676fc99dc20d
γₑ = EmpiricalVariogram(𝒮, :Z, nlags=17, maxlag=200.)

# ╔═╡ d8ab2dd0-1841-11eb-3116-0fd157960d7e
plot(γₑ)

# ╔═╡ dd031640-1841-11eb-2fb7-f750c7369758
md"""## Theoretical variogram

We will exploit the Interact.jl package again to fit a theoretical variogram to the data. Various theoretical variogram models are available in GeoStats.jl, including a composite additive model that can be used to combine different variogram types. Please consult the documentation for more details. Here we will use a simple spherical variogram:"""

# ╔═╡ 8be1d270-1868-11eb-09ec-697b71485822
md"s = $(@bind s Slider(range(0, stop=.1, length=50); default=0.05, show_value=true))"

# ╔═╡ 50706cf0-1869-11eb-375f-9112e5f435e9
md"r = $(@bind r Slider(range(0, stop=100., length=50); default=50.0, show_value=true))"

# ╔═╡ 608a9660-1869-11eb-170f-699e74bd5a47
md"n = $(@bind n Slider(range(0, stop=.1, length=50); default=0.05, show_value=true))"

# ╔═╡ 824e62e0-1869-11eb-1d9b-4b650c21564f
begin
	γₜ1 = SphericalVariogram(sill=Float64(s), range=Float64(r), nugget=Float64(n))
	
    plot(γₑ, label="empirical")
    plot!(γₜ1, label="theoretical")
end

# ╔═╡ f067afc0-1841-11eb-3d48-e9e39922429c
md"""After tuning the parameters of the theoretical variogram interactively, we plot the final result:"""

# ╔═╡ f7a13450-1841-11eb-0829-0d76cd22c6e4
begin
	γₜ2 = SphericalVariogram(sill=0.083, range=55., nugget=0.02)
	plot(γₑ, label="empirical")
	plot!(γₜ2, label="theoretical")
end

# ╔═╡ 124f09d0-1842-11eb-308f-d36056ffdc9a
md"""## Automatic fitting

An alternative option to interactive variogram modeling is automatic fitting. We can fit a specific theoretical variogram model with:"""

# ╔═╡ 1a702f40-1842-11eb-1ec9-4faf98e29997
begin
	γₜ3 = fit(SphericalVariogram, γₑ)
	
	plot(γₑ, label="empirical")
	plot!(γₜ3, label="theoretical")
end

# ╔═╡ 2f3988e0-1842-11eb-0663-2bfccebc6320
md"""or let GeoStats.jl pick the model with minimum weighted least squares error by passing the super type `Variogram`:"""

# ╔═╡ 3e5d445e-1842-11eb-05bc-d98aa1f73728
begin
	γₜ4 = fit(Variogram, γₑ)
	
	plot(γₑ, label="empirical")
	plot!(γₜ4, label="theoretical")
end

# ╔═╡ 49fabf52-1842-11eb-2266-7d8c574bb98a
md"""which in this example turns out to be:"""

# ╔═╡ 519109e0-1842-11eb-38d5-f150bb9ef89b
γₜ4

# ╔═╡ 583e05e0-1842-11eb-23fb-c3630190751c
md"""## Conclusions

- Variogram modeling is an important step in classical geostatistics. Some practioners prefer to model variograms interactively to enforce a specific type of spatial continuity, whereas others prefer to use automatic fitting procedures, which are guaranteed to minimize a given loss function.

- Both approaches are available in GeoStats.jl, and are useful under different circumstances."""

# ╔═╡ Cell order:
# ╠═a89ee010-1840-11eb-1516-e11f57e7213a
# ╠═b9a194c0-1840-11eb-250d-e5c7477e7899
# ╟─bd2c9040-1840-11eb-3b28-7f9d8681a235
# ╟─cb75ca40-1840-11eb-07b6-65ffdf8b1480
# ╠═11cd3600-184a-11eb-31ee-0b05bb76e7da
# ╠═198427f0-184a-11eb-029f-995776ad5462
# ╠═1afce542-184a-11eb-048e-85fea03cbfca
# ╠═206af8a0-184a-11eb-2223-616291fcd27e
# ╠═241f9c30-184a-11eb-04ef-2139b4670992
# ╟─f90ceec2-1840-11eb-02e8-777f5d33918f
# ╠═13459440-1841-11eb-3f1b-ef74d30c10b2
# ╟─1d0fd5d0-1841-11eb-271f-ab7825a69fec
# ╠═2c733260-1841-11eb-19c7-c1934da53ce3
# ╟─2480c170-186a-11eb-021f-b7a99755a0d4
# ╠═c0c8ad50-1841-11eb-18f7-37e327fd0ae5
# ╟─c74aefd0-1841-11eb-1865-dfaae510218d
# ╠═d35b50d0-1841-11eb-0814-676fc99dc20d
# ╠═d8ab2dd0-1841-11eb-3116-0fd157960d7e
# ╟─dd031640-1841-11eb-2fb7-f750c7369758
# ╟─8be1d270-1868-11eb-09ec-697b71485822
# ╟─50706cf0-1869-11eb-375f-9112e5f435e9
# ╟─608a9660-1869-11eb-170f-699e74bd5a47
# ╠═824e62e0-1869-11eb-1d9b-4b650c21564f
# ╟─f067afc0-1841-11eb-3d48-e9e39922429c
# ╠═f7a13450-1841-11eb-0829-0d76cd22c6e4
# ╟─124f09d0-1842-11eb-308f-d36056ffdc9a
# ╠═1a702f40-1842-11eb-1ec9-4faf98e29997
# ╟─2f3988e0-1842-11eb-0663-2bfccebc6320
# ╠═3e5d445e-1842-11eb-05bc-d98aa1f73728
# ╟─49fabf52-1842-11eb-2266-7d8c574bb98a
# ╠═519109e0-1842-11eb-38d5-f150bb9ef89b
# ╟─583e05e0-1842-11eb-23fb-c3630190751c
