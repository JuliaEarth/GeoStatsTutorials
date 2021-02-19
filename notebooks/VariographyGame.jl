### A Pluto.jl notebook ###
# v0.12.20

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

# ╔═╡ d9955fc8-20e2-11eb-22cf-4133b4a9ae4b
begin
	using Distributed
	pids = [myid()]
	
	md"""
	Running on processes: $pids
	
	Use `pids = addprocs(n)` to run the notebook with `n` parallel processes.
	"""
end

# ╔═╡ 07255bc8-20e3-11eb-1600-07e828113814
@everywhere pids begin
	using Pkg; Pkg.activate(@__DIR__)
	Pkg.instantiate(); Pkg.precompile()
end

# ╔═╡ 0be25d28-20e3-11eb-0c18-99a30c8b9291
@everywhere pids begin
	# packages used in this notebook
	using Statistics
	using GeoStats
	using Images
	using PlutoUI
	using JSON
	using Printf
	
	# default plot settings
	using Plots; gr(size=(650,400), ms=1, c=:cividis)
	
	# make sure that results are reproducible
	using Random; Random.seed!(2020)
end

# ╔═╡ 2a67644a-20e9-11eb-02f2-71def5bd3fe5
md"""
# Variography

Suppose we are given a real-world satellite image, and are asked to describe its spatial correlation. Although this task seems simple, many people get confused by the concept.

In this tutorial we will address two issues:

1. *What is spatial correlation?*
2. *How to describe it properly?*
"""

# ╔═╡ 56478930-20e3-11eb-1d50-97c7bf7499b6
begin
	json = download(
		"https://raw.githubusercontent.com/limhenry/earthview/master/earthview.json"
	)
	
	data = JSON.parse(read(json, String))
	
	nimgs = length(data)
	
	md"""
	Image ID: $(@bind id Slider(1:nimgs, default=68, show_value=true))
	
	Suggested: 68, 810, 487, 214
	"""
end

# ╔═╡ 32899148-20e6-11eb-06e4-25de9b79df3a
begin
	item = data[id]
	country = item["country"]
	region = item["region"]
	imgurl = item["image"]
	gmaps = item["map"]
	
	location = isempty(region) ? country : country*", "*region
	
	file = download(imgurl)
	img  = load(file)
	pix  = join(size(img), "×")
	
	HTML("""
	<a href=$gmaps>$location</a> ($pix pixels)
	<a href=$gmaps><img src=$imgurl></a>
	<p align="right">Source: <a href=https://github.com/limhenry/earthview>https://github.com/limhenry/earthview</a></p>
	""")
end

# ╔═╡ bcc7d906-20e8-11eb-3d30-0920f09d6d9a
md"""
## What is spatial correlation?

What does it mean to say that an image looks *"noisy"* or *"smooth"*? How are these concepts connected to the concept of *"spatial correlation"* that we hear all the time?

First, let's recall the definition of the Pearson (a.k.a. linear) correlation coefficient between two random variables ``X`` and ``Y`` that have been sampled jointly ``\{(x_1,y_1),(x_2,y_2),\ldots,(x_n,y_n)\}``:

$$\newcommand{\cov}{\text{cov}}
\newcommand{\xbar}{\bar{x}}
\newcommand{\ybar}{\bar{y}}
\rho_{XY} = \frac{\cov(X,Y)}{\sigma_X\sigma_Y}
            \approx \frac{\sum_1^n (x_i-\xbar)(y_i-\ybar)}
               {\sqrt{\sum_1^n (x_i-\xbar)^2}\sqrt{\sum_1^n (y_i-\ybar)^2}}$$

Imagine that ``X`` and ``Y`` are different layers of the satellite image, and that the index ``i`` refers to a pixel. For example, consider that ``X`` and ``Y`` are the red and blue channels of the image:
"""

# ╔═╡ 7d5dc4fa-20fb-11eb-0d94-c587609d715e
md"""
X = $(@bind xselection Select(["red", "green", "blue"], default="red"))
Y = $(@bind yselection Select(["red", "green", "blue"], default="blue"))
"""

# ╔═╡ c77808d2-20fd-11eb-204b-c9c637194c59
begin
	channel = Dict("red"=>red, "green"=>green, "blue"=>blue)
	
	X = channel[xselection].(img)
	Y = channel[yselection].(img)
	
	Gray.([X Y])
end

# ╔═╡ 2e6609ca-2102-11eb-0265-b75e14f95935
md"""
The corresponding samples are obtained by looping over the pixels:
"""

# ╔═╡ 39d6774e-2105-11eb-3dac-0bfd44d54bc5
begin
	I = georef((X=X, Y=Y))
	S = sample(I, 5000)
	hscatter(S, :X, :Y)
end

# ╔═╡ e3ebc034-2105-11eb-3b1c-078a372c0e22
ρ = cor(vec(X), vec(Y))

# ╔═╡ c05b6ba0-2115-11eb-2a70-65e56ea77fd5
md"""
The correlation coefficient tells us very little about the relationship between the variables, but at least it gives us an idea of *association* at a single pixel.

### h-scatter

Now, consider samples of ``X`` and ``Y`` obtained from two (possibly) different pixels ``i \ne j``:

$$\{(x_i,y_j)\}_{i,j=1,2,\ldots,n}$$

In particular, consider pixels that are aligned with a direction of interest (e.g. West-East), and that are ``h`` units apart. Slide an arrow ``i\underbrace{\longrightarrow}_{h \text{ units}} j`` of fixed length over the image to produce samples:
"""

# ╔═╡ 4af354f6-212f-11eb-39a8-f11b5d90d517
Lx, Ly = 200, 200;

# ╔═╡ 1755b762-2125-11eb-0399-1f8e5c9042ec
md"""
Arrow length: $(@bind h Slider(1:Lx-1, default=Lx÷2))
"""

# ╔═╡ c1b263f4-212a-11eb-15c2-dbf154f593cf
md"""
Slide controls:
$(@bind ix Slider(0:Lx-h-1, default=5))
$(@bind iy Slider(0:Ly-1, default=Ly-8))
"""

# ╔═╡ ded79e0c-2124-11eb-180e-7da2f0cef9fe
begin
	XX = X[1:Lx,1:Ly]
	YY = Y[1:Lx,1:Ly]
	
	xi = @sprintf "%.2f" X[ix+1,iy+1]
	xj = @sprintf "%.2f" X[ix+h+1,iy+1]
	yi = @sprintf "%.2f" Y[ix+1,iy+1]
	yj = @sprintf "%.2f" Y[ix+h+1,iy+1]
	
	p1 = heatmap(XX, yflip=true, aspect_ratio=:equal, axis=false, ticks=false,
		         colorbar=false, legend=false, c=:grays,
		         title="(xᵢ=$xi, xⱼ=$xj)")
	plot!([(ix+2,iy),(ix+h-2,iy)], yflip=false, lw=3, arrow=true, c=:cyan)
	annotate!(ix, iy, text("i", 18, :white))
	annotate!(ix+h, iy, text("j", 18, :white))
	
	p2 = heatmap(YY, yflip=true, aspect_ratio=:equal, axis=false, ticks=false,
		         colorbar=false, legend=false, c=:grays,
		         title="(yᵢ=$yi, yⱼ=$yj)")
	plot!([(ix+2,iy),(ix+h-2,iy)], yflip=false, lw=3, arrow=true, c=:cyan)
	annotate!(ix, iy, text("i", 18, :white))
	annotate!(ix+h, iy, text("j", 18, :white))
	
	plot(p1, p2)
end

# ╔═╡ 19eda198-212e-11eb-023b-59ea6cec0547
md"""
Denote the variable at the head of the arrow by ``H`` (either ``X`` or ``Y``), and the variable at the tail of the arrow by ``T`` (either ``X`` or ``Y``). Let's take a look at samples for a given arrow length ``h`` without restricting ourselves to any particular direction. The following plot is known as the *h-scatter* plot:
"""

# ╔═╡ 0274f2b0-21c8-11eb-3ffe-3b2ea8ba1002
𝒟 = georef((X=X, Y=Y))

# ╔═╡ 81317d2e-2134-11eb-3f53-fb7f0b56178d
𝒮 = sample(𝒟, 5000)

# ╔═╡ bcae9fb2-2134-11eb-1f6f-3d3d65209557
md"""
H = $(@bind hselection Select(["X","Y"], default="X"))
T = $(@bind tselection Select(["X","Y"], default="Y"))
lag = $(@bind lag Slider(0:500,show_value=true))
"""

# ╔═╡ 7ae0aa56-21c7-11eb-2908-cbdbccf37b39
begin
	H = Symbol(hselection)
	T = Symbol(tselection)
end;

# ╔═╡ 0f2a7ce4-2139-11eb-16e7-eb14a0e1688e
hscatter(𝒮, H, T, lag=lag)

# ╔═╡ c029a4b0-2202-11eb-1eef-274833b44bba
md"""
### Correlogram vs. Variogram

The *spatial correlation* is the Pearson correlation ``\rho_{HT}`` between ``H`` and ``T``. When studied as a function of the arrow length (or lag) ``h``, it is known as the *correlogram*:
"""

# ╔═╡ be9e549c-22e8-11eb-0c3b-a3037ecb2655
plot(h->exp(-h), 0, 10, xlabel="h", ylabel="ρ(h)", label="correlogram", size=(700,300))

# ╔═╡ b53b130e-22e8-11eb-046f-f550407335f7
md"""
- The correlogram ``\rho(h)`` is often a non-increasing function.
- At ``h=0`` we recover the non-spatial (often maximum) correlation.
- In practice we have ``\rho(h) \to 0`` as ``h \to \infty``.

The terms *auto-correlogram* (``H = T``) and *cross-correlogram* (``H \ne T``) are also encountered in the literature to differentiate the various spatial correlations in the multivariate case. Correspondingly, the terms *auto-covariance* and *cross-covariance* are encountered for the covariance ``\cov(h)``.

We now consider an alternative statistic of association known as the (auto-)variogram:

$$\gamma_X(h) \approx \frac{1}{2|N(h)|}\sum_{N(h)}(x_i-x_j)^2$$

where ``N(h) = \Big\{(i,j): i\underbrace{\longrightarrow}_{h \text{ units}} j\Big\}``. And similarly, in the multivariate case the (cross-)variogram:

$$\gamma_{XY}(h) \approx \frac{1}{2|N(h)|}\sum_{N(h)}(x_i-x_j)(y_i-y_j)$$
"""

# ╔═╡ 42ab6054-22e9-11eb-0b60-8f7fd5bc94cb
plot(h->1-exp(-h), 0, 10, xlabel="h", ylabel="γ(h)", label="variogram", size=(700,300))

# ╔═╡ 37684bc6-22e9-11eb-35ea-5354c205dd19
md"""
- The variogram ``\gamma_*(h)`` does not involve the means ``\xbar`` and ``\ybar``.
- It is a measure of "spread" of the h-scatter plot.
- At ``h = 0`` we have that ``\gamma_*(h) = 0``.
- In practice we have ``\gamma_*(h) \to \sigma^2`` as ``h \to \infty``.

The term *semi-variogram* is also encountered in the literature to emphasize the ``\frac{1}{2}`` term in the formula, but various authors consider it to be outdated terminology.

In many practical cases the variogram is an inverted version of the covariance:

$$\gamma(h) = \cov(0) - \cov(h)$$

but it can also be used with random processes that are only intrinsically 2nd-order stationary. For more details, please refer to [Chapter 4 - Chilès & Delfiner 2012](https://onlinelibrary.wiley.com/doi/abs/10.1002/9781118136188.ch4).
"""

# ╔═╡ 2fdc9392-2214-11eb-2bd5-974b00d0d9eb
md"""
## The four elements of the variogram

Now that we have explained what is spatial correlation and how it is connected to the variogram, we can build intuition with theoretical models. Widely used models have at least four elements controlling spatial variability:

1. The *range* controls the average radius (or correlation length) of the "blobs" in the image.
2. The *sill* controls the amplitude (or height) of the "blobs", i.e. mountains and valleys.
3. The *nugget* controls an additional variability at scales smaller than the pixel.
4. The *model* type controls the behavior near the origin, i.e. short-scale variability.
"""

# ╔═╡ 294d57f6-2290-11eb-3383-857c0013d5a9
md"""
range = $(@bind r Slider(1:25, default=10, show_value=true))
"""

# ╔═╡ c95ece92-228f-11eb-28c5-7571c47dea7a
md"""
sill = $(@bind s Slider(0.5:0.1:1, default=0.7, show_value=true))
"""

# ╔═╡ 3fd1c186-2290-11eb-2045-3b2196cb0692
md"""
nugget = $(@bind n Slider(0:0.05:0.2, default=0.1, show_value=true))
"""

# ╔═╡ 85fe6280-2290-11eb-26ef-618d06d9eb0a
md"""
model = $(@bind m Select(["Gaussian","Spherical","Exponential"]))
"""

# ╔═╡ 0ec036f4-2299-11eb-2344-75a39828fde1
begin
	xs = rand(0.0:1.0:99.0, 100)
	ys = rand(0.0:1.0:24.0, 100)
	zs = randn(100)
		
	sdata = georef((X=zs,), collect(zip(xs,ys)))
end;

# ╔═╡ 939c95dc-2291-11eb-33e0-23ddb7b27f65
begin
	model = Dict("Spherical"=>SphericalVariogram,
		         "Gaussian"=>GaussianVariogram,
		         "Exponential"=>ExponentialVariogram)
	
	g = model[m](sill=Float64(s), range=Float64(r), nugget=Float64(n))
	
	gplot = plot(g, 0, 25, c=:black, ylim=(0,1),
		         legend=:topright, size=(650,300))
	vline!([r], c=:grey, ls=:dash, primary=false)
	annotate!(r-2, 1, "range")
	hline!([s], c=:brown, ls=:dash, primary=false)
	annotate!(23, s+0.05, "sill")
	if n > 0
		hline!([n], c=:orange, ls=:dash, primary=false)
		annotate!(23, n+0.05, "nugget")
	end
	gplot
end

# ╔═╡ 2069cac4-2294-11eb-0f72-ef9c871b6b40
begin
	P   = SimulationProblem(sdata, RegularGrid(100,25), :X, 1)
	
	LU  = LUGS(:X => (variogram=g,))
	
	sol = solve(P, LU)
	
	plot(sol, c=:cividis, clim=(-3,3), size=(700,200))
	plot!(sdata, markersize=2, markershape=:square,
		  markerstrokecolor=:white, markerstrokewidth=3)
end

# ╔═╡ 50bc2970-22a4-11eb-34d1-dff6778b93c7
md"""
## The variography game

Given an image, can you describe its spatial correlation?
"""

# ╔═╡ ba0ee796-22a7-11eb-1cba-c7ad39e63390
@bind sampled Button("SAMPLE")

# ╔═╡ 78218b52-22ab-11eb-2fc9-c53dcf9ff6c7
md"Show answer: $(@bind answer CheckBox(default=false))"

# ╔═╡ f36e26de-22a7-11eb-0dd3-11b98dbd1e8c
begin
	sampled
	
	Model  = rand([GaussianVariogram, SphericalVariogram, ExponentialVariogram])
	range  = rand(0.0:1.0:100.)
	sill   = rand(0.1:0.1:1.0)
	nugget = rand(0.0:0.1:sill)
	γ = Model(range=range, sill=sill, nugget=nugget)
end;

# ╔═╡ 8342770a-22a9-11eb-2f9d-ff9b0b0cb550
begin
	problem  = SimulationProblem(RegularGrid(600,300), :X=>Float64, 1)
	
	solver   = FFTGS(:X => (variogram=γ,))
	
	solution = solve(problem, solver)
	
	plot(solution)
end

# ╔═╡ 321e88fc-22ac-11eb-1fc7-834c9282f423
answer ? γ : nothing

# ╔═╡ bc12255e-22ad-11eb-27b3-b191a5455c4c
md"""
## Remarks

- Understanding spatial correlation is key in geospatial applications.
- If you work with geospatial data, we highly recommend the variography game.
- [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) compiles a great list of variogram models to explore:
"""

# ╔═╡ ede8f78a-22c4-11eb-1ad6-4991aaf69684
begin
	vplot = plot()
	for V in subtypes(Variogram)
		if V ≠ NestedVariogram
			plot!(V(), c=:auto)
		end
	end
	vplot
end

# ╔═╡ Cell order:
# ╟─d9955fc8-20e2-11eb-22cf-4133b4a9ae4b
# ╟─07255bc8-20e3-11eb-1600-07e828113814
# ╠═0be25d28-20e3-11eb-0c18-99a30c8b9291
# ╟─2a67644a-20e9-11eb-02f2-71def5bd3fe5
# ╟─32899148-20e6-11eb-06e4-25de9b79df3a
# ╟─56478930-20e3-11eb-1d50-97c7bf7499b6
# ╟─bcc7d906-20e8-11eb-3d30-0920f09d6d9a
# ╟─7d5dc4fa-20fb-11eb-0d94-c587609d715e
# ╟─c77808d2-20fd-11eb-204b-c9c637194c59
# ╟─2e6609ca-2102-11eb-0265-b75e14f95935
# ╟─39d6774e-2105-11eb-3dac-0bfd44d54bc5
# ╠═e3ebc034-2105-11eb-3b1c-078a372c0e22
# ╟─c05b6ba0-2115-11eb-2a70-65e56ea77fd5
# ╟─4af354f6-212f-11eb-39a8-f11b5d90d517
# ╟─c1b263f4-212a-11eb-15c2-dbf154f593cf
# ╟─1755b762-2125-11eb-0399-1f8e5c9042ec
# ╟─ded79e0c-2124-11eb-180e-7da2f0cef9fe
# ╟─19eda198-212e-11eb-023b-59ea6cec0547
# ╠═0274f2b0-21c8-11eb-3ffe-3b2ea8ba1002
# ╠═81317d2e-2134-11eb-3f53-fb7f0b56178d
# ╟─bcae9fb2-2134-11eb-1f6f-3d3d65209557
# ╟─7ae0aa56-21c7-11eb-2908-cbdbccf37b39
# ╠═0f2a7ce4-2139-11eb-16e7-eb14a0e1688e
# ╟─c029a4b0-2202-11eb-1eef-274833b44bba
# ╟─be9e549c-22e8-11eb-0c3b-a3037ecb2655
# ╟─b53b130e-22e8-11eb-046f-f550407335f7
# ╟─42ab6054-22e9-11eb-0b60-8f7fd5bc94cb
# ╟─37684bc6-22e9-11eb-35ea-5354c205dd19
# ╟─2fdc9392-2214-11eb-2bd5-974b00d0d9eb
# ╟─294d57f6-2290-11eb-3383-857c0013d5a9
# ╟─c95ece92-228f-11eb-28c5-7571c47dea7a
# ╟─3fd1c186-2290-11eb-2045-3b2196cb0692
# ╟─85fe6280-2290-11eb-26ef-618d06d9eb0a
# ╟─0ec036f4-2299-11eb-2344-75a39828fde1
# ╟─939c95dc-2291-11eb-33e0-23ddb7b27f65
# ╟─2069cac4-2294-11eb-0f72-ef9c871b6b40
# ╟─50bc2970-22a4-11eb-34d1-dff6778b93c7
# ╟─ba0ee796-22a7-11eb-1cba-c7ad39e63390
# ╟─78218b52-22ab-11eb-2fc9-c53dcf9ff6c7
# ╟─f36e26de-22a7-11eb-0dd3-11b98dbd1e8c
# ╟─8342770a-22a9-11eb-2f9d-ff9b0b0cb550
# ╟─321e88fc-22ac-11eb-1fc7-834c9282f423
# ╟─bc12255e-22ad-11eb-27b3-b191a5455c4c
# ╟─ede8f78a-22c4-11eb-1ad6-4991aaf69684
