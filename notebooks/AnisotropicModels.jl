### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# ╔═╡ e08cb0c6-a704-4ab8-8c51-484e5cb770e4
using Pkg; Pkg.instantiate(); Pkg.precompile()

# ╔═╡ 20409989-41a8-4366-88ad-c562fff4dfa1
using Random; Random.seed!(2000);

# ╔═╡ 7ff34075-f185-4b75-809b-54710f3c9722
begin
	using GeoStats
	using Distances
	using Plots
	gr(size=(700,400))
end;

# ╔═╡ 3f7c2cb0-4847-4581-82fa-c1e577142c99
md"""
# Anisotropic models

In this tutorial, we demonstrate how to perform estimation with anisotropic variograms.
"""

# ╔═╡ 9abef139-b525-4d4a-97f4-d4ea9185af0c
md"""
## Ellipsoid distance

Anisotropy can be thought of as a deformation of space with an ellipsoid distance. The semiaxes of the ellipsoid determine the preferential directions of the field and their lengths characterize the anisotropy ratio. In GeoStats.jl, all variogram models (empirical and theoretical) support a custom distance function that
can be used to model anisotropy.

$
\newcommand{x}{\mathbf{x}}
\newcommand{y}{\mathbf{y}}
$

A variogram object $\gamma$ can be evaluated as an isotropic model $\gamma(h)$ or as a (possibly) anisotropic model $\gamma(\x,\y)$. For the Euclidean distance (the default), these two operations match $\gamma(\x,\y) = \gamma(h)$ in all directions:
"""

# ╔═╡ 740bd290-19bc-11eb-0b38-a19251f5c4a1
begin
	γₑ = GaussianVariogram()
	γₑ([1.,0.], [0.,0.]) ≈ γₑ(1.)
end

# ╔═╡ 1b4bc522-c693-4a94-898c-b444e087f001
md"""
If instead of an Euclidean ball, we use an ellipsoid with different semiaxes, the operation $\gamma(x,y)$ becomes a function of the direction $x-y$. For example, we can create an ellipsoid distance aligned with the coordinate system where the major semiaxis has twice the size of the minor semiaxis:
"""

# ╔═╡ ac70ed00-19bc-11eb-08d5-5d35918b469c
begin
	γₐ = GaussianVariogram(distance=Ellipsoidal([2.,1.],[0.]))
	γₐ([1.,0.],[0.,0.]) ≠ γₐ([0.,1.],[0.,0.])
end

# ╔═╡ 55647561-3f1b-475a-b457-bcc4dabc223f
md"""
## Effects on estimation

Now that we know how to construct anisotropic variograms, we can investigate the effect of varying the anisotropy ratio and alignement angle on estimation results.

We start by generating some random data:
"""

# ╔═╡ 322e0710-19be-11eb-32c9-9920da6d416a
begin
	dim, nobs = 2, 50
	X = 100*rand(dim, nobs)
	z = rand(nobs)
	𝒮 = georef((z=z,), X)
	plot(𝒮)
end

# ╔═╡ 81077b35-85dd-4f76-8986-1223c1965d08
md"""
and by defining an estimation problem:
"""

# ╔═╡ 08efc0b2-242a-4eab-accc-17f90949469b
𝒟 = RegularGrid(100, 100)

# ╔═╡ 4ca2a420-19be-11eb-2761-37be6ab6bdd2
𝒫 = EstimationProblem(𝒮, 𝒟, :z)

# ╔═╡ e6c5a087-373d-4d67-9b5e-24b8ce0e70c7
md"""
First, we vary the anisotropy ratio with an ellipsoid that is aligned with the coordinate system:
"""

# ╔═╡ 374d8b78-544b-4a59-99e3-f39996297edc
anim = @animate for r in range(1, stop=10., length=10)
    d = Ellipsoidal([r,1.], [0.])
    
    γ = GaussianVariogram(range=5., distance=d)
    
    s = solve(𝒫, Kriging(:z => (variogram=γ,)))
    
    plot(s, size=(800,400))
end;

# ╔═╡ a597d6f0-19c7-11eb-1d7b-8fc5d3bed81a
gif(anim, "figs/anisotropy_ratio.gif", fps=1)

# ╔═╡ b0ab7f62-914c-4cb9-a56e-1ee1ffb3295c
md"""
Second, we fix the anisotropy ratio and vary the alignment angle:
"""

# ╔═╡ 9cd702b9-ed26-445e-ac01-8eba860753da
anim1 = @animate for θ in range(0, stop=2π, length=10)
	d = Ellipsoidal([10.,1.], [θ])

	γ = GaussianVariogram(range=5., distance=d)

	s = solve(𝒫, Kriging(:z => (variogram=γ,)))

	plot(s, size=(800,400))
end;

# ╔═╡ f89daa10-1a20-11eb-339a-11253411a516
gif(anim1, "figs/anisotropy_angle.gif", fps=1)

# ╔═╡ 81a64650-f46f-41aa-b0aa-5b4e20765317
md"""
This experiment can be extended to 3D with the only difference being that ellipsoids therein are defined by 3 semiaxes and 3 angles. For example, the Euclidean distance in 3D can be recovered with a degenerated ellipsoid with equal semiaxes (i.e. sphere):
"""

# ╔═╡ 6d347750-1a21-11eb-2d8f-81cdd3950b1a
begin
	d₁ = Ellipsoidal([1.,1.,1.],[0.,0.,0.])
	d₂ = Euclidean()
	
	a, b = rand(3), rand(3)
	
	evaluate(d₁, a, b) ≈ evaluate(d₂, a, b)
end

# ╔═╡ a13ad039-5389-4723-8848-0a6fecada557
md"""
## Conclusions

- Geometric anisotropy can be easily modeled with an `Ellipsoidal` distance

- GeoStats.jl recognizes any distance following the [Distances.jl](https://github.com/JuliaStats/Distances.jl) API
"""

# ╔═╡ Cell order:
# ╟─e08cb0c6-a704-4ab8-8c51-484e5cb770e4
# ╟─20409989-41a8-4366-88ad-c562fff4dfa1
# ╠═7ff34075-f185-4b75-809b-54710f3c9722
# ╟─3f7c2cb0-4847-4581-82fa-c1e577142c99
# ╟─9abef139-b525-4d4a-97f4-d4ea9185af0c
# ╠═740bd290-19bc-11eb-0b38-a19251f5c4a1
# ╟─1b4bc522-c693-4a94-898c-b444e087f001
# ╠═ac70ed00-19bc-11eb-08d5-5d35918b469c
# ╟─55647561-3f1b-475a-b457-bcc4dabc223f
# ╠═322e0710-19be-11eb-32c9-9920da6d416a
# ╟─81077b35-85dd-4f76-8986-1223c1965d08
# ╠═08efc0b2-242a-4eab-accc-17f90949469b
# ╠═4ca2a420-19be-11eb-2761-37be6ab6bdd2
# ╟─e6c5a087-373d-4d67-9b5e-24b8ce0e70c7
# ╠═374d8b78-544b-4a59-99e3-f39996297edc
# ╠═a597d6f0-19c7-11eb-1d7b-8fc5d3bed81a
# ╟─b0ab7f62-914c-4cb9-a56e-1ee1ffb3295c
# ╠═9cd702b9-ed26-445e-ac01-8eba860753da
# ╠═f89daa10-1a20-11eb-339a-11253411a516
# ╟─81a64650-f46f-41aa-b0aa-5b4e20765317
# ╠═6d347750-1a21-11eb-2d8f-81cdd3950b1a
# ╟─a13ad039-5389-4723-8848-0a6fecada557
