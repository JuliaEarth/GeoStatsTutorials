### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 20c69a8e-1fa2-11eb-3f1e-ef154de99450
begin
  using Distributed
  pids = [myid()]

  md"""
  Running on processes: $pids

  Use `pids = addprocs(n)` to run the notebook with `n` parallel processes.
  """
end

# ╔═╡ 2399f800-1fa2-11eb-185d-53d05516dacf
@everywhere pids begin
  using Pkg; Pkg.activate(@__DIR__)
  Pkg.instantiate(); Pkg.precompile()
end

# ╔═╡ 2c2e7df6-1fa2-11eb-29a8-3b59d382267e
@everywhere pids begin
  # packages used in this notebook
  using GeoStats
  using FileIO
  using GslibIO

  # default plot settings
  using Plots; gr(size=(700,400))
end


# ╔═╡ 524c65d0-2ea9-11eb-0450-5ff97535406a
md"""
# GSLIB comparison

This tutorial shows how to run a Kriging estimation problem that homologates GSLIB kriging program *kt3d*.

## Case study

The case study used in this tutorial uses synthetic drillcore samples adapted from (https://github.com/exepulveda/geomet\_datasets/tree/master/datasets/porphyry_01).
There are 6,817 samples with clay, chalcosite, bornite and chalcopyrite content.
For simplicity, we will estimate clay content at a specific elevation (2200.5 m) to make the results easy to visualize in a two-dimensional view.

The variogram has been normalized to have a total sill of 1.0.
The nugget effect accounts for 20% of the total sill and a spherical structure does for the 80% with 
rotation angles (45.0, 10.0, -5.0) and ranges (50.0, 40.0, 10.0).

We use a homotopic search radius of 500.0 m and a maximum of 50 samples.

At this moment, GeoStats does not support block kriging, therefore, ordinary kriging at point support
is compared and a block discretization of 1,1,1 is used.
"""

# ╔═╡ 524c65d0-2ea9-11eb-0450-5ff97535406b
md"""
## GSLIB setup

GSLIB *kt3d* uses a parameter file that contains all parameter values for performing kriging. 
We have stored the samples and all GSLIB required files: parameter, inputs and results;
you can reproduce the results by downloading the GSLIB executables (source code also available)
from (http://www.gslib.com/main_gd.html) and executing the kt3d program with the provided parameter file.

For this case the GSLIB parameter file is as follows:

```
                  Parameters for KT3D
                  *******************

START OF PARAMETERS:
data3D_samples.gslib              -file with data
0  1  2  3  4  0                 -   columns for DH,X,Y,Z,var,sec var
-1.0e21   1.0e21                 -   trimming limits
0                                -option: 0=grid, 1=cross, 2=jackknife
data3D_grid.gslib                -file with jackknife data
1   2   3    0    0              -   columns for X,Y,Z,vr and sec var
3                                -debugging level: 0,1,2,3
kt3d_case.dbg                         -file for debugging output
kt3d_output.gslib                         -file for kriged output
100    0.5    1.0                  -nx,xmn,xsiz
100    0.5    1.0                  -ny,ymn,ysiz
1   2200.5    1.0                  -nz,zmn,zsiz
1    1      1                    -x,y and z block discretization
4    16                          -min, max data for kriging
0                                -max per octant (0-> not used)
250.0  250.0  250.0              -maximum search radii
0.0   0.0   0.0                 -angles for search ellipsoid
1     2.302                      -0=SK,1=OK,2=non-st SK,3=exdrift
0 0 0 0 0 0 0 0 0                -drift: x,y,z,xx,yy,zz,xy,xz,zy
0                                -0, variable; 1, estimate trend
extdrift.dat                     -gridded file with drift/mean
4                                -  column number in gridded file
1    0.2                         -nst, nugget effect
1    0.8  45.0   10.0   -5.0        -it,cc,ang1,ang2,ang3
         50.0  40.0  10.0     -a_hmax, a_hmin, a_vert
```

For running kt3d, just save the above content into a text file, for example 'kriging.par' and run it as:

```
$ kt3d kriging.par

 KT3D Version: 3.000

  data file = data3D_samples.gslib                    
  columns =            0           1           2           3           4
           0
  trimming limits =  -1.0000000E+21  1.0000000E+21
  kriging option =            0
  jackknife data file = data3D_grid.gslib                       
  columns =            1           2           3           0           0
  debugging level =            3
  debugging file = kt3d_xy_ani.dbg                         
  output file = kt3d_xy_ani.gslib                       
  nx, xmn, xsiz =          100  0.5000000       1.000000    
  ny, ymn, ysiz =          100  0.5000000       1.000000    
  nz, zmn, zsiz =            1   2200.500       1.000000    
  block discretization:           1           1           1
  ndmin,ndmax =            1          50
  max per octant =            0
  search radii =    500.0000       500.0000       500.0000    
  search anisotropy angles =   0.0000000E+00  0.0000000E+00  0.0000000E+00
  ktype, skmean =           1   2.302000    
  drift terms =            0           0           0           0           0
           0           0           0           0
  itrend =            0
  external drift file = extdrift.dat                            
  variable in external drift file =            4
  nst, c0 =            1  0.2000000    
  it,cc,ang[1,2,3];            1  0.8000000       45.00000       10.00000    
  -5.000000    
  a1 a2 a3:    50.00000       40.00000       10.00000    
 Data for KT3D: Variable number            4
   Number   =         6817
   Average  =    3.631277    
   Variance =    21.16069    
 Setting up rotation matrices for variogram and search
 Setting up super block search strategy
 
 Working on the kriging 
   currently on estimate      1000
   currently on estimate      2000
   currently on estimate      3000
   currently on estimate      4000
   currently on estimate      5000
   currently on estimate      6000
   currently on estimate      7000
   currently on estimate      8000
   currently on estimate      9000
   currently on estimate     10000

Estimated      10000 blocks 
  average    3.5655079
  variance  0.39872074


 KT3D Version: 3.000 Finished

```

This execution should take about 10 seconds for 10,000 blocks.
"""

# ╔═╡ 524c65d0-2ea9-11eb-0450-5ff97535406c
md"""
## Homologation

In order to homologate GSLIB to GeoStats, we need to work with an anisotropic distance metric for the variogram model, due to GeoStats does not implement anisotropic searches yet.

Let's write the GeoStats code.

First, we load and plot the dataset.

"""

# ╔═╡ 32e5c724-2ea9-11eb-087c-4343d7cd04f1
begin
  # load samples
  𝒮 = readgeotable("data/data3D_samples.csv", coordnames=(:x,:y,:z));
  plot(𝒮)
end

# ╔═╡ 8da256e0-2f95-11eb-2e49-77bc738016dc
md"""Second, we display basic statistics for clay."""


# ╔═╡ 9ddcfb6e-2f95-11eb-19dd-23fddd5ad1a5
begin
  clay_mean = mean(𝒮[:clay])
  clay_var = var(𝒮[:clay])
  "Statistics for clay: mean=$clay_mean and var=$clay_var"
end

# ╔═╡ 6f03023c-2eac-11eb-3500-8d4560423cdc
md"""Let's define the regular grid (100 x 100 x 1) for the specific elevation at 2,200.5 m
"""

# ╔═╡ 6e03f314-2eac-11eb-38c1-4ff4d9bc802d
begin
  𝒟 = RegularGrid((100, 100, 1),(0.5, 0.5, 2200.5),(1.0, 1.0, 1.0))
end

# ╔═╡ b50ca8e0-2eb2-11eb-1ad1-ed726e10868c
md"""and the variogram model with anisotropy:"""

# ╔═╡ e3733c0e-2eac-11eb-19cc-41c2d5ec46ee
begin
  vmodel_range = 50.0
  anisotropy = [1.0, 40.0/vmodel_range, 10.0/vmodel_range]
  angles = [45.0, 10.0, -5.0]

  gslib_distance_vmodel = aniso2distance(anisotropy, angles, convention=:GSLIB)

  nugget = 0.2
  cc = 0.8

  γ = SphericalVariogram(range=vmodel_range,nugget=nugget, sill=nugget+cc, distance=gslib_distance_vmodel)
end

# ╔═╡ 0a18aa14-2eb3-11eb-3975-0b419b7524a2
md"""The estimation problem and the solution:"""

# ╔═╡ 27c72838-2eb3-11eb-29f4-ad9b1d5ce93e
begin
  problem = EstimationProblem(𝒮, 𝒟, :clay)

  solver = Kriging(
    :clay => (variogram=γ, minneighbors=1, maxneighbors=50, neighborhood=BallNeighborhood(500.0))
  )

  solution = solve(problem, solver)

  μ, σ² = solution[:clay]

  # get the right slice
  S_μ = georef((clay=reshape(μ, 100, 100, 1)[:,:,1],), RegularGrid(100,100))
  plot(S_μ)
end

# ╔═╡ 52196222-2eb3-11eb-3141-41c53f8cae56
md"""Let's now compare the stored results by GSLIB with the results by GeoStats:"""

# ╔═╡ 71b88180-2eb3-11eb-3937-8d8f6447d7b4
begin
  #load computed results from GSLIB
  grid_gslib = GslibIO.load_legacy("data/kt3d_output.gslib", (100, 100, 1))

  μ_gslib = grid_gslib[:Estimate]
  σ²_gslib = grid_gslib[:EstimationVariance]

  # get the right slice
  S_μ_gslib = georef((clay=reshape(μ_gslib, 100, 100, 1)[:,:,1],), RegularGrid(100,100))
  plot(S_μ_gslib)
end

# ╔═╡ 0bfe2786-2eb9-11eb-3248-335c1647c09f
begin
  # compute MSE 
  μ_mse = mean((μ-μ_gslib).^2)
  σ²_mse = mean((σ²-σ²_gslib).^2)

  md"""Plots look the same and mean error for estimation is $μ_mse and for variance is $σ²_mse.
  There are almost zero due to small differences in precision used in the GSLIB output file,
  therefore we can say that we have obtained the same results.
  """
end

# ╔═╡ Cell order:
# ╟─20c69a8e-1fa2-11eb-3f1e-ef154de99450
# ╟─2399f800-1fa2-11eb-185d-53d05516dacf
# ╠═2c2e7df6-1fa2-11eb-29a8-3b59d382267e
# ╠═524c65d0-2ea9-11eb-0450-5ff97535406a
# ╠═524c65d0-2ea9-11eb-0450-5ff97535406b
# ╠═524c65d0-2ea9-11eb-0450-5ff97535406c
# ╠═32e5c724-2ea9-11eb-087c-4343d7cd04f1
# ╠═8da256e0-2f95-11eb-2e49-77bc738016dc
# ╠═9ddcfb6e-2f95-11eb-19dd-23fddd5ad1a5
# ╠═6f03023c-2eac-11eb-3500-8d4560423cdc
# ╠═6e03f314-2eac-11eb-38c1-4ff4d9bc802d
# ╠═b50ca8e0-2eb2-11eb-1ad1-ed726e10868c
# ╠═e3733c0e-2eac-11eb-19cc-41c2d5ec46ee
# ╠═0a18aa14-2eb3-11eb-3975-0b419b7524a2
# ╠═27c72838-2eb3-11eb-29f4-ad9b1d5ce93e
# ╠═52196222-2eb3-11eb-3141-41c53f8cae56
# ╠═71b88180-2eb3-11eb-3937-8d8f6447d7b4
# ╠═0bfe2786-2eb9-11eb-3248-335c1647c09f
