# GeoStatsTutorials

[![gitter](https://img.shields.io/badge/chat-on%20gitter-bc0067.svg)](https://gitter.im/JuliaEarth/GeoStats.jl)

Tutorials in the form of [Pluto](https://github.com/fonsp/Pluto.jl) notebooks for the [GeoStats.jl](https://github.com/juliohm/GeoStats.jl) framework.

## Notebooks

- Estimation problems
- Variogram modeling
- Anisotropic models
- Directional variograms
- Variography game
- Declustered statistics
- Two-point statistics
- Gaussian processes
- Image quilting
- Cookie-cutter

To run the notebooks locally, install Pluto:

```julia
julia> ]
pkg> add Pluto
```

and launch it from the `notebooks` folder:

```julia
using Pluto

cd("notebooks")

Pluto.run()
```

## Contributing

Contributions are very welcome, please submit a pull request or open an issue with an example that you feel is missing.
