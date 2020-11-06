# GeoStatsTutorials

[![gitter](https://img.shields.io/badge/chat-on%20gitter-bc0067.svg)](https://gitter.im/JuliaEarth/GeoStats.jl)

Tutorials in the form of Pluto notebooks for the [GeoStats.jl](https://github.com/juliohm/GeoStats.jl) framework.

## Notebooks

- Estimation problems
- Variogram modeling
- Anisotropic models
- Directional variograms
- Declustered statistics
- Two-point statistics
- Gaussian processes
- Image quilting
- Cookie-cutter

They can be run locally with [Pluto](https://github.com/fonsp/Pluto.jl):

Run Julia and add the package:
```
julia> ]
pkg> add Pluto
```

To run the notebook server:
```
using Pluto

cd("notebooks")

Pluto.run()
```
Pluto will open in your browser, and you can get started!

## Contributing

Contributions are very welcome, please submit a pull request or open an issue with an example that you feel is missing.
