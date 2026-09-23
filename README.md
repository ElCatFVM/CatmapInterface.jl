[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://elcatfvm.github.io/CatmapInterface.jl/dev/)
[![Build status](https://github.com/ElCatFVM/CatmapInterface.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ElCatFVM/CatmapInterface.jl/actions/workflows/CI.yml)

# CatmapInterface: [CatMap](https://github.com/SUNCAT-Center/catmap) based microkinetic models in Julia


CatmapInterface implements the functionality of a subset of the Python package [CatMAP](https://catmap.readthedocs.io) that is used in the Python package [CatINT](https://catint.readthedocs.io).
Whereas [CatINT](https://catint.readthedocs.io) uses an iterative approach for combining the Poisson-Nernst-Planck transport model with a microkinetic model of the surface reactions at the electrode, the output of CatmapInterface can be directly plugged into the functionality of [LiquidElectrolytes](https://j-fu.github.io/LiquidElectrolytes.jl) to solve the coupled system.

## Installation

### Version >=0.4
Starting with version 0.4, the package can be installed with the Julia package manager from the Julia General Registry

### Older versions
Package versions up to v0.3.1 are registered in the julia package registry [https://github.com/j-fu/PackageNursery](https://github.com/j-fu/PackageNursery)
To add the registry (needed only once), and to install the package, 
from the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> registry add https://github.com/j-fu/PackageNursery
pkg> add https://github.com/ElcatFVM/CatmapInterface.jl
```

Please be aware that adding a registry to your Julia installation requires to
trust the registry maintainer for handling things in a correct way. In particular,
the registry should not register higher versions of packages which are already
registered in the Julia General Registry. One can check this by visiting the above mentionend
github repository URL and inspecting the contents.

