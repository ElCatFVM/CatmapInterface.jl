# CatmapInterface

| **Documentation**                 | **Build Status**        |
|:---------------------------------:|:-----------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][GHA-img]][GHA-url] |

## Installation

The package can be installed with the Julia package manager.
For the time being, it is registered in the julia package registry [https://github.com/j-fu/PackageNursery](https://github.com/j-fu/PackageNursery)
To add the registry (needed only once), and to install the package, 
from the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> registry add https://github.com/j-fu/PackageNursery
pkg> add https://github.com/smaasz/CatmapInterface.jl
```

Please be aware that adding a registry to your Julia installation requires to
trust the registry maintainer for handling things in a correct way. In particular,
the registry should not register higher versions of packages which are already
registered in the Julia General Registry. One can check this by visiting the above mentionend
github repository URL and inspecting the contents.


## Documentation

- [**DEVEL**][docs-dev-url] &mdash; *documentation of the in-development version.*


[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://smaasz.github.io/CatmapInterface.jl/dev/

[GHA-img]: https://github.com/smaasz/CatmapInterface.jl/actions/workflows/CI.yml/badge.svg?branch=main
[GHA-url]: https://github.com/smaasz/CatmapInterface.jl/actions/workflows/CI.yml?query=branch%3Amain

