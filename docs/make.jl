using Documenter
using DocumenterMermaid
using ExampleJuggler, PlutoStaticHTML
using CairoMakie
using CatmapInterface
using Catalyst
using Pkg


function mkdocs()

    ExampleJuggler.verbose!(true)

    cleanexamples()
    notebookdir = joinpath(@__DIR__, "..", "notebooks")
    notebooks = ["CO2 reduction" => "CO2R.jl"]
    notebook_examples = @docplutonotebooks(notebookdir, notebooks, iframe = false, append_build_context = false)

    DocMeta.setdocmeta!(CatmapInterface, :DocTestSetup, :(using CatmapInterface, Catalyst); recursive = true)

    return makedocs(
        sitename = "CatmapInterface.jl",
        modules = [CatmapInterface],
        format = Documenter.HTML(
            size_threshold = nothing,
            mathengine = MathJax3(
                Dict(
                    :tex => Dict(
                        "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                        "tags" => "ams",
                        "packages" => ["base", "ams", "autoload", "mhchem"],
                    )
                )
            )
        ),
        clean = false,
        doctest = false,
        draft = false,
        authors = "S. Maaß",
        repo = "https://github.com/ElCatFVM/CatmapInterface.jl",
        pages = [
            "Home" => "index.md",
            "Guide Outline" => "guide.md",
            "Public API" => "public.md",
            "Information flow" => "package-flow.md",
            "Internal API" => "internal.md",
            "Notebooks" => notebook_examples,
        ]
    )
end


mkdocs()

if !isinteractive()
    deploydocs(repo = "github.com/ElCatFVM/CatmapInterface.jl.git", devbranch = "main")
end
