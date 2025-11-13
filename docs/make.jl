using Documenter
using CatmapInterface
using Catalyst
using PlutoStaticHTML
using Pkg


function mkdocs()
    # TODO: add notebooks with the help of ExampleJuggler 

    DocMeta.setdocmeta!(CatmapInterface, :DocTestSetup, :(using CatmapInterface, Catalyst); recursive=true)

    makedocs(
        sitename    = "CatmapInterface.jl",
        modules     = [CatmapInterface],
        format      = Documenter.HTML(
            size_threshold  = nothing,
            mathengine      = MathJax3(
                Dict(
                    :tex => Dict(
                        "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
                        "tags" => "ams",
                        "packages" => ["base", "ams", "autoload", "mhchem"],
                    )
                )
            )
        ),
        clean       = false,
        doctest     = false,
        draft       = false,
        authors     = "S. Maaß",
        repo        = "https://github.com/ElCatFVM/CatmapInterface.jl",
        pages       = [
            "Home"      => "index.md",
            "Guide"     => "guide.md",
            "Public"    => "public.md",
            "Internal"  => "internal.md",
            #"Notebooks" => notebooks,
        ]
    )
end

mkdocs()

if !isinteractive()
#    deploydocs(repo = "github.com/ElCatFVM/CatmapInterface.jl.git", devbranch = "main")
end
