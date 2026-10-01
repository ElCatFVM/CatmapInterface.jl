using Documenter
using DocumenterMermaid
using ExampleJuggler, PlutoStaticHTML
using CairoMakie
using CatmapInterface
using Catalyst
using Pkg


# from https://github.com/Sienna-Platform/PowerSystems.jl/pull/1814
# Pin Mermaid to 11.16.1 (issue #1812). Mermaid 11.17.0 started bundling fastdom, which
# registers with Documenter's RequireJS instead of exporting itself: no diagram renders, and in
# Safari KaTeX breaks too. DocumenterMermaid hard-codes `mermaid@11`, so override its script.
# Remove this override once a Mermaid release includes the fix (mermaid-js/mermaid#8154).
function Documenter.HTMLWriter.domify(
        ::Documenter.HTMLWriter.DCtx,
        ::DocumenterMermaid.MarkdownAST.Node,
        ::DocumenterMermaid.MermaidScriptBlock,
    )
    Documenter.DOM.@tags script
    return script[:type => "module"](
        """
        import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@11.16.1/dist/mermaid.esm.min.mjs';
        mermaid.initialize({
            startOnLoad: true,
            theme: "neutral"
        });
        """
    )
end

function mkdocs()

    ExampleJuggler.verbose!(true)

    cleanexamples()
    notebookdir = joinpath(@__DIR__, "..", "notebooks")
    notebooks = [
        "CO2 reduction" => "CO2R.jl",
        "CO2R with raw API" => "CO2R_raw.jl",
    ]

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
