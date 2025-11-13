using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ArgParse


function main(args)
    s = ArgParseSettings(description = "Create and deploy documentation")
    @add_arg_table! s begin
        "--makedocs"
        action = :store_true
        help = "Create documentation"

    end
    parsed_args = parse_args(s) # the result is a Dict{String,Any}
    makedocs = parsed_args["makedocs"]
    deploydir = "gate.wias-berlin.de:/server/www/ROOT/people/fuhrmann/_projects/CatmapInterface"
    docsdir = joinpath(@__DIR__, "..", "docs")
    if makedocs
        Pkg.activate(docsdir)
        Pkg.develop(path = joinpath(@__DIR__, ".."))
        Pkg.resolve()
        include(joinpath(docsdir, "make.jl"))
    end
    return run(`rsync -avu $(joinpath(docsdir,"build"))/ $(deploydir)`)
end

main(ARGS)
