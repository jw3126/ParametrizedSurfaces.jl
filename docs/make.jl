using ParametrizedSurfaces
using Documenter

DocMeta.setdocmeta!(ParametrizedSurfaces, :DocTestSetup, :(using ParametrizedSurfaces); recursive=true)

makedocs(;
    modules=[ParametrizedSurfaces],
    authors="Jan Weidner <jw3126@gmail.com> and contributors",
    repo="https://github.com/jw3126/ParametrizedSurfaces.jl/blob/{commit}{path}#{line}",
    sitename="ParametrizedSurfaces.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jw3126.github.io/ParametrizedSurfaces.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jw3126/ParametrizedSurfaces.jl",
    devbranch="main",
)
