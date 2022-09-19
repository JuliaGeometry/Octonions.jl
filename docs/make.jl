using Octonions
using Documenter

DocMeta.setdocmeta!(Octonions, :DocTestSetup, :(using Octonions); recursive=true)

makedocs(;
    modules=[Octonions],
    authors="Seth Axen <seth@sethaxen.com> and contributors",
    repo="https://github.com/JuliaGeometry/Octonions.jl/blob/{commit}{path}#{line}",
    sitename="Octonions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGeometry.github.io/Octonions.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGeometry/Octonions.jl",
    devbranch="main",
)
