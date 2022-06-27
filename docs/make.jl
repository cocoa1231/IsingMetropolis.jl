using IsingMetropolis
using Documenter

DocMeta.setdocmeta!(IsingMetropolis, :DocTestSetup, :(using IsingMetropolis); recursive=true)

makedocs(;
    modules=[IsingMetropolis],
    authors="Cocoa <cocoathepenguin@protonmail.com> and contributors",
    repo="https://github.com/cocoa1231/IsingMetropolis.jl/blob/{commit}{path}#{line}",
    sitename="IsingMetropolis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cocoa1231.github.io/IsingMetropolis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cocoa1231/IsingMetropolis.jl",
    devbranch="main",
)
