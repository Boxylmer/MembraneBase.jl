using Documenter, MembraneBase

DocMeta.setdocmeta!(MembraneBase, :DocTestSetup, :(using MembraneBase); recursive=true)

makedocs(;
    modules=[MembraneBase],
    authors="Will <sypheren@gmail.com> and contributors",
    repo="https://github.com/Boxylmer/MembraneBase.jl/blob/{commit}{path}#{line}",
    sitename="MembraneBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Boxylmer.github.io/MembraneBase.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Boxylmer/MembraneBase.jl.git",
    devbranch="master",
)
