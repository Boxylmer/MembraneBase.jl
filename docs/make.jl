using Documenter, MembraneBase

DocMeta.setdocmeta!(MembraneBase, :DocTestSetup, :(using MembraneBase); recursive=true)

makedocs(;
    modules=[MembraneBase],
    authors="Will <william.joseph.box@gmail.com> and contributors",
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

ENV["GITHUB_REPOSITORY"] = "Boxylmer/MembraneBaseDocs"
deploydocs(;
    repo="github.com/Boxylmer/MembraneBaseDocs.git",
    devbranch="master",
)
