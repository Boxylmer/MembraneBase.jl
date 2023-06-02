using Documenter, MembraneBase

makedocs(sitename="MembraneBase.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ))


deploydocs(
    repo = "github.com/Boxylmer/MembraneBase.jl.git",
)