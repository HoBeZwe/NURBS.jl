using Nurbs
using Documenter

DocMeta.setdocmeta!(Nurbs, :DocTestSetup, :(using Nurbs); recursive=true)

makedocs(;
    modules=[NURBS],
    authors="Bernd Hofmann <bernd.hofmann@tum.de> and contributors",
    sitename="Nurbs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://HoBeZwe.github.io/Nurbs.jl", edit_link="main", assets=String[]
    ),
    pages=[
        "Introduction" => "index.md",
        "Definitions" => Any["Bases" => "basis_def.md", "Curves" => "curves_def.md", "Surfaces" => "surfaces_def.md"],
        "Manual" => Any["Bases" => "basis.md", "Curves" => "curves.md", "Surfaces" => "surfaces.md"],
        "Utils" => "utils.md",
        "Contributing" => "contributing.md",
        "API Reference" => "apiref.md",
    ],
)

deploydocs(; repo="github.com/HoBeZwe/Nurbs.jl", devbranch="main")
