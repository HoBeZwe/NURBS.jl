using NURBS
using Documenter

DocMeta.setdocmeta!(NURBS, :DocTestSetup, :(using NURBS); recursive=true)

makedocs(;
    modules=[NURBS],
    authors="Bernd Hofmann <bernd.hofmann@tum.de> and contributors",
    repo="https://github.com/HoBeZwe/NURBS.jl/blob/{commit}{path}#{line}",
    sitename="NURBS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://HoBeZwe.github.io/NURBS.jl", edit_link="main", assets=String[]
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/HoBeZwe/NURBS.jl", devbranch="main")
