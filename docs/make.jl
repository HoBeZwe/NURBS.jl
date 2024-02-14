using NURBS
using Documenter

makedocs(;
    modules=[NURBS],
    authors="Bernd Hofmann <bernd.hofmann@tum.de> and contributors",
    sitename="NURBS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://HoBeZwe.github.io/NURBS.jl", edit_link="main", assets=String[]
    ),
    pages=[
        "Introduction" => "index.md",
        "Definitions" => Any["Bases" => "basis_def.md", "Curves" => "curves_def.md", "Surfaces" => "surfaces_def.md"],
        "Fundamental Operations" => Any[
            "Knot Insertion" => "knotInsertion.md",
            "Knot Removal" => "knotRemoval.md",
            "Degree Elevation" => "degreeElevation.md",
            "Degree Reduction" => "degreeReduction.md",
        ],
        "Manual" => Any[
            "Bases" => "basis.md",
            "Curves" => "curves.md",
            "Surfaces" => "surfaces.md",
            #"Fundamental Operations" => "fundOperations.md",
        ],
        "File I/O" => "fileio.md",
        "Utils" => "utils.md",
        "Contributing" => "contributing.md",
        "API Reference" => "apiref.md",
    ],
)

deploydocs(;
    repo="github.com/HoBeZwe/NURBS.jl",
    target="build",
    devbranch="main",
    push_preview=true,
    forcepush=true,
    versions=["stable" => "v^", "v#.#", "v0.5.1", "v0.5.0", "v0.4.1", "v0.4.0", "v0.3.0", "dev" => "dev"],
)
