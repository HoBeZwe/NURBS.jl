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
        "Definitions & Terminology" => Any[
            "B-Spline Bases" => "basis_def.md",
            "Manifolds" => Any["Curves" => "curves_def.md", "Surfaces" => "surfaces_def.md"],
            "Fundamental Operations" => Any[
                "Knot Insertion" => "knotInsertion.md",
                "Knot Removal" => "knotRemoval.md",
                "Degree Elevation" => "degreeElevation.md",
                "Degree Reduction" => "degreeReduction.md",
            ],
            "Connectivity" => Any["Patch Interfaces" => "interfaces.md", "Bezier Mesh" => "bezier.md"],
        ],
        "Manual" => Any[
            "Bases Evaluation" => "basis.md",
            "Curves" => Any["Defining & Evaluating" => "curves.md", "Manipulating" => "curves_man.md"],
            "Surfaces" => Any["Defining & Evaluating" => "surfaces.md", "Manipulating" => "surfaces_man.md"],
            "Connectivity" => "meshes.md",
            "Transformations" => "trafos.md",
            "File I/O" => "fileio.md",
            "Utils" => "utils.md",
        ],
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
