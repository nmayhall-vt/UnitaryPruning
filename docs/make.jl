# Inside make.jl
push!(LOAD_PATH,"../src/")
using UnitaryPruning
using Documenter
pages = [
         "Home" => "index.md",
         "Functions" => Any[
                            "DFS" => "dfs.md",
                            "Table of Contents" => "toc.md"
                           ],
        ]
makedocs(
         sitename = "UnitaryPruning.jl",
         #source = "src",
         #build = "build",
         clean = true,
         modules  = [UnitaryPruning],
         authors="Nick Mayhall <nmayhall@vt.edu> and contributors",
         pages=pages)

deploydocs(
    repo="github.com/nmayhall-vt/UnitaryPruning.jl",
    branch = "gh-pages",
    devbranch = "main",
    target= "build",
    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math")
)
