using Rigid
using Documenter

makedocs(;
    modules=[Rigid],
    authors="Hetao Z.",
    repo="https://github.com/HetaoZ/Rigid.jl/blob/{commit}{path}#L{line}",
    sitename="Rigid.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HetaoZ.github.io/Rigid.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HetaoZ/Rigid.jl",
)
