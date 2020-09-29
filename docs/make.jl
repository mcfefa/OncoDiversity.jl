using OncoDiversity
using Documenter

makedocs(;
    modules=[OncoDiversity],
    authors="Meghan Ferrall-Fairbanks <meghan.ferrall.fairbanks@gmail.com> and contributors",
    repo="https://github.com/mcfefa/OncoDiversity.jl/blob/{commit}{path}#L{line}",
    sitename="OncoDiversity.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mcfefa.github.io/OncoDiversity.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mcfefa/OncoDiversity.jl",
)
