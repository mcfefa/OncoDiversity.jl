using Documenter
# using Literate

#const literate_dir = joinpath(@__DIR__, "literate")
#const generated_dir = joinpath(@__DIR__, "src", "generated")

@info "Loading OncoDiversity.jl"
using OncoDiversity

# # XXX: Work around old LaTeX distribution in GitHub CI.
# @info "Building Literate.jl docs"
# if haskey(ENV, "GITHUB_ACTIONS")
#     import TikzPictures
#     TikzPictures.standaloneWorkaround(true)
#   end

# # Set Literate.jl config if not being compiled on recognized service.
# config = Dict{String,String}()
# if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
#   config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/mcfefa/OncoDiversity.jl/blob/gh-pages/dev"
#   config["repo_root_url"] = "https://github.com/mcfefa/OncoDiversity.jl/blob/master/docs"
# end

# for (root, dirs, files) in walkdir(literate_dir)
#   out_dir = joinpath(generated_dir, relpath(root, literate_dir))
#   for file in files
#     if last(splitext(file)) == ".jl"
#       Literate.markdown(joinpath(root, file), out_dir;
#         config=config, documenter=true, credit=false)
#       Literate.notebook(joinpath(root, file), out_dir;
#         execute=true, documenter=true, credit=false)
#     end
#   end
# end

@info "Building Documenter.jl docs"

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


@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/mcfefa/OncoDiversity.jl.git",
  branch = "gh-pages"
)
