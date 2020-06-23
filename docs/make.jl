using Documenter
using Gadfly
using FFAST
using Weave

weaveit(name) = weave(joinpath("src", "$name"), out_path=joinpath("src", "$(splitext(name)[1]).md"), doctype="github")

names = ( "example.jmd", )

weaveit.(names)

makedocs(
    modules = [FFAST],
    sitename = "FFAST.jl",
    pages = [ "Home" => "index.md", "Using FFAST" => "example.md" ]
)

map(name->rm(joinpath("src","$(splitext(name)[1]).md")), names)

# deploydocs(repo = "github.com/NicholasWMRitchie/FFAST.jl.git")
