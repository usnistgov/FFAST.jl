using Documenter
using Gadfly
using FFAST
using Weave

weave("src/example.jmd", out_path="src/example.html")

makedocs(modules = [FFAST], sitename = "FFAST.jl")

# deploydocs(repo = "github.com/NicholasWMRitchie/FFAST.jl.git")
