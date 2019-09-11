using Documenter
using Gadfly
using FFAST

makedocs(modules = [FFAST], sitename = "FFAST.jl")

deploydocs(repo = "github.com/NicholasWMRitchie/FFAST.jl.git")
