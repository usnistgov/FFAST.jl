using Documenter
using Gadfly
using FFAST
using Weave

function addNISTHeaders(htmlfile::String)
    # read HTML
    html = transcode(String,read(htmlfile))
    # Find </head>
    i = findfirst(r"</[Hh][Ee][Aa][Dd]>", html)
    # Already added???
    j = findfirst("nist-header-footer", html)
    if isnothing(j) && (!isnothing(i))
        # Insert nist-pages links right before </head>
        res = html[1:i.start-1]*
            "<link rel=\"stylesheet\" href=\"https://pages.nist.gov/nist-header-footer/css/nist-combined.css\">\n"*
            "<script src=\"https://pages.nist.gov/nist-header-footer/js/jquery-1.9.0.min.js\" type=\"text/javascript\" defer=\"defer\"></script>\n"*
            "<script src=\"https://pages.nist.gov/nist-header-footer/js/nist-header-footer.js\" type=\"text/javascript\" defer=\"defer\"></script>\n"*
            html[i.start:end]
        write(htmlfile, res)
        println("Inserting NIST header/footer into $htmlfile")
    end
    return htmlfile
end

weaveit(name) = weave(joinpath("src", "$name"), out_path=joinpath("src", "$(splitext(name)[1]).md"), doctype="github")

names = ( "example.jmd", )

weaveit.(names)

makedocs(
    modules = [FFAST],
    sitename = "FFAST.jl",
    pages = [ "Home" => "index.md", "Using FFAST" => "example.md" ]
)

map(name->rm(joinpath("src","$(splitext(name)[1]).md")), names)
addNISTHeaders(joinpath(@__DIR__, "build","index.html"))
addNISTHeaders.(map(name->joinpath(@__DIR__, "build", splitext(name)[1], "index.html"), names))

# deploydocs(repo = "github.com/NicholasWMRitchie/FFAST.jl.git")
