using HTTP
using Gumbo
using CSV
using DataFrames

sciNot = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
numNot = r"[-+]?[0-9]*\.?[0-9]+"
shellNot = r"[K-Q] (XII|XI|X|IX|VIII|VII|VI|V|IV|III|II|I)?"

shellmap = (
    "K ", "L I", "L II", "L III",
    "M I", "M II", "M III", "M IV", "M V",
    "N I", "N II", "N III", "N IV", "N V",
    "N VI","N VII", "O I", "O II", "O III",
    "O IV", "O V", "P I", "P II", "P III",
)


function parseEdges(dd)
    lines = split(dd,"\n")
    res = Vector{Any}(missing, length(shellmap))
    for line in lines
        if !isnothing(findfirst("24 edges", line))
            continue
        end
        shells, energies = [], []
        for (sh, num) in zip( eachmatch(shellNot, line), eachmatch(sciNot,line))
            # println("Shell = $(sh.match), Number = $(num.match)")
            res[findfirst(n->n==sh.match,shellmap)] = 1000.0*parse(Float64, num.match)
        end
        @assert length(shells) == length(energies)
    end
    return res
end

"""
    scrape(z::Int)
Scrapes FFAST data from the NIST website by downloading by element the
associated HTML page and parsing the content to extract the various different
values.  Not exported as part of FFAST but provided to assist with downloading
updated values if they should become available.
"""
function scrape(z::Int)
    parseSN(dd) = parse(Float64, match(sciNot, dd).match)
    parseN(dd) = parse(Float64, match(numNot, dd).match)
    multiSN(dd) = [ parse(Float64, m.match) for m in eachmatch(sciNot, dd) ]
    html = HTTP.get("https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z=$(z)&Formula=&gtype=4&range=S&lower=0.001&upper=1000&density=&frames=no&htmltable=1");
    doc = parsehtml(String(html.body))
    headerNode = doc.root[2][1]
    if z <= 2
        shells = Vector{Any}(missing, length(shellmap))
        shells[1]= 1000.0*multiSN(headerNode[4][1].text)[2]
    else
        shells = parseEdges(headerNode[4][2][1].text)
    end
    data0 = Array{Any}(missing, 1, 7 + length(shells))
    data0[8:end] = shells
    data0[1,1] = parseN(headerNode[1][1][3].text) # A '=\n  54.93800 g mol'
    if z == 32
        data0[1,2] = parseSN(headerNode[2][7].text)
    else
        data0[1,2] = parseSN(headerNode[2][11].text) #  xsec ')  ×  9.12268E+01'
    end
    data0[1,3] = parseSN(headerNode[1][6].text) # density ')  = 7.4300'
    if z == 32
        data0[1,4] = parseSN(headerNode[3][13].text) # EV ')  ×  7.65960E+05'
    else
        data0[1,4] = parseSN(headerNode[3][17].text) # EV ')  ×  7.65960E+05'
    end
    if z <= 2
        rc = multiSN(headerNode[5][4].text)
        data0[1, 5], data0[1,6] = rc[length(rc) - 1:end] # # Reativistic correction '(H82,3/5CL) = -9.8594E-02, -6.1200E-02'
        data0[1, 7] = parseSN(headerNode[6][length(headerNode[6].children) - 3].text)
    else
        rc = multiSN(headerNode[4][6].text)
        data0[1, 5], data0[1, 6] = rc[length(rc) - 1:5] # # Relativistic correction '(H82,3/5CL) = -9.8594E-02, -6.1200E-02'
        data0[1, 7] = parseSN(headerNode[5][length(headerNode[5].children) - 3].text)
    end
    shSym = [ :K, :L1, :L2, :L3, :M1, :M2, :M3, :M4, :M5,
        :N1, :N2, :N3, :N4, :N5, :N6, :N7,
        :O1, :O2, :O3, :O4, :O5, :P1, :P2, :P3   ]
    df0 = DataFrame(data0, [ :A, :xsec, :density, :ev, :rc1, :rc2, :nt, shSym...])
    CSV.write("c:\\Users\\nritchie\\.julia\\dev\\FFAST\\data\\data[$(z)].csv", df0)

    tableNode = doc.root[2][8][1][1]
    data = zeros(Float64, length(tableNode.children) - 1, 8);
    for row in 2:length(tableNode.children)
        tr = tableNode[row]
        idx = 1
        for col in 1:length(tr.children)
            nums = multiSN(tr[col][1].text)
            for num in nums
                data[row - 1,idx] = num
                idx += 1
            end
        end
    end
    df1 = DataFrame(data)
    names!(df1, [:E,:f1,:f2,:μρpe,:μρci,:μρtot,:μρK,:λ])
    CSV.write("c:\\Users\\nritchie\\.julia\\dev\\FFAST\\data\\mac[$(z)].csv", df1)
end

for z in 1:92 scrape(z) end

df0=DataFrame( A=[],xsec=[],density=[],ev=[],rc1=[],rc2=[],nt=[],K=[],
    L1=[], L2=[], L3=[], M1=[], M2=[], M3=[], M4=[], M5=[],
    N1=[], N2=[], N3=[], N4=[], N5=[], N6=[], N7=[],
    O1=[], O2=[], O3=[], O4=[], O5=[], P1=[], P2=[], P3=[]  )
all = mapreduce(z->z == 0 ? df0 : CSV.read("c:\\Users\\nritchie\\.julia\\dev\\FFAST\\data\\data[$(z)].csv"),append!,0:92)
CSV.write("c:\\Users\\nritchie\\.julia\\dev\\FFAST\\data\\shelldata.csv",all)
