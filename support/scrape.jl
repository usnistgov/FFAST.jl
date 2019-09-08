using HTTP
using Gumbo
using CSV
using DataFrames

function scrape(z::Int)
    sciNot = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
    numNot = r"[-+]?[0-9]*\.?[0-9]+"
    parseSN(dd) = parse(Float64,match(sciNot,dd).match)
    parseN(dd)=parse(Float64,match(numNot,dd).match)
    multiSN(dd) = [ parse(Float64, m.match) for m in eachmatch(sciNot, dd) ]
    html = HTTP.get("https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z=$(z)&Formula=&gtype=4&range=S&lower=0.001&upper=1000&density=&frames=no&htmltable=1");
    doc = parsehtml(String(html.body))
    headerNode=doc.root[2][1]
    if z<=2
        shells = multiSN(headerNode[4][1].text)[2:end]
    else
        shells = multiSN(headerNode[4][2][1].text)
    end
    data0 = zeros(Float64,1,7+length(shells))
    data0[1,1] = parseN(headerNode[1][1][3].text) # A '=\n  54.93800 g mol'
    if z==32
        data0[1,2] = parseSN(headerNode[2][7].text)
    else
        data0[1,2] = parseSN(headerNode[2][11].text) #  xsec ')  ×  9.12268E+01'
    end
    data0[1,3] = parseN(headerNode[1][6].text) # density ')  = 7.4300'
    if z==32
        data0[1,4] = parseSN(headerNode[3][13].text) # EV ')  ×  7.65960E+05'
    else
        data0[1,4] = parseSN(headerNode[3][17].text) # EV ')  ×  7.65960E+05'
    end
    if z<=2
        rc = multiSN(headerNode[5][4].text)
        data0[1,5], data0[1,6] = rc[length(rc)-1:end] # # Reativistic correction '(H82,3/5CL) = -9.8594E-02, -6.1200E-02'
        data0[7] = parseSN(headerNode[6][length(headerNode[6].children)-3].text)
    else
        rc = multiSN(headerNode[4][6].text)
        data0[1,5], data0[1,6] = rc[length(rc)-1:5] # # Relativistic correction '(H82,3/5CL) = -9.8594E-02, -6.1200E-02'
        data0[7] = parseSN(headerNode[5][length(headerNode[5].children)-3].text)
    end
    for i in eachindex(shells) data0[1,7+i] = 1000.0*shells[i] end # in eV
    shSym = [ :K, :L1, :L2, :L3, :M1, :M2, :M3, :M4, :M5, :N1, :N2, :N3, :N4, :N5, :N6, :N7, :O1, :O2, :O3, :O4, :O5, :O6, :O7, :O8, :O9   ]
    df0 = DataFrame(data0,[ :A, :xsec, :density, :ev, :rc1, :rc2, :nt, shSym[1:length(shells)]...])
    CSV.write("c:\\Users\\nritchie\\.julia\\dev\\FFAST\\data\\data[$(z)].csv", df0)

    tableNode=doc.root[2][8][1][1]
    data = zeros(Float64,length(tableNode.children)-1,8);
    for row in 2:length(tableNode.children)
        tr = tableNode[row]
        idx=1
        for col in 1:length(tr.children)
            nums = multiSN(tr[col][1].text)
            for num in nums
                data[row-1,idx]=num
                idx+=1
            end
        end
    end
    df1 = DataFrame(data)
    names!(df1,[:E,:f1,:f2,:μρpe,:μρci,:μρtot,:μρK,:λ])
    CSV.write("c:\\Users\\nritchie\\.julia\\dev\\FFAST\\data\\mac[$(z)].csv", df1)
end

for z in 1:92 scape(z) end
