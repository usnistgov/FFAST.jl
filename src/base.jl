using DataFrames
using CSV

struct FFASTElement
    energy::Vector{Float64} # in eV
    data::DataFrame
    macs::DataFrame
end

function loadFFAST()::Vector{FFASTElement}
    path = dirname(pathof(@__MODULE__))
    res=Vector{FFASTElement}()
    for z in 1:92
        println("Loading z = $(z)")
        data = CSV.read("$(path)\\..\\data\\data[$(z)].csv")
        macs = CSV.read("$(path)\\..\\data\\mac[$(z)].csv")
        push!(res,FFASTElement(1000.0*macs[!,1],data,macs))
    end
    return res
end

FFASTData = loadFFAST()
ffastShells = ( :K, :L1, :L2, :L3, :M1, :M2, :M3, :M4, :M5, :N1, :N2, :N3, :N4, :N5, :N6, :N7, :O1, :O2, :O3, :O4, :O5, :O6, :O7, :O8, :O9 )


ffastEdgeCount(zz::Int) =
    length(FFASTData[z].data)-7

ffastEdgeAvailable(zz::Int, shell::Int) =
    shell<=ffastEdgeCount(zz)

ffastEdgeEnergy(z::Int, shell::Int) =
    FFASTData[z].data[1,ffastShells[shell]]

ffastAtomicWeight(z::Int) =
    FFASTData[z].data[1,:A]

ffastCrossSectionFactor(z::Int) =
    FFASTData[z].data[1,:A]*10e-24

ffastDensity(z::Int) =
    FFASTData[z].data[1,:density]

ffastEV(z::Int) =
    FFASTData[z].data[1,:ev]

ffastRelativisticCorrections(z::Int) =
    ( FFASTData[z].data[1,:rc1], FFASTData[z].data[1,:rc2] )

ffastNuclearThompsonCorrection(z::Int) =
    FFASTData[z].data[1,:nt]

function binarySearch(array,value)
    start, stop = 1, length(array)
    while stop-start>1
        next=(start+stop)÷2
        if array[next]==value return next end
        (start,stop) = array[next]<value ? ( next,stop ) : ( start,next )
    end
    return start
end

function ffastIndex(z::Int, energy::Float64)
    df = FFASTData[z].macs
    @assert(energy>=1000.0*df[1,:E], "Energy too small in ffastMAC")
    @assert(energy<=1000.0*df[end,:E], "Energy too large in ffastMAC")
    return binarySearch(FFASTData[z].energy, energy)
end

linearInterp(x1, x2, y1, y2, x) =
    ((y2-y1)/(x2-x1))*(x-x1) + y1

loglogInterp(x0, x1, y0, y1, x) =
    exp(log(y1/y0)/log(x1/x0)*log(x/x0)+log(y0))

function ffastF(z::Int, energy::Float64)
    idx = ffastIndex(z,energy)
    ffd=FFASTData[z]
    return ( linearInterp(ffd.energy[idx],ffd.energy[idx+1],FFASTData[z].macs[idx,:f1],FFASTData[z].macs[idx+1,:f1],energy),
             linearInterp(ffd.energy[idx],ffd.energy[idx+1],FFASTData[z].macs[idx,:f2],FFASTData[z].macs[idx+1,:f2],energy) )
end

function ffastMACpe(z::Int, energy::Float64)
    idx = ffastIndex(z,energy)
    ffd=FFASTData[z]
    return linearInterp(ffd.energy[idx],ffd.energy[idx+1],FFASTData[z].macs[idx,:μρpe],FFASTData[z].macs[idx+1,:μρpe],energy)
end

function ffastMACci(z::Int, energy::Float64)
    idx = ffastIndex(z,energy)
    ffd=FFASTData[z]
    return linearInterp(ffd.energy[idx],ffd.energy[idx+1],FFASTData[z].macs[idx,:μρci],FFASTData[z].macs[idx+1,:μρci],energy)
end

function ffastMACtot(z::Int, energy::Float64)
    idx = ffastIndex(z,energy)
    ffd=FFASTData[z]
    return loglogInterp(ffd.energy[idx],ffd.energy[idx+1],FFASTData[z].macs[idx,:μρtot],FFASTData[z].macs[idx+1,:μρtot],energy)
end

function ffastMACμρK(z::Int, energy::Float64)
    idx = ffastIndex(z,energy)
    ffd=FFASTData[z]
    return linearInterp(ffd.energy[idx],ffd.energy[idx+1],FFASTData[z].macs[idx,:μρK],FFASTData[z].macs[idx+1,:μρK],energy)
end
