using DataFrames
using CSV

struct FFASTElement
    energy::Vector{Float64} # in eV
    data::DataFrame # Supplemental data
    macs::DataFrame # Form factor and MAC data
    function FFASTElement(ev, sup, mac)
        @assert all(ty->isequal(ty, Float64), eltypes(sup))
        @assert all(ty->isequal(ty, Float64), eltypes(mac))
        return new(ev, sup, mac)
    end
end


function loadFFAST()::Vector{FFASTElement}
    path = dirname(pathof(@__MODULE__))
    res = Vector{FFASTElement}()
    for z in 1:92
        # println("Loading z = $(z)")
        data = CSV.read("$(path)\\..\\data\\data[$(z)].csv")
        data = mapcols(x->convert.(Float64, x), data) # ensure Float64
        macs = CSV.read("$(path)\\..\\data\\mac[$(z)].csv")
        # mapcols(x->convert.(Float64,x), macs) # ensure Float64
        push!(res, FFASTElement(1000.0 * macs[!,1], data, macs))
    end
    return res
end

FFASTData = loadFFAST()
ffastShells = (:K, :L1, :L2, :L3, :M1, :M2, :M3, :M4, :M5, :N1, :N2, :N3, :N4, :N5, :N6, :N7, :O1, :O2, :O3, :O4, :O5, :O6, :O7, :O8, :O9)

"""
    ffastElementCount()

The atomic number of the last element supported by the FFAST database.
"""
ffastElementCount() =
    length(FFASTData)


"""
    ffastEdgeCount(z::Int)

Returns the number of edge energy values available for the element identified by atomic number, z.
"""
ffastEdgeCount(z::Int) =
    length(FFASTData[z].data[1,:]) - 7

"""
    ffastEdgeAvailable(z::Int, shell::Int)

Is a value available for the specific shell's edge energy for the element identified by atomic number, z.
"""
ffastEdgeAvailable(z::Int, shell::Int) =
    shell <= ffastEdgeCount(z)

"""
    ffastEdgeEnergy(z::Int, shell::Int)

The edge energy (in eV) for the specific element and shell.  Chantler references these sources for the values
  1) Bearden, J.A., Rev. Mod. Phys. 39, 78-124 (1967).
  2) Bearden, J.A., Burr, A.F., Rev. Mod. Phys. 39, 125-142 (1967).
"""
ffastEdgeEnergy(z::Int, shell::Int) =
    FFASTData[z].data[1,ffastShells[shell]]

"""
    ffastAtomicWeight(z::Int)

The mean atomic weight for the specified element.
"""
ffastAtomicWeight(z::Int) =
    FFASTData[z].data[1,:A]

"""
    ffastCrossSectionFactor(z::Int)

The constant factor to convert [μ/ρ] to cross section in cm²/atom.
"""
ffastCrossSectionFactor(z::Int) =
    FFASTData[z].data[1,:xsec] * 1.0e-24

"""
    ffastDensity(z::Int)

Nominal value of the density of the element.
"""
ffastDensity(z::Int) =
    FFASTData[z].data[1,:density]

"""
    ffastEV(z::Int)

E(eV) [μ/ρ](cm²/g) = f2(e/atom)  ×  ffastEV(z)

Nominal value of the density of the element.
"""
ffastEV(z::Int) =
    FFASTData[z].data[1,:ev]

"""
    ffastRelativisticCorrections(z::Int)

Relativistic correction estimates fᵣₑₗ(H82,3/5CL) in e/atom.
Returns a tuple with two values.
"""
ffastRelativisticCorrections(z::Int) =
    (FFASTData[z].data[1,:rc1], FFASTData[z].data[1,:rc2])

"""
    ffastNuclearThompsonCorrection(z::Int)

Nuclear Thomson correction fₙₜ</sub> in e/atom.
"""
ffastNuclearThompsonCorrection(z::Int) =
    FFASTData[z].data[1,:nt]

function binarySearch(array, value)
    start, stop = 1, length(array)
    while stop - start > 1
        next = (start + stop) ÷ 2
        if array[next] == value return next end
        (start, stop) = array[next] < value ? (next, stop) : (start, next)
    end
    return start
end

function ffastIndex(z::Int, energy::Float64)
    df = FFASTData[z].macs
    @assert(energy >= 1000.0 * df[1,:E], "Energy too small in ffastMAC")
    @assert(energy <= 1000.0 * df[end,:E], "Energy too large in ffastMAC")
    return binarySearch(FFASTData[z].energy, energy)
end

linearInterp(x1, x2, y1, y2, x) =
    ((y2 - y1) / (x2 - x1)) * (x - x1) + y1

loglogInterp(x0, x1, y0, y1, x) =
    exp(log(y1 / y0) / log(x1 / x0) * log(x / x0) + log(y0))

"""
    ffastF(z::Int, energy::Float64)

Returns a tuple containing the form factors for the specified element and X-ray energy (in eV).
"""
function ffastFF(z::Int, energy::Float64)
    idx = ffastIndex(z, energy)
    ffd = FFASTData[z]
    return (loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], FFASTData[z].macs[idx,:f1], FFASTData[z].macs[idx + 1,:f1], energy),
             loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], FFASTData[z].macs[idx,:f2], FFASTData[z].macs[idx + 1,:f2], energy))
end

"""
    ffastMACpe(z::Int, energy::Float64)

Returns the photoelectric attenuation coefficient in cm²/g for the specified element and X-ray energy (in eV).
"""
function ffastMACpe(z::Int, energy::Float64)
    idx = ffastIndex(z, energy)
    ffd = FFASTData[z]
    return loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], FFASTData[z].macs[idx,:μρpe], FFASTData[z].macs[idx + 1,:μρpe], energy)
end

"""
    ffastMACci(z::Int, energy::Float64)

Returns the coherent/incoherent attenuation coefficient in cm²/g for the specified element and X-ray energy (in eV).
"""
function ffastMACci(z::Int, energy::Float64)
    idx = ffastIndex(z, energy)
    ffd = FFASTData[z]
    return loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], FFASTData[z].macs[idx,:μρci], FFASTData[z].macs[idx + 1,:μρci], energy)
end

"""
    ffastMACtot(z::Int, energy::Float64)

Returns the total attenuation coefficient in cm²/g for the specified element and X-ray energy (in eV).
"""
function ffastMACtot(z::Int, energy::Float64)
    idx = ffastIndex(z, energy)
    ffd = FFASTData[z]
    return loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], FFASTData[z].macs[idx,:μρtot], FFASTData[z].macs[idx + 1,:μρtot], energy)
end

"""
    ffastMACK(z::Int, energy::Float64)

Returns the K-shell only attenuation coefficient in cm²/g for the specified element and X-ray energy (in eV).
"""
function ffastMACK(z::Int, energy::Float64)
    idx = ffastIndex(z, energy)
    ffd = FFASTData[z]
    return loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], FFASTData[z].macs[idx,:μρK], FFASTData[z].macs[idx + 1,:μρK], energy)
end
