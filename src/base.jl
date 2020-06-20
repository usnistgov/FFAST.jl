using DataFrames
using CSV


struct FFASTElement
    energy::Vector{Float64} # in eV
    data::Dict{Symbol, Float64}
    edgeEnergies::Dict{Int, Float64}
    macs::DataFrame # Form factor and MAC data
    function FFASTElement(ev, sup, mac)
        ffastShells = ( :K, :L1, :L2, :L3, :M1, :M2, :M3, :M4, :M5, :N1, :N2, :N3, :N4,
                        :N5, :N6, :N7, :O1, :O2, :O3, :O4, :O5, :P1, :P2, :P3 )
        allshells = ( :K, :L1, :L2, :L3, :M1, :M2, :M3, :M4, :M5, :N1, :N2, :N3, :N4,
                      :N5, :N6, :N7, :O1, :O2, :O3, :O4, :O5, :O6, :O7, :O8, :O9, :P1,
                      :P2, :P3, :P4, :P5, :P6, :P7, :P8, :P9, :P10, :P11 )
        data = Dict( ( sh, convert(Float64, sup[1,sh]) ) for sh in ( :A, :xsec, :density, :ev, :rc1, :rc2, :nt ) )
        ee = Dict{Int, Float64}()
        for sh in ffastShells
            if !ismissing(sup[1,sh])
                ee[findfirst(s->s==sh,allshells)] = convert(Float64, sup[1,sh])
            end
        end
        return new(ev, data, ee, mac)
    end
end


function loadFFAST()::Vector{FFASTElement}
    path = dirname(pathof(@__MODULE__))
    res = Vector{FFASTElement}()
    data = CSV.read("$(path)\\..\\data\\shelldata.csv")
    for z in 1:92
        # println("Loading z = $(z)")
        macs = CSV.read("$(path)\\..\\data\\mac[$(z)].csv")
        # mapcols(x->convert.(Float64,x), macs) # ensure Float64
        push!(res, FFASTElement(1000.0 * macs[!,1], data[[z],:], macs))
    end
    return res
end

FFASTData = loadFFAST()

"""
    ffastElementRange()

The atomic number of the last element supported by the FFAST database.
"""
ffastElementRange() =
    1:length(FFASTData)


"""
    ffastEdges(z::Int)

Returns a set containing the shells for which there is an edge energy in the database
for the specified element.
"""
ffastEdges(z::Int) =
    keys(FFASTData[z].edgeEnergies)

"""
    ffastEdgeAvailable(z::Int, shell::Int)

Is a value available for the specific shell's edge energy for the element identified by atomic number, z.
"""
ffastEdgeAvailable(z::Int, shell::Int) =
    haskey(FFASTData[z].edgeEnergies, shell)

"""
    ffastEdgeEnergy(z::Int, shell::Int)

The edge energy (in eV) for the specific element and shell.  Chantler references these sources for the values
  1) Bearden, J.A., Rev. Mod. Phys. 39, 78-124 (1967).
  2) Bearden, J.A., Burr, A.F., Rev. Mod. Phys. 39, 125-142 (1967).
"""
ffastEdgeEnergy(z::Int, shell::Int) =
    FFASTData[z].edgeEnergies[shell]

"""
    ffastAtomicWeight(z::Int)

The mean atomic weight for the specified element.
"""
ffastAtomicWeight(z::Int) =
    FFASTData[z].data[:A]

"""
    ffastCrossSectionFactor(z::Int)

The constant factor to convert [μ/ρ] to cross section in cm²/atom.
"""
ffastCrossSectionFactor(z::Int) =
    FFASTData[z].data[:xsec] * 1.0e-24

"""
    ffastDensity(z::Int)

Nominal value of the density of the element.
"""
ffastDensity(z::Int) =
    FFASTData[z].data[:density]

"""
    ffastEV(z::Int)

E(eV) [μ/ρ](cm²/g) = f2(e/atom)  ×  ffastEV(z)

Nominal value of the density of the element.
"""
ffastEV(z::Int) =
    FFASTData[z].data[:ev]

"""
    ffastRelativisticCorrections(z::Int)

Relativistic correction estimates fᵣₑₗ(H82,3/5CL) in e/atom.
Returns a tuple with two values.
"""
ffastRelativisticCorrections(z::Int) =
    ( FFASTData[z].data[:rc1], FFASTData[z].data[:rc2] )

"""
    ffastNuclearThompsonCorrection(z::Int)

Nuclear Thomson correction fₙₜ</sub> in e/atom.
"""
ffastNuclearThompsonCorrection(z::Int) =
    FFASTData[z].data[:nt]

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
    @assert energy >= 0.0 "energy < 0.0 in ffastMAC"
    @assert energy <= FFASTData[z].energy[end] "Energy too large in ffastMAC"
    return binarySearch(FFASTData[z].energy, max(FFASTData[z].energy[1], energy))
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


"""
    ffastJumpRatio(z::Int, shell::Int)

Returns the jump-ratio for the specified element and shell.  This uses a simple
algorithm in which we look just above and just below the shell edge to
determine the height of the edge.  Returns zero if the edge isn't available.
"""
function ffastJumpRatio(z::Int, shell::Int)
    if ffastEdgeAvailable(z, shell)
        ee = ffastEdgeEnergy(z, shell)
        if ee > FFASTData[z].energy[1]
            if ffastEdgeAvailable(z, shell+1) && (ee == ffastEdgeEnergy(z, shell+1))
                # Clean up some ugliness!!!
                if shell == 3
                    return 1.41
                elseif shell == 4
                    return 1.16
                end
            end
            idx = ffastIndex(z, ee)
            ffd = FFASTData[z]
            return ffd.macs[idx+1,:μρpe] / ffd.macs[idx,:μρpe]
        end
    end
    return 0.0
end
