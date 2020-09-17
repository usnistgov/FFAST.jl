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

function loadFFAST()::NTuple{92,FFASTElement}
    path = dirname(pathof(@__MODULE__))
    res = Vector{FFASTElement}()
    shelldata = CSV.read(joinpath(path,"..","data","shelldata.csv"), DataFrame, header=1)
    for z in 1:92
        # println("Loading z = $(z)")
        macs = CSV.read(joinpath(path,"..","data","mac[$(z)].csv"), DataFrame, header=1)
        # mapcols(x->convert.(Float64,x), macs) # ensure Float64
        push!(res, FFASTElement(1000.0 * macs[!,1], shelldata[[z],:], macs))
    end
    return tuple(res...)
end

let FFASTData = loadFFAST()
    global getFFASTData() = FFASTData
end

"""
    eachelement(::Type{FFAST})

The range of available elements.
"""
eachelement() = eachindex(getFFASTData())


"""
    eachedge(z::Int)::Set{Integer}

Returns a set containing the shells for which there is an edge energy in the database
for the specified element.
"""
eachedge(z::Int) = keys(getFFASTData()[z].edgeEnergies)

"""
    hasedge(z::Int, shell::Int)::Bool

Is a value available for the specific shell's edge energy for the element identified by atomic number, z.
"""
hasedge(z::Int, shell::Int) = haskey(getFFASTData()[z].edgeEnergies, shell)

"""
    edgeenergy(z::Int, shell::Int)::Float64

The edge energy (in eV) for the specific element and shell.  Chantler references these sources for the values
  1) Bearden, J.A., Rev. Mod. Phys. 39, 78-124 (1967).
  2) Bearden, J.A., Burr, A.F., Rev. Mod. Phys. 39, 125-142 (1967).
"""
edgeenergy(z::Int, shell::Int) = getFFASTData()[z].edgeEnergies[shell]

"""
    atomicweight(z::Int)::Float64

The mean atomic weight for the specified element.
"""
atomicweight(z::Int) = getFFASTData()[z].data[:A]

"""
    crosssectionfactor(z::Int)::Float64

The constant factor to convert [μ/ρ] to cross section in cm²/atom.
"""
crosssectionfactor(z::Int) = getFFASTData()[z].data[:xsec] * 1.0e-24

"""
    density(z::Int)::Float64

Nominal value of the density of the element in g/cm³.
"""
density(z::Int) = getFFASTData()[z].data[:density]

"""
    EV(z::Int)::Float64

E(eV) [μ/ρ] in cm²/g = f₂(e/atom) ⋅ EV(z)

Conversion factor for the f₂ form factor.
"""
EV(z::Int) = getFFASTData()[z].data[:ev]

"""
    relativisticcorrection(z::Int)::Tuple{Float64,Float64}

Relativistic correction estimates fᵣₑₗ(H82,3/5CL) in e/atom.
Returns a tuple with two values.
"""
relativisticcorrection(z::Int) = ( getFFASTData()[z].data[:rc1], getFFASTData()[z].data[:rc2] )

"""
    nuclearthompsoncorrection(z::Int)::Float64

Nuclear Thomson correction fₙₜ</sub> in e/atom.
"""
nuclearthompsoncorrection(z::Int) = getFFASTData()[z].data[:nt]

function findindex(z::Int, energy::Float64)
    function binarysearch(array, value)
        start, stop = 1, length(array)
        while stop - start > 1
            next = (start + stop) ÷ 2
            if array[next] == value return next end
            (start, stop) = array[next] < value ? (next, stop) : (start, next)
        end
        return start
    end
    df = getFFASTData()[z].macs
    @assert energy >= 0.0 "energy < 0.0 in ffastMAC"
    @assert energy <= getFFASTData()[z].energy[end] "Energy too large in ffastMAC"
    return binarysearch(getFFASTData()[z].energy, max(getFFASTData()[z].energy[1], energy))
end

linearInterp(x1, x2, y1, y2, x) =
    ((y2 - y1) / (x2 - x1)) * (x - x1) + y1

loglogInterp(x0, x1, y0, y1, x) =
    exp(log(y1 / y0) / log(x1 / x0) * log(x / x0) + log(y0))

"""
    formfactors(z::Int, energy::Float64)::Tuple{Float64,Float64}

Returns a tuple containing the form factors for the specified element and X-ray energy (in eV).
"""
function formfactors(z::Int, energy::Float64)
    idx = findindex(z, energy)
    ffd = getFFASTData()[z]
    return (loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], getFFASTData()[z].macs[idx,:f1], getFFASTData()[z].macs[idx + 1,:f1], energy),
             loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], getFFASTData()[z].macs[idx,:f2], getFFASTData()[z].macs[idx + 1,:f2], energy))
end

"""
    mac(::Type{PhotoElectricMAC}, z::Int, energy::Float64)::Float64

Returns the photoelectric attenuation coefficient in cm²/g for the specified element and X-ray energy (in eV).
"""
function mac(::Type{PhotoElectricMAC}, z::Int, energy::Float64)
    idx = findindex(z, energy)
    ffd = getFFASTData()[z]
    return loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], getFFASTData()[z].macs[idx,:μρpe], getFFASTData()[z].macs[idx + 1,:μρpe], energy)
end

"""
    mac(::Type{CoherentIncoherentMAC}, z::Int, energy::Float64)::Float64

Returns the coherent/incoherent attenuation coefficient in cm²/g for the specified element and X-ray energy (in eV).
"""
function mac(::Type{CoherentIncoherentMAC}, z::Int, energy::Float64)
    idx = findindex(z, energy)
    ffd = getFFASTData()[z]
    return loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], getFFASTData()[z].macs[idx,:μρci], getFFASTData()[z].macs[idx + 1,:μρci], energy)
end

"""
    mac(::Type{TotalMAC}, z::Int, energy::Float64)::Float64

Returns the total attenuation coefficient in cm²/g for the specified element and X-ray energy (in eV).
"""
function mac(::Type{TotalMAC} ,z::Int, energy::Float64)
    idx = findindex(z, energy)
    ffd = getFFASTData()[z]
    return loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], getFFASTData()[z].macs[idx,:μρtot], getFFASTData()[z].macs[idx + 1,:μρtot], energy)
end

"""
    mac(::Type{KMAC}, z::Int, energy::Float64)::Float64

Returns the K-shell only attenuation coefficient in cm²/g for the specified element and X-ray energy (in eV).
"""
function mac(::Type{KMAC} ,z::Int, energy::Float64)
    idx = findindex(z, energy)
    ffd = getFFASTData()[z]
    return energy > edgeenergy(z,1) ?
        loglogInterp(ffd.energy[idx], ffd.energy[idx + 1], getFFASTData()[z].macs[idx,:μρK], getFFASTData()[z].macs[idx + 1,:μρK], energy) :
        0.0
end

"""
    jumpratio(z::Int, shell::Int)::Float64

Returns the jump-ratio for the specified element and shell.  This implementation attempts
to use the FFAST MAC data to extract the jump-ratio. It uses a simple
algorithm in which we look just above and just below the shell edge to
determine the height of the edge.  Returns zero if the edge isn't available.
"""
function jumpratio( z::Int, shell::Int)::Float64
    if hasedge(z, shell)
        ee = edgeenergy(z, shell)
        if ee > getFFASTData()[z].energy[1]
            if hasedge(z, shell+1) && (ee == edgeenergy(z, shell+1))
                # Clean up some ugliness!!!
                if shell == 3
                    return 1.41
                elseif shell == 4
                    return 1.16
                end
            end
            idx = findindex(z, ee)
            ffd = getFFASTData()[z]
            return ffd.macs[idx+1,:μρpe] / ffd.macs[idx,:μρpe]
        end
    end
    return 0.0
end
