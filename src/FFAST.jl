module FFAST

using Requires
using DataFrames
using CSV

include("base.jl")

export ffastEdgeCount
export ffastEdgeAvailable
export ffastEdgeEnergy
export ffastAtomicWeight
export ffastCrossSectionFactor
export ffastDensity
export ffastEV
export ffastRelativisticCorrections
export ffastNuclearThompsonCorrection
export ffastFF # Tuple of two form factors in e/atom
export ffastMACpe # photoelectric mass absorption coefficient
export ffastMACci # coherent/incoherent mass attenuation coefficient
export ffastMACtot # sum mass attenuation coefficient
export ffastMACK # K-shell contribution to the total mass attenuation coefficient

include("uncertainties.jl")
export ffastUncertaintiesSolidLiquid # Fractional uncertainties for a solid or liquid
export ffastUncertaintiesMonatomicGas # Fractional uncertainties for a monatomic gas

function __init__()
    @require Gadfly="c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflysupport.jl")
end

export plotFfastMACs # Plot [μ/ρ].  Only available if Gadfly is loaded.
export plotFfastFF # Plot form factors.  Only available if Gadfly is loaded.
export plotFfastEdgeCount # Plot edge/shell count.  Only available if Gadfly is loaded.
export plotFfastAtomicWeight # Plot nominal atomic weight.  Only available if Gadfly is loaded.
export plotFfastCrossSectionFactor # Plot crossection factor.  Only available if Gadfly is loaded.
export plotFfastDensity # Plot nominal density.  Only available if Gadfly is loaded.
export plotFfastRelativisticCorrections # Plot the relativistic corretion factors.  Only available if Gadfly is loaded.
export plotFfastNuclearThompsonCorrection # Plot the nuclear Thompson correction.  Only available if Gadfly is loaded.

end
