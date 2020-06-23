module FFAST

using Requires
using DataFrames
using CSV

# This type is used to identify the algorithm implementation to facilitate replacing one
# or more of the method implementations.
struct FFASTMAC end

# This type identifies the class of MAC to compute
abstract type MassAbsorptionCoefficientType end
# photoelectric mass absorption coefficient
struct PhotoElectricMAC <: MassAbsorptionCoefficientType end
# coherent/incoherent mass attenuation coefficient
struct CoherentIncoherentMAC <: MassAbsorptionCoefficientType end
# sum mass attenuation coefficient
struct TotalMAC <: MassAbsorptionCoefficientType end
# K-shell contribution to the total mass attenuation coefficient
struct KMAC <: MassAbsorptionCoefficientType end

# Used by the fractionaluncertainty function
abstract type MACUncertaintyType end
# The uncertainties associated with the MAC due to a solid or liquid
struct SolidLiquid <: MACUncertaintyType end
# The uncertainties associated with the MAC due to a monatomic gas (typically smaller)
struct MonatomicGas <: MACUncertaintyType end

include("base.jl")

export FFASTMAC
export eachelement
export eachedge
export hasedge
export edgeenergy
export atomicweight
export crosssectionfactor
export density
export EV
export relativisticcorrection
export nuclearthompsoncorrection
export formfactors # Tuple of two form factors in e/atom
export jumpratio # Ratio of the height after the edge over before the edge

export mac
export MassAbsorptionCoefficientType, PhotoElectricMAC, CoherentIncoherentMAC, TotalMAC, KMAC


include("uncertainties.jl")
export MACUncertaintyType, SolidLiquid, MonatomicGas
export fractionaluncertainty # Fractional uncertainties for a SolidLiquid or MonatomicGas

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflysupport.jl")
end

export plotMACs # Plot [μ/ρ].  Only available if Gadfly is loaded.
export plotFormFactors # Plot form factors.  Only available if Gadfly is loaded.
export plotEdgeCount # Plot edge/shell count.  Only available if Gadfly is loaded.
export plotAtomicWeight # Plot nominal atomic weight.  Only available if Gadfly is loaded.
export plotCrossSectionFactor # Plot crossection factor.  Only available if Gadfly is loaded.
export plotDensity # Plot nominal density.  Only available if Gadfly is loaded.
export plotRelativisticCorrection # Plot the relativistic corretion factors.  Only available if Gadfly is loaded.
export plotNuclearThompsonCorrection # Plot the nuclear Thompson correction.  Only available if Gadfly is loaded.
export plotJumpRatios # Plot the jump ratios for the specified shell

end
