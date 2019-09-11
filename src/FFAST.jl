module FFAST

using Requires
using DataFrames
using CSV

"""
# FFAST.jl

This is a Julia 1.X implementation of the NIST X-Ray Form Factor, Attenuation, and Scattering Tables (FFAST) database.
It is designed to
  * only expose data from the NIST FFAST database,
  * not depend upon any other physics-related modules,
  * not pollute the namespace, and
  * impose minimal requirements on the library user.
    * Elements are identified by atomic number (Int)
    * Shells/edges are identified by index (Int) where 1->K, 2->L1, 3->L2, ...., 25->O9

If you use this library please be sure to adequately reference:
  * Chantler, C.T., Olsen, K., Dragoset, R.A., Chang, J., Kishore, A.R., Kotochigova, S.A., and Zucker, D.S. (2005), *X-Ray Form Factor, Attenuation and Scattering Tables (version 2.1)*.
[Online](http://physics.nist.gov/ffast) Downloaded: 10-Sep-2019. National Institute of Standards and Technology, Gaithersburg, MD.
Originally published as [Chantler, C.T., J. Phys. Chem. Ref. Data 29(4), 597-1048 (2000);](https://physics.nist.gov/PhysRefData/FFast/Text2000/contents2000.html) and
[Chantler, C.T., J. Phys. Chem. Ref. Data 24, 71-643 (1995)](https://physics.nist.gov/PhysRefData/FFast/Text1995/contents1995.html).

  * This library is based on tables scraped from the NIST website on 10-Sep-2019.
  * The library implements functions to log-log interpolate the tabulated values for

    + form factors (f₁ and f₂)
    + mass attenuation/absorption coefficients μ/ρ(photoelectric), μ/ρ(coherent/incoherent), μ/ρ(total), μ/ρ(K-only)

  * In addition, the library implements functions to access the supplemental values

    + mean atomic weight
    + edge energies
    + nominal density
    + relativistic corrections
    + nuclear Thompson correction
"""

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
