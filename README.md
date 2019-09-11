# FFAST.jl
## X-ray Mass Attenuation and Absorption Factors

Implements in Julia an interpolation scheme for Chantler's FFAST X-Ray Form Factor, Attenuation, and Scattering Tables as downloaded from the [NIST FFAST web site](https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables)

Users should refer to the NIST website for a detailed discussion of the provenance of these values. In general, however, they are based on a self-consistent Dirac-Hartree-Fock framework.  Interpolated values of the f<sub>1</sub> and f<sub>2</sub> components of the form factors, the photoelectric attenuation coefficient for the atom, [µ/ρ]<sub>pe</sub>, the coherent and incoherent attenuation coefficient, [μ/ρ]<sub>ci</sub>, the total mass absorption coefficient [μ/ρ]<sub>tot</sub> and the value for the K-shell, [µ/ρ]<sub>K</sub> are made available through exported functions.

In addition, exported functions provide access to other elemental data items including the mean atomic weight, the cross-section factors, the nominal density, E(V), the relativistic correction factors and the nuclear Thompson correction.

The tabulations are discussed in two references:
* Theoretical Form Factor, Attenuation and Scattering Tabulation for Z=1-92 from E=1-10 eV to E=0.4-1.0 MeV
J. Phys. Chem. Ref. Data 1995 (https://physics.nist.gov/PhysRefData/FFast/Text1995/contents1995.html)
* Detailed Tabulation of Atomic Form Factors, Photoelectric Absorption and Scattering Cross Section, and Mass Attenuation Coefficients in the Vicinity of Absorption Edges in the Soft X-Ray (Z = 30-36, Z = 60-89, E = 0.1 keV-10 keV), Addressing Convergence Issues of Earlier Work
J. Phys. Chem. Ref. Data 2000 (https://physics.nist.gov/PhysRefData/FFast/Text2000/contents2000.html)
