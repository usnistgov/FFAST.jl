# FFAST.jl

This is a Julia 1.X implementation of the [NIST X-Ray Form Factor, Attenuation, and Scattering Tables (FFAST) database](https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables).
It is designed to
  * only expose data from the NIST FFAST database,
  * not depend upon any other physics-related modules,
  * not pollute the namespace, and
  * impose minimal requirements on the library user.
    * Elements are identified by atomic number (Int)
    * Shells/edges are identified by index (Int) where 1->K, 2->L1, 3->L2, ...., 25->O9

If you use this library please be sure to adequately reference:
  * Chantler, C.T., Olsen, K., Dragoset, R.A., Chang, J., Kishore, A.R., Kotochigova, S.A., and Zucker, D.S. (2005), *X-Ray Form Factor, Attenuation and Scattering Tables (version 2.1)*. [Online](http://physics.nist.gov/ffast) Downloaded: 10-Sep-2019. National Institute of Standards and Technology, Gaithersburg, MD.
  * Originally published as
    * [Chantler, C.T., J. Phys. Chem. Ref. Data 29(4), 597-1048 (2000)](https://physics.nist.gov/PhysRefData/FFast/Text2000/contents2000.html); and
    * [Chantler, C.T., J. Phys. Chem. Ref. Data 24, 71-643 (1995)](https://physics.nist.gov/PhysRefData/FFast/Text1995/contents1995.html).

  * This library is based on tables scraped from the NIST website on 10-Sep-2019.
  * The library implements functions to log-log interpolate the tabulated values for
    * form factors (f₁ and f₂)
    * mass attenuation/absorption coefficients μ/ρ(photoelectric), μ/ρ(coherent/incoherent), μ/ρ(total), μ/ρ(K-only)

  * In addition, the library implements functions to access the supplemental values
    * mean atomic weight
    * edge energies
    * nominal density
    * relativistic corrections
    * nuclear Thompson correction

```@autodocs
Modules = [FFAST]
```
