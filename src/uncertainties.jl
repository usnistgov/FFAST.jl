
"""
    ffastUncertaintiesMonatomicGas(z::Integer, energy)

Determines from the element and energy, the approximate range of fractional uncertainties to
associate with the total and photoelectric components of the mass attenuation coefficients
for monatomic gas samples.
Based on [this table](https://physics.nist.gov/PhysRefData/FFast/Text2000/sec06.html#tab2).
"""
function ffastUncertaintiesMonatomicGas(z::Integer, energy)
    low, high = 0.01, 0.01
    if energy<200.0
        low, high = 0.5, 1.0
    elseif energy<500.0
        low, high = 0.20, 0.30
    elseif energy<1.0
        low, high = 0.03, 0.10
    end
    distance = ( (energy-ffastEdgeEnergy(z,sh)) / energy for sh in 1:ffastEdgeCount(z) )
    m = minimum.(abs.(distance))
    if m<0.001
        low, high = max(low, 0.2), max(high, 0.3)
    end
    u = ( energy/ffastEdgeEnergy(z, sh) for sh in 1:ffastEdgeCount(z) )
    if (u[1]>1.0) && (u[1]<1.1)
        low, high = max(low, 0.1), max(high, 0.1)
    elseif (u[1]>=1.1) && (u[1]<1.2)
        low, high = max(low, 0.03), max(high, 0.03)
    end
    for sh in filter(sh -> get(u, sh, 0.0)>=1.0, ( 2, 4, 6, 7))
        if u[sh]<1.15
            low, high = max(low, 0.15), max(high, 0.15)
        elseif u[sh]<1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    for sh in filter(sh -> get(u, sh, 0.0)>=1.0, (3, 4, 8, 9))
        if u[sh]<1.15
            low, high = max(low, 0.20), max(high, 0.20)
        elseif u[sh]<1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    if energy>200.0e3
        low, high = max(low, 0.02), max(high, 0.03)
    end
    return ( low, high )
end

"""
    ffastUncertaintiesMonatomicGas(z::Integer, energy)

Determines from the element and energy, the approximate range of fractional uncertainties to
associate with the total and photoelectric components of the mass attenuation coefficients
for solids and liquids.
Based on [this table](https://physics.nist.gov/PhysRefData/FFast/Text2000/sec06.html#tab2).
"""
function ffastUncertaintiesSolidLiquid(z::Integer, energy)
    low, high = 0.01, 0.01
    if energy<200.0
        low, high = 1.0, 2.0
    elseif energy<500.0
        low, high = 0.50, 1.0
    elseif energy<1.0
        low, high = 0.05, 0.20
    end
    distance = ( (energy-ffastEdgeEnergy(z,sh)) / energy for sh in 1:ffastEdgeCount(z) )
    m = minimum.(abs.(distance))
    if m<0.001
        low, high = max(low, 0.5), max(high, 0.5)
    end
    u = ( energy/ffastEdgeEnergy(z,sh) for sh in 1:ffastEdgeCount(z) )
    if (u[1]>1.0) && (u[1]<1.1)
        low, high = max(low, 0.1), max(high, 0.2)
    elseif (u[1]>=1.1) && (u[1]<1.2)
        low, high = max(low, 0.03), max(high, 0.03)
    end
    for sh in filter(sh -> get(u, sh, 0.0)>=1.0, ( 2, 4, 6, 7))
        if u[sh]<1.15
            low, high = max(low, 0.15), max(high, 0.30)
        elseif u[sh]<1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    for sh in filter(sh -> get(u, sh, 0.0)>=1.0, (3, 4, 8, 9))
        if u[sh]<1.15
            low, high = max(low, 0.20), max(high, 0.40)
        elseif u[sh]<1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    if energy>200.0e3
        low, high = max(low, 0.02), max(high, 0.03)
    end
    return (low, high)
end
