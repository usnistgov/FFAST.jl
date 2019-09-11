

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
    delta = ( energy/ffastEdgeEnergy(z,sh) for sh in 1:ffastEdgeCount(z))
    if (delta[1]>1.0) && (delta[1]<1.1)
        low, high = max(low, 0.1), max(high, 0.1)
    elseif (delta[1]>=1.1) && (delta[1]<1.2)
        low, high = max(low, 0.03), max(high, 0.03)
    end
    for sh in filter(sh -> ffastEdgeAvailable(z,sh) && (delta[sh]>=1.0), ( 2, 4, 6, 7))
        if delta[sh]<1.15
            low, high = max(low, 0.15), max(high, 0.15)
        elseif delta[sh]<1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    for sh in filter(sh -> ffastEdgeAvailable(z,sh) && (delta[sh]>=1.0), (3, 4, 8, 9))
        if delta[sh]<0.15
            low, high = max(low, 0.20), max(high, 0.20)
        elseif mlm1<0.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    if energy>200.0e3
        low, high = max(low, 0.02), max(high, 0.03)
    end
    return (low, high)
end

function ffastUncertaintiesSolidLiquid(z::Integer, energy)
    low, high = 0.01, 0.01
    if energy<200.0
        low, high = 1.0, 2.0
    elseif energy<500.0
        low, high = 0.50, 1.0
    elseif energy<1.0
        low, high = 0.05, 0.20
    end
    m = 1.0/max(0.00001, maximum( ( energy/(energy-ffastEdgeEnergy(z,sh)) for sh in 1:ffastEdgeCount(z)) ))
    if m<0.001
        low, high = max(low, 0.5), max(high, 0.5)
    end
    mk = (energy-ffastEdgeEnergy(z,1))/energy
    if (mk>1.0) && (mk<1.1)
        low, high = max(low, 0.1), max(high, 0.2)
    elseif (mk>=1.1) && (mk<1.2)
        low, high = max(low, 0.03), max(high, 0.03)
    end
    mlm1 = 1.0 / max(0.00001, maximum( ( energy/(energy-ffastEdgeEnergy(z,sh)) for sh in (2, 5, 6, 7) ) ))
    if mlm1<0.15
        low, high = max(low, 0.15), max(high, 0.30)
    elseif mlm1<0.4
        low, high = max(low, 0.04), max(high, 0.04)
    end
    mlm2 = 1.0 / max(0.00001, maximum( ( energy/(energy-ffastEdgeEnergy(z,sh)) for sh in (3, 4, 8, 9) ) ))
    if mlm1<0.15
        low, high = max(low, 0.20), max(high, 0.20)
    elseif mlm1<0.4
        low, high = max(low, 0.04), max(high, 0.04)
    end
    if energy>200.0e3
        low, high = max(low, 0.02), max(high, 0.03)
    end
    return (low, high)
end
