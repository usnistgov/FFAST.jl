using .Gadfly

function plotMACs(z::Integer)
    colors, name = [ "red", "green", "blue", "orange" ], [ "[μ/ρ]ₚₑ", "[μ/ρ]ci", "[μ/ρ]ₜₒₜ", "[μ/ρ]ₖ" ]
    lpe = layer(loge->log10(mac(PhotoElectricMAC, z, 10.0^loge)), log10(50.0), log10(4.0e5), Theme(default_color = colors[1]))
    lci = layer(loge->log10(mac(CoherentIncoherentMAC, z, 10.0^loge)), log10(50.0), log10(4.0e5), Theme(default_color = colors[2]))
    ltot = layer(loge->log10(mac(TotalMAC, z, 10.0^loge)), log10(50.0), log10(4.0e5), Theme(default_color = colors[3]))
    lK = layer(loge->log10(mac(KMAC,z, 10.0^loge)), log10(50.0), log10(4.0e5), Theme(default_color = colors[4]))
    plot(lpe, lci, ltot, lK,
        Guide.title("Z = $z"),
        Guide.manual_color_key("Type", name, colors),
        Guide.xlabel("log₁₀(E) (eV)"), Guide.ylabel("log₁₀([μ/ρ]) (g/cm²)"),
        Coord.cartesian(xmin = log10(50.0), xmax = log10(4.0e5), ymin = -2.0))
end

function plotFormFactors(z::Integer)
    colors, name = [ "red", "blue" ], [ "f₁(E)", "f₂(E)" ]
    lf1 = layer(loge->formfactors(z, 10.0^loge)[1], log10(50.0), log10(4.0e5), Theme(default_color = colors[1]))
    lf2 = layer(loge->formfactors(z, 10.0^loge)[2], log10(50.0), log10(4.0e5), Theme(default_color = colors[2]))
    plot(lf1, lf2,
        Guide.title("Z = $z"),
        Guide.manual_color_key("Type", name, colors),
        Guide.xlabel("log₁₀(E) (eV)"), Guide.ylabel("log₁₀(f(E)) (e/atom)"),
        Coord.cartesian(xmin = log10(50.0), xmax = log10(4.0e5), ymin = 0.0))
end


function plotEdgeCount()
    colors, name = [ "red", "blue" ], [ "Count", "Last" ]
    l1 = layer(x = eachelement(), y = length.(eachedge.(eachelement())), Theme(default_color = colors[1]))
    l2 = layer(x = eachelement(), y = maximum.(eachedge.(eachelement())), Theme(default_color = colors[2]))
    plot(l1, l2, Guide.title("Edges"), Guide.manual_color_key("Quantity", name, colors), Guide.xlabel("Z"), Guide.ylabel("Quantity"))
end

plotAtomicWeight() =
    plot(x = eachelement(), y = atomicweight.(eachelement()),
    Guide.title("Nominal Atomic Weight"),
    Guide.xlabel("Z"), Guide.ylabel("g/mole"))


plotCrossSectionFactor() =
    plot(x = eachelement(), y = 1.0e24 * crosssectionfactor.(eachelement()),
        Guide.xlabel("Z"), Guide.ylabel("barns/atom"),
        Guide.title("Cross Section Factor"),
        Coord.cartesian(xmin = 0, xmax = 92, ymin = 0.0, ymax = 400.0))

plotDensity() =
    plot(x = eachelement(), y = density.(eachelement()), #
    Guide.title("Nominal Density"), #
    Guide.xlabel("Z"), Guide.ylabel("g/cm³"))

function plotRelativisticCorrection()
    l1 = layer(x = eachelement(), y = map(z->relativisticcorrection(z)[1], eachelement()), Theme(default_color = "red"))
    l2 = layer(x = eachelement(), y = map(z->relativisticcorrection(z)[2], eachelement()), Theme(default_color = "blue"))
    plot(l1, l2, Guide.xlabel("Z"), Guide.ylabel("e/atom"), #
        Guide.manual_color_key("Type", ["1", "2" ], [ "red", "blue"]), #
        Guide.title("fᵣₑₗ(H82,3/5CL)"), #
        Coord.cartesian(xmin = 0, xmax = 92, ymin = -2.6, ymax = 0.0))
end


plotNuclearThompsonCorrection() =
    plot(x = eachelement(), y = nuclearthompsoncorrection.(eachelement()),
        Guide.title("Nuclear Thompson Correction"),
        Guide.xlabel("Z"), Guide.ylabel("e/atom"))


plotJumpRatios(shell::Int) =
    plot(x = eachelement(), y = jumpratio.(eachelement(), shell),
        Guide.title("Jump Ratio"),
        Guide.xlabel("Z"), Guide.ylabel("Ratio"))
