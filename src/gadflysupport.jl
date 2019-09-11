using .Gadfly

function plotFfastMACs(z::Integer)
    colors = [ "red", "green", "blue", "orange" ]
    name = [ "[μ/ρ]pe", "[μ/ρ]ci", "[μ/ρ]tot", "[μ/ρ]K" ]
    lpe = Gadfly.layer(loge->log10(ffastMACpe(z,10.0^loge)), log10(50.0), log10(4.0e5),Gadfly.Theme(default_color=colors[1]))
    lci = Gadfly.layer(loge->log10(ffastMACci(z,10.0^loge)), log10(50.0), log10(4.0e5),Gadfly.Theme(default_color=colors[2]))
    ltot = Gadfly.layer(loge->log10(ffastMACtot(z,10.0^loge)), log10(50.0), log10(4.0e5),Gadfly.Theme(default_color=colors[3]))
    lK = Gadfly.layer(loge->log10(ffastMACK(z,10.0^loge)), log10(50.0), log10(4.0e5),Gadfly.Theme(default_color=colors[4]))
    Gadfly.plot(lpe, lci, ltot, lK,
        Gadfly.Guide.title("Mass Attenuation Coefficients"),
        Gadfly.Guide.manual_color_key("Type", name, colors),
        Gadfly.Guide.xlabel("log(E) (eV)"), Guide.ylabel("log([μ/ρ]) (g/cm²)"),
        Gadfly.Coord.cartesian(xmin=log10(50.0), xmax=log10(4.0e5), ymin=-2.0))
end

function plotFfastFF(z::Integer)
    colors = [ "red", "blue" ]
    name = [ "f₁(E)", "f₂(E)" ]
    lf1 = Gadfly.layer(loge->ffastFF(z,10.0^loge)[1], log10(50.0), log10(4.0e5),Gadfly.Theme(default_color=colors[1]))
    lf2 = Gadfly.layer(loge->ffastFF(z,10.0^loge)[2], log10(50.0), log10(4.0e5),Gadfly.Theme(default_color=colors[2]))
    Gadfly.plot(lf1, lf2,
        Gadfly.Guide.title("Form Factors"),
        Gadfly.Guide.manual_color_key("Type", name, colors),
        Gadfly.Guide.xlabel("log(E) (eV)"), Guide.ylabel("log(f(E)) (e/atom)"),
        Gadfly.Coord.cartesian(xmin=log10(50.0), xmax=log10(4.0e5), ymin=0.0))
end


plotFfastEdgeCount() =
    Gadfly.plot(x=1:92, y=ffastEdgeCount.(1:92),
    Gadfly.Guide.title("Edge/Shell Count"),
    Gadfly.Guide.xlabel("Z"), Guide.ylabel("Shell count") )


plotFfastAtomicWeight() =
    Gadfly.plot(x=1:92, y=ffastAtomicWeight.(1:92),
    Gadfly.Guide.title("Nominal Atomic Weight"),
    Gadfly.Guide.xlabel("Z"), Guide.ylabel("g/mole") )


plotFfastCrossSectionFactor() =
    Gadfly.plot(x=1:92, y=1.0e24*ffastCrossSectionFactor.(1:92),
        Gadfly.Guide.xlabel("Z"), Guide.ylabel("barns/atom"),
        Gadfly.Guide.title("Crosssection Factor"),
        Gadfly.Coord.cartesian(xmin=0, xmax=92, ymin=0.0, ymax=400.0) )

plotFfastDensity() =
    Gadfly.plot(x=1:92, y=ffastDensity.(1:92),
    Gadfly.Guide.title("Nominal Density"),
    Gadfly.Guide.xlabel("Z"), Guide.ylabel("g/cm³") )

function plotFfastRelativisticCorrections()
    l1 = Gadfly.layer(x=1:92, y=map(x->x[1],ffastRelativisticCorrections.(1:92)), Gadfly.Theme(default_color="red"))
    l2 = Gadfly.layer(x=1:92, y=map(x->x[2],ffastRelativisticCorrections.(1:92)), Gadfly.Theme(default_color="blue"))
    plot(l1, l2,
        Gadfly.Guide.xlabel("Z"), Guide.ylabel("e/atom"),
        Gadfly.Guide.manual_color_key("Type", ["1", "2" ], [ "red", "blue"]),
        Gadfly.Guide.title("fᵣₑₗ(H82,3/5CL)"),
        Gadfly.Coord.cartesian(xmin=0, xmax=92, ymin=-2.6, ymax=0.0) )
end


plotFfastNuclearThompsonCorrection() =
    Gadfly.plot(x=1:92, y=ffastNuclearThompsonCorrection.(1:92),
        Gadfly.Guide.title("Nuclear Thompson Correction"),
        Gadfly.Guide.xlabel("Z"), Guide.ylabel("e/atom") )
