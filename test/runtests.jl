using Test, FFAST

@testset "Iron"
    @test ffastEdgeCount(26)==9

    @test ffastEdgeAvailable(26,1)
    @test ffastEdgeAvailable(26,2)
    @test ffastEdgeAvailable(26,3)
    @test ffastEdgeAvailable(26,4)
    @test ffastEdgeAvailable(26,5)
    @test ffastEdgeAvailable(26,6)
    @test ffastEdgeAvailable(26,7)
    @test ffastEdgeAvailable(26,8)
    @test ffastEdgeAvailable(26,9)
    @test !ffastEdgeAvailable(26,10)

    @test ffastEdgeEnergy(26,1)==7112.0
    @test ffastEdgeEnergy(26,9)==3.6
    @test ffastEdgeEnergy(26,4)≈708.1

    @test isapprox(ffastAtomicWeight(26),55.845,atol=0.003)

    @test ffastDensity(26)≈7.860

    @test ffastEV(26)≈7.53493E+05

    @test ffastCrossSectionFactor(26) ≈ 9.27362E+01*1.0e-24

    @test ffastRelativisticCorrections(26)[1] ≈ -1.0812E-01
    @test ffastRelativisticCorrections(26)[2] ≈ -6.7800E-02

    @test ffastNuclearThompsonCorrection(26) ≈ -6.6403E-03
end
