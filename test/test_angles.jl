# Pairwise angles between three unit vectors of the plane, at right angles and
# at 45 degrees.

@testset "angles                                 " begin
    Γ = couple([1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
               [1.0, 1.0, 0.0] ./ sqrt(2))

    @test angles(Γ) ≈ [0    π/2  π/4;
                       π/2  0    π/4;
                       π/4  π/4  0  ]  atol=1e-7
end
