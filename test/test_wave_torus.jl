module TestWaveTorus

include("test_triximanifolds.jl")

@trixi_testset "elixir_wave_torus_short" begin
    @test_trixi_include(joinpath(@__DIR__, "elixir_wave_torus_short.jl"),
                        l2=[0.005452328126366846, 3.4027041766292076e-18,
                            0.02855329879820269],
                        linf=[0.009383782249574741, 1.879879915184092e-17,
                            0.050748360894784],
                        atol=5.0e-4)
end

end
