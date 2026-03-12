using Test

const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")

@time @testset verbose=true showtiming=true "TrixiManifolds.jl tests" begin
    @time if TRIXI_TEST == "all" || TRIXI_TEST == "wave_torus"
        include("test_wave_torus.jl")
    end
end
