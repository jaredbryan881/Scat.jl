@testset "Scattering Transform Definition" begin
Transform = ScatteringTransform1d(2, [12, 1], [8, 8])
@test Scat.depth(Transform)==2
end
