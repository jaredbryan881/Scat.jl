@testset "scatter1d assertions" begin
N=10
x=ones(2^N)
Q=[12, 1]
J=[8, 8]
σ0=0.1
FilterBanks = get_FilterBanks(N, Q, J, σ0)
Transform = ScatteringTransform1d(3, Q, J)
@test_throws AssertionError scat=scatter1d(x, Transform, FilterBanks)
S=Scattered1d(FilterBanks)
FilterBanks_changed = get_FilterBanks(N, [8,4], J, σ0)
@test_throws AssertionError Scat.scatter1d_layer(S, 1, Vector{Complex{Float64}}(undef,10), Transform, FilterBanks_changed)
end
