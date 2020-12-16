@testset "scatter1d assertions" begin
N=10
x=ones(2^N)
Q=[12, 1]
J=[8, 8]
σ0=0.1
# test difference in depth and number of filter banks for FilterBank1d
FilterBanks = get_FilterBanks(N, Q, J, σ0)
Transform = ScatteringTransform1d(3, Q, J)
@test_throws AssertionError scat=Scat.scatter1d(x, Transform, FilterBanks)

# test difference in depth and number of filter banks for FilterBank1dBlock
FilterBankBlocks=get_FilterBanks(N, Q, J, σ0, typ=FilterBank1dBlock)
@test_throws AssertionError scat=Scat.scatter1d(x, Transform, FilterBankBlocks)

# test transform depth
FilterBankBlocks=get_FilterBanks(N, [12,12,1], [8,8,8], σ0, typ=FilterBank1dBlock)
Transform = ScatteringTransform1d(3, Q, J)
@test_throws AssertionError scat=Scat.scatter1d(x, Transform, FilterBankBlocks)
end
