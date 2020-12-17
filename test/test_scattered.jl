@testset "Scattered1d definition" begin
# define filterbanks
N=12
Q=[12,1]
J=[8,8]
σ0=[0.1, 0.1]

# test construction with FilterBank1d
FilterBanks = get_FilterBanks(N,Q,J,σ0)
# test inner constructor
S = Scattered1d(FilterBanks)

@test S.T==2^N
@test S.nFilters==[1, length(FilterBanks[1].Λ), length(FilterBanks[2].Λ)]

# test size of second layer
j1s = [ψ1.j for ψ1 in FilterBanks[1].Λ]
j2s = [ψ2.j for ψ2 in FilterBanks[2].Λ]
nNodes=sum([length([j2 for j2 in j2s if j2>j1]) for j1 in j1s])
@test nNodes==size(S.Coeff[3],2)

# test construction with FilterBank1dBlock
FilterBanks = get_FilterBanks(N,Q,J,σ0,typ=FilterBank1dBlock)
# test inner constructor
S = Scattered1d(FilterBanks)

@test S.T==2^N
@test S.nFilters==[1, size(FilterBanks[1].Λ,2), size(FilterBanks[2].Λ,2)]

# test size of second layer
nNodes=sum([length([j2 for j2 in FilterBanks[2].js if j2>j1]) for j1 in FilterBanks[1].js])
@test nNodes==size(S.Coeff[3],2)
end
