@testset "Scattered1d definition" begin
# define filterbanks
N=12
Q=[12,1]
J=[8,4]
σ0=[0.5, 0.1]
FilterBanks = get_FilterBanks(N,Q,J,σ0)
# test inner constructor
S = Scattered1d(FilterBanks)
@test S.T==2^N
@test S.nFilters==[1, length(FilterBanks[1].Λ), length(FilterBanks[2].Λ)]
@test (size(S.Coeff[2])==(2^N, length(FilterBanks[1].Λ))) &
      (size(S.Coeff[3])==(2^N, length(FilterBanks[2].Λ)))
end
