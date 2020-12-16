# test spectral parameters for the filters
@testset "Spectral filter parameters" begin
# test ξmax
Q=12 # number of wavelets per octave
ξ=Scat.ξmax(Q)
@test ξ==1/(1+2^(0.25))
Q=1 # number of wavelets per octave
ξ=Scat.ξmax(Q)
@test ξ==0.35

# test σmin
J=1
σ0=0.1
σ = Scat.σmin(σ0, J)
@test σ==σ0/2.0

# test σmax
ξ = 1.0
rψ = 2.718281828459045^(-1/2)
Q = 1
σ = Scat.σmax(ξ, Q, rψ)
@test σ==1/3

# test dyadic subsampling
α=1.0
ξ=1.0
σ=1.0
mds = Scat.max_dyadic_subsampling(ξ, σ, α)
@test mds==0
α=0.0
ξ=0.25
σ=1.0
mds = Scat.max_dyadic_subsampling(ξ, σ, α)
@test mds==1

σ_min = 0.5
σ_max = 2.0
Q = 1
nds = Scat.n_dyadic_steps(Q, σ_min, σ_max)
@test nds==2

# test filter bank parameters
# single Q
Q=1
σ_min=0.05
rψ=sqrt(0.5)
α=5.0
σs, ξs, js = get_filter_params(Q, σ_min, rψ=rψ, α=α)
ξ_max = Scat.ξmax(Q)
σ_max = Scat.σmax(ξ_max, Q, rψ)
nds = Scat.n_dyadic_steps(Q, σ_min, σ_max)
collect(0.0:nds+Q-2)
@test length(σs)==nds+Q-1
@test length(ξs)==nds+Q-1
@test length(js)==nds+Q-1
# array of Q
Q=[12, 1]
J=12
σ0=0.1
rψ=sqrt(0.5)
α=5.0
σs, ξs, js = get_filter_params(Q, J, σ0, rψ=rψ, α=α)
@test length(σs)==length(Q)
@test length(ξs)==length(Q)
@test length(js)==length(Q)
end

# test FilterBank1d struct
@testset "Morlet wavelet filter bank" begin
# test inner constructor
N=12
Q=12
J=8
σ0=0.1
rψ=sqrt(0.5)
α=5.0

F = FilterBank1d(N, Q, J, σ0, rψ=rψ, α=α)
# test frequency vector
ω=FFTW.fftfreq(2^N)
@test ω==F.ω
# test number of wavelets in the filter bank
σs, ξs, js = get_filter_params(Q, σ0/2^J, rψ=rψ, α=α)
@test length(F.Λ)==length(σs)
end

# test get_FilterBanks
@testset "Convenience methods for FilterBank1d creation" begin
N=12
Q=[12,1]
J=[8,4]
σ0=[0.5, 0.1]
FilterBanks = get_FilterBanks(N,Q,J,σ0)
@test length(FilterBanks)==length(Q)
@test (FilterBanks[1].Q==Q[1]) & (FilterBanks[2].Q==Q[2])
@test (FilterBanks[1].J==J[1]) & (FilterBanks[2].J==J[2])

# common J and σ0
N=12
Q=[12,1]
J=8
σ0=0.1
FilterBanks = get_FilterBanks(N,Q,J,σ0)
@test length(FilterBanks)==length(Q)
@test (FilterBanks[1].Q==Q[1]) & (FilterBanks[2].Q==Q[2])
@test (FilterBanks[1].J==J) & (FilterBanks[2].J==J)

# common J
N=12
Q=[12,1]
J=8
σ0=[0.5, 0.1]
FilterBanks = get_FilterBanks(N,Q,J,σ0)
@test length(FilterBanks)==length(Q)
@test (FilterBanks[1].Q==Q[1]) & (FilterBanks[2].Q==Q[2])
@test (FilterBanks[1].J==J) & (FilterBanks[2].J==J)

# common σ0
N=12
Q=[12,1]
J=[8, 4]
σ0=0.1
FilterBanks = get_FilterBanks(N,Q,J,σ0)
@test length(FilterBanks)==length(Q)
@test (FilterBanks[1].Q==Q[1]) & (FilterBanks[2].Q==Q[2])
@test (FilterBanks[1].J==J[1]) & (FilterBanks[2].J==J[2])
end

# test FilterBank1dBlock struct
@testset "Morlet wavelet filter bank block" begin
# test inner constructor
N=12
Q=12
J=8
σ0=0.1
F = FilterBank1d(N, Q, J, σ0)
Fblock = FilterBank1dBlock(N, Q, J, σ0)
# TODO: Find a better way to test the equality of the FilterBank1d and FilterBank1dBlock fields
testBoolΛ=true
testBoolσ=true
testBoolξ=true
testBoolj=true
for i=1:length(F.Λ)
    if F.Λ[i].ψ != Fblock.Λ[:,i]
        testBoolΛ=false
    end
    if F.Λ[i].σ != Fblock.σs[i]
        testBoolσ=false
    end
    if F.Λ[i].ξ != Fblock.ξs[i]
        testBoolξ=false
    end
    if F.Λ[i].j != Fblock.js[i]
        testBoolj=false
    end
end
# test whether morlet wavelets are the same
@test testBoolΛ=true
@test testBoolσ=true
@test testBoolξ=true
@test testBoolj=true
# test whether gaussian filter is the same
@test F.ϕ.ϕ==Fblock.ϕ
@test F.ϕ.σ==Fblock.σg
end
