# test spectral parameters for the filters
@testset "Spectral filter parameters" begin
# test ξmax
Q=12 # number of wavelets per octave
ξ=ξmax(Q)
@test ξ==1/(1+2^(0.25))
Q=1 # number of wavelets per octave
ξ=ξmax(Q)
@test ξ==0.35

# test σmin
J=1
σ0=0.1
σ = σmin(σ0, J)
@test σ==σ0/2.0

# test σmax
ξ = 1.0
rψ = 2.718281828459045^(-1/2)
Q = 1
σ = σmax(ξ, Q, rψ)
@test σ==1/3

# test dyadic subsampling
α=1.0
ξ=1.0
σ=1.0
mds = max_dyadic_subsampling(ξ, σ, α)
@test mds==0
α=0.0
ξ=0.25
σ=1.0
mds = max_dyadic_subsampling(ξ, σ, α)
@test mds==1

σ_min = 0.5
σ_max = 2.0
Q = 1
nds = n_dyadic_steps(Q, σ_min, σ_max)
@test nds==2

# test filter bank parameters
# single Q
Q=1
σ_min=0.05
rψ=sqrt(0.5)
α=5.0
σs, ξs, js = get_filter_params(Q, σ_min, rψ=rψ, α=α)
ξ_max = ξmax(Q)
σ_max = σmax(ξ_max, Q, rψ)
nds = n_dyadic_steps(Q, σ_min, σ_max)
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
