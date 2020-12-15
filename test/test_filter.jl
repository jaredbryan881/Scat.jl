@testset "Gabor wavelet" begin
# there are two methods and one struct for calculating Gabor wavelets
# calculate Gabor wavelet given ω
N=4
ω=zeros(2^N)
σ=1.0
ξ=0.0
γ = Scat.gabor1d(ω, σ, ξ)
@test length(γ)==length(ω)
@test γ==ones(length(ω))

# calculate Gabor wavelet without ω
N=4
σ=1.0/sqrt(2.0)
ξ=0.0
ω=FFTW.fftfreq(N)
γ = Scat.gabor1d(N, σ, ξ)
@test length(γ)==length(ω)
@test exp.(-ω.^2)≈γ

# test inner constructor for GaborWavelet struct
N=4
σ=1.0
ξ=0.0
j=1
γ = GaborWavelet(N, σ, ξ, j)
@test length(γ.γ)==N
@test γ.γ == Scat.gabor1d(N, σ, ξ)
end

@testset "Gaussian filter" begin
# there are two methods and one struct for calculating Gaussian lowpass filters
# calculate Gaussian lowpass filter given ω
N=4
ω=zeros(2^N)
σ=1.0
G = Scat.gauss1d(ω, σ)
@test length(G)==length(ω)
@test G==ones(length(ω))

# calculate Gaussian lowpass filter without ω
N=4
σ=1.0/sqrt(2.0)
ω=collect(0.0:N-1)./N
G = Scat.gauss1d(N, σ)
@test length(G)==length(ω)
@test exp.(-ω.^2)≈G

# test inner constructor for GaussianFilter struct
N=4
σ=1.0
ξ=0.0
j=1
G = GaussianFilter(N, σ)
@test length(G.ϕ)==N
@test G.ϕ == Scat.gauss1d(N, σ)
end

@testset "Morlet wavelet" begin
# there are two methods and one struct for calculating Morlet wavelets
# calculate Morlet wavelet given Gabor wavelet and Gaussian lowpass filter
N=4
σ=1.0/sqrt(2.0)
ξ=0.0
G=ones(N)
γ=ones(N)
ψ=Scat.morlet1d(γ, G)
@test length(ψ)==N
@test ψ==zeros(N)

# calculate Morlet wavelet without ω
N=4
σ=1.0/sqrt(2.0)
ξ=0.0
ψ=Scat.morlet1d(N, σ, ξ)
@test length(G)==N

# test both inner constructors for MorletWavelet struct
N=4
σ=1.0
ξ=0.0
j=1
# get GaborWavelet and GaussianFilter structs to build MorletWavelet
G = GaussianFilter(N, σ)
γ = GaborWavelet(N, σ, ξ, j)
ψ = MorletWavelet(γ, G)
@test length(ψ.ψ)==N
@test ψ.ψ == Scat.morlet1d(γ.γ, G.ϕ)
end
