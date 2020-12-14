export GaussianFilter, GaborWavelet, MorletWavelet

abstract type AbstractFilter end

"""
    gabor1d(ω, ξ, σ)

Calculate the Fourier transform of a gabor wavelet.

# Arguments
 - `ω::Vector{Float64}`: Frequencies over which to calculate the filter.
 - `σ::Float64`: Bandwidth of the filter
 - `ξ::Float64`: Central frequency
"""
function gabor1d(ω::Vector{Float64}, σ::Float64, ξ::Float64)
    return exp.(-(ω.-ξ).^2 ./ (2*σ^2))
end
"""
    gabor1d(N, ξ, σ)

Convenience method for calculating the Fourier transform of a gabor wavelet.

# Arguments
 - `ω::Vector{Float64}`: Frequencies over which to calculate the filter.
 - `σ::Float64`: Bandwidth of the filter
 - `ξ::Float64`: Central frequency
"""
function gabor1d(N::Int64, σ::Float64, ξ::Float64)
    ω = FFTW.fftfreq(N)
    return gabor1d(ω[:], σ, ξ)
end

"""
    gauss1d(ω, σ)

Return a Fourier transformed Gaussian lowpass filter.

# Arguments
 - `ω::Vector{Float64}`: Frequencies over which to calculate the filter.
 - `σ::Float64`: Bandwidth of lowpass filter.
"""
function gauss1d(ω::Vector{Float64}, σ::Float64)
    return exp.(-(ω.^2)./(2*σ^2))
end
"""
    gauss1d(N, σ)

Convenience method for calculate the Fourier transform of a Gaussian lowpass filter.

# Argument
 - `N::Int64`: Length of filter in the time domain.
 - `σ::Float64`: Bandwidth of the lowpass filter.
"""
function gauss1d(N::Int64, σ::Float64)
    ω = collect(0.0:N-1)./N
    return gauss1d(ω[:], σ)
end

"""
    morlet1d(gab, lowpass)

Calculate the Fourier transform of a Morlet wavelet.

# Arguments
 - `gab::Vector{Float64}`: Fourier transformed Gabor wavelet
 - `lowpass::Vector{Float64}`: Fourier transformed Gaussian lowpass filter
"""
function morlet1d(gab::Vector{Float64}, lowpass::Vector{Float64})
    # summation factor to ensure first value is 0
    κ = gab[1]/lowpass[1]
    return gab .- κ.*lowpass
end
"""
    morlet1d(N, ξ, σ)

Calculate the Fourier transform of a Morlet wavelet.

# Arguments
 - `N::Int64`: Length of filter in the time domain.
 - `σ::Float64`: Bandwidth of the lowpass filter.
 - `ξ::Float64`: Central frequency.
"""
function morlet1d(N::Int64, σ::Float64, ξ::Float64)
    gab = gabor1d(N, σ, ξ)
    lowpass = gauss1d(N, σ)
    return morlet1d(gab, lowpass)
end

"""
    GaussianFilter

A structure for Fourier transformed Gaussian lowpass Filters

## Fields: GaussianFilter
 | **Field** | **Description** |
 |:----------|:----------------|
 | :σ        | Filter Bandwidth |
 | :ϕ        | Fourier-transformed Gaussian lowpass filter |
"""
struct GaussianFilter <: AbstractFilter
    σ::Float64 # filter bandwidth

    ϕ::Vector{Float64} # Fourier-transformed lowpass Gaussian filter

    function GaussianFilter(N::Int64, σ::Float64)
        ϕ = gauss1d(N, σ)
        return new(σ, ϕ)
    end
end

"""
    GaborWavelet

A structure for Fourier transformed Gabor Wavelets

## Fields: GaborWavelet
 | **Field** | **Description** |
 |:----------|:----------------|
 | :σ        | Filter Bandwidth |
 | :ξ        | Central Frequency |
 | :j        | Maximal subsampling |
 | :γ        | Fourier-transformed Gabor wavelet |
"""
struct GaborWavelet <: AbstractFilter
    σ::Float64 # filter bandwidth
    ξ::Float64 # central frequency
    j::Int64 # 2^j is the maximal subsampling possible without aliasing

    γ::Vector{Float64} # Fourier-transformed Gabor wavelet

    function GaborWavelet(N::Int64, σ::Float64, ξ::Float64, j::Int64)
        γ = gabor1d(N, σ, ξ)
        return new(σ, ξ, j, γ)
    end
end

"""
    MorletWavelet

A structure for Fourier transformed Morlet Wavelets

## Fields: MorletWavelet
 | **Field** | **Description** |
 |:----------|:----------------|
 | :σ        | Filter Bandwidth |
 | :ξ        | Central Frequency |
 | :j        | Maximal subsampling |
 | :ψ        | Fourier-transformed Morlet wavelet |
"""
struct MorletWavelet <: AbstractFilter
    σ::Float64 # filter bandwidth
    ξ::Float64 # central frequency
    j::Int64 # 2^j is the maximal subsampling possible without aliasing

    ψ::Vector{Float64} # Fourier-transformed Morlet wavelet

    function MorletWavelet(N::Int64, σ::Float64, ξ::Float64, j::Int64)
        ψ = morlet1d(N, σ, ξ)
        return new(σ, ξ, j, ψ)
    end
    function MorletWavelet(γ::GaborWavelet, ϕ::GaussianFilter)
        ψ = morlet1d(γ.γ, ϕ.ϕ)
        return new(γ.σ, γ.ξ, γ.j, ψ)
    end
end
