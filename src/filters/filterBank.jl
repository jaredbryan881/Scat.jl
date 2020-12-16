export get_filter_params, FilterBank1d, get_FilterBanks, FilterBank1dBlock

"""
    ξmax(Q)

Calculate the maximum central frequency of a Morlet wavelet.

# Arguments
 - `Q::Int64`: Number of wavelets per octave.
"""
function ξmax(Q::Int64)
    return max(1/(1+2^(3/Q)), 0.35)
end

"""
    σmin(σ0, J)

Calculate the bandwidth of the lowpass filter ϕ.

# Arguments
 - `σ0::Float64`: Bandwidth of lowpass filter ϕ at scale 0.
 - `J::Int64`: Maximum wavelet scale.
"""
function σmin(σ0::Float64, J::Int64)
    return σ0/2^J
end

"""
    σmax(ξmax, Q, rψ)

Calculate the bandwidth of a Morlet wavelet in a family of Q per octave.

# Arguments
 - `ξmax::Float64`: Maximum central frequency.
 - `Q::Int64`: Number of wavelets per octave.
 - `rψ::Float64`: Parameter controlling frequency width of filters.
"""
function σmax(ξmax::Float64, Q::Int64, rψ::Float64)
    return ξmax * ((1-2^(-1/Q))/(1+2^(-1/Q))) / (2*log(1/rψ))^(1/2)
end

"""
    max_dyadic_subsampling(ξ, σ, α)

Compute the maximal subsampling possible for a Gabor wavelet without aliasing.

# Arguments
 - `ξ::Float64`: Central frequency of the wavelet
 - `σ::Float64`: Bandwidth of the wavelet
 - `α::Float64`: Parameter controlling allowable error, with error∝1/α
"""
function max_dyadic_subsampling(ξ::Float64, σ::Float64, α::Float64)
    upper_bound = min(ξ+α*σ, 0.5)
    j = convert(Int, floor(-log2(upper_bound))) - 1
    return j
end

"""
    n_dyadic_steps(Q, σ_min, σ_max)

Get the number of dyadic steps between σ_min and σ_max given Q wavelets per octave.

# Arguments
 - `Q::Int64`: Number of wavelets per octave.
 - `σ_min::Float64`: Bandwidth of the lowpass filter
 - `σ_max::Float64`: Maximum filter bandwidth
"""
function n_dyadic_steps(Q::Int64, σ_min::Float64, σ_max::Float64)
    return convert(Int, ceil(Q*log2(σ_max/σ_min)))
end

"""
    get_filter_params(Q, J, σ0, rψ=rψ, α=α)

Get the central frequencies, bandwidths, and maximum subsampling for a family of Morlet wavelets.

# Arguments
 - `Q::Int64`: Number of wavelets per octave.
 - `J::Int64`: Maximum scale of the filters.
 - `σ_min::Float64`: Bandwidth of the low-pass filter.
 - `rψ`::Float64: Parameter controlling frequency width of filters.
 - `α::Float64`: Parameter controlling allowable error, with error∝1/α
"""
function get_filter_params(Q::Int64, σ_min::Float64; rψ::Float64=sqrt(0.5), α::Float64=5.0)
    @assert Q>=1 "At least one wavelet per octave is required."

    # maximum central frequency
    ξ_max = ξmax(Q)

    # maximum filter bandwidth
    σ_max = σmax(ξ_max, Q, rψ)

    # number of dyadic steps between σ_min and σ_max
    n = n_dyadic_steps(Q, σ_min, σ_max)

    # initialize vectors
    σ = collect(0.0:n+Q-2)
    ξ = collect(0.0:n+Q-2)
    j = Vector{Int64}(undef, n+Q-1)

    # σ_i = σ_max/2^(i/Q)
    σ[1:n] .= σ_max./2.0.^(σ[1:n]./Q)
    σ[n+1:end] .= σ_min
    # ξ_i = ξ_max/2^(i/Q)
    ξ[1:n] .= ξ_max./2.0.^(ξ[1:n]./Q)
    ξ[n+1:end] .= ξ[n].*((Q.-collect(1:Q-1))./Q)

    # js
    j .= max_dyadic_subsampling.(ξ, σ, α)

    return σ, ξ, j
end

"""
    get_filter_params(Q, J, σ0, rψ=rψ, α=α)

Get the central frequencies, bandwidths, and maximum subsampling for families of
Morlet wavelets for each layer of the scattering transform.

# Arguments
- `Q::Array{Int64}`: Number of wavelets per octave for each layer.
- `J::Int64`: Maximum scale of the filters.
- `σ0::Float64`: Bandwidth of the low-pass filter.
- `rψ`::Float64: Parameter controlling frequency width of filters.
- `α::Float64`: Parameter controlling allowable error, with error∝1/α
"""
function get_filter_params(Q::Vector{Int64}, J::Int64, σ0::Float64; rψ::Float64=sqrt(0.5), α::Float64=5.0)
    # minimum filter bandwidth
    σ_min = σ0/2^J

    # initialize arrays for the parameters for each layer
    ξs = Vector{Vector{Float64}}(undef, length(Q))
    σs = Vector{Vector{Float64}}(undef, length(Q))
    js = Vector{Vector{Float64}}(undef, length(Q))

    # get filter parameters for each layer
    for i=1:length(Q)
        σs[i], ξs[i], js[i] = get_filter_params(Q[i], σ_min, rψ=rψ, α=α)
    end

    return σs, ξs, js
end


abstract type AbstractFilterBank end

"""
    FilterBank1d

A structure for collections of 1D Morlet Wavelets and a Gaussian lowpass filter

## Fields: FilterBank1d
 | **Field** | **Description** |
 |:----------|:----------------|
 | :N        | log2(signal length) |
 | :J        | log2(max filter scale) |
 | :Q        | Number of wavelets per octave |
 | :ω        | Frequencies |
 | :ϕ        | Fourier transformed low-pass Gaussian filter |
 | :Λ        | Fourier transformed Morlet wavelets |
"""
struct FilterBank1d <: AbstractFilterBank
    N::Int # log2 length of the signal in the time domain
    J::Int # log2 maximum scale of filters
    Q::Int # Number of wavelets per octave

    ω::Vector{Float64} # frequencies
    ϕ::GaussianFilter # low-pass filter
    Λ::Vector{MorletWavelet} # Fourier-transformed Morlet wavelets stored in a vector

    function FilterBank1d(N::Int64, Q::Int64, J::Int64, σ0::Float64; rψ::Float64=sqrt(0.5), α::Float64=5.0)
        σs, ξs, js = get_filter_params(Q, σ0/2^J, rψ=rψ, α=α)
        ω = FFTW.fftfreq(2^N)

        Λ = Vector{MorletWavelet}(undef, length(σs))
        ϕ = GaussianFilter(2^N, σ0/2^J)
        for i=1:length(σs)
            γ = GaborWavelet(2^N, σs[i], ξs[i], js[i])
            Λ[i] = MorletWavelet(γ, ϕ)
        end

        return new(N, J, Q, ω, ϕ, Λ)
    end
end

"""
    get_FilterBanks(N, Q, J, σ0, rψ=rψ, α=α)

Get filter banks for every layer.

# Arguments
- `N::Int64`: log2(signal length)
- `Q::Vector{Int64}`: Number of wavelets per octave for each layer.
- `J::Vector{Int64}`: Maximum scale of the filters for each layer
- `σ0::Vector{Float64}`: Bandwidth of the low-pass filter for each layer
- `rψ`::Float64: Parameter controlling frequency width of filters.
- `α::Float64`: Parameter controlling allowable error, with error∝1/α
"""
function get_FilterBanks(N::Int64, Q::Vector{Int64}, J::Vector{Int64}, σ0::Vector{Float64}; rψ::Float64=sqrt(0.5), α::Float64=5.0)
    FilterBanks=Vector{FilterBank1d}(undef, length(Q))
    for i=1:length(Q)
        FilterBanks[i] = FilterBank1d(N, Q[i], J[i], σ0[i], rψ=rψ, α=α)
    end
    return FilterBanks
end
"""
    get_FilterBanks(N, Q, J, σ0, rψ=rψ, α=α)

Convenience method fo getting filter banks for each layer with common σ0.
"""
function get_FilterBanks(N::Int64, Q::Vector{Int64}, J::Vector{Int64}, σ0::Float64; rψ::Float64=sqrt(0.5), α::Float64=5.0)
    σ0s = repeat([σ0], length(Q))
    get_FilterBanks(N, Q, J, σ0s; rψ=rψ, α=α)
end
"""
    get_FilterBanks(N, Q, J, σ0, rψ=rψ, α=α)

Convenience method fo getting filter banks for each layer with common J.
"""
function get_FilterBanks(N::Int64, Q::Vector{Int64}, J::Int64, σ0::Vector{Float64}; rψ::Float64=sqrt(0.5), α::Float64=5.0)
    Js = repeat([J], length(Q))
    get_FilterBanks(N, Q, Js, σ0; rψ=rψ, α=α)
end
"""
    get_FilterBanks(N, Q, J, σ0, rψ=rψ, α=α)

Convenience method fo getting filter banks for each layer with common σ0 and J.
"""
function get_FilterBanks(N::Int64, Q::Vector{Int64}, J::Int64, σ0::Float64; rψ::Float64=sqrt(0.5), α::Float64=5.0)
    σ0s = repeat([σ0], length(Q))
    Js = repeat([J], length(Q))
    get_FilterBanks(N, Q, Js, σ0s; rψ=rψ, α=α)
end

"""
    FilterBank1dBlock

A structure for collections of 1D Morlet Wavelets and a Gaussian lowpass filter,
both stored as arrays.

## Fields: FilterBank1d
 | **Field** | **Description** |
 |:----------|:----------------|
 | :N        | log2(signal length) |
 | :J        | log2(max filter scale) |
 | :Q        | Number of wavelets per octave |
 | :σs       | Bandwidth of filters |
 | :ξs       | Central frequencies of filters |
 | :js       | Maximum acceptable subsampling |
 | :σg       | Bandwidth of gaussian lowpass filter |
 | :ω        | Frequencies |
 | :ϕ        | Fourier transformed low-pass Gaussian filter |
 | :Λ        | Fourier transformed Morlet wavelets |
"""
struct FilterBank1dBlock <: AbstractFilterBank
    N::Int64 # log2 length of filter in time domain
    J::Int # log2 maximum scale of filters
    Q::Int # Number of wavelets per octave

    # these would have been in the MorletWavelet objects
    σs::Vector{Float64}
    ξs::Vector{Float64}
    js::Vector{Int64}
    # this would have been in the GaussianFilter object
    σg::Float64

    ω::Vector{Float64} # frequencies
    ϕ::Vector{Float64} # low-pass filter
    Λ::Array{Float64,2} # Fourier-transformed Morlet wavelets stored in a block

    function FilterBank1dBlock(N::Int64, Q::Int64, J::Int64, σ0::Float64; rψ::Float64=sqrt(0.5), α::Float64=5.0)
        σs, ξs, js = get_filter_params(Q, σ0/2^J, rψ=rψ, α=α)
        ω = FFTW.fftfreq(2^N)

        Λ = Array{Float64,2}(undef, 2^N, length(σs))
        ϕ = GaussianFilter(2^N, σ0/2^J)
        for i=1:length(σs)
            γ = GaborWavelet(2^N, σs[i], ξs[i], js[i])
            Λ[:, i] = MorletWavelet(γ, ϕ).ψ
        end

        return new(N, J, Q, σs, ξs, js, ϕ.σ, ω, ϕ.ϕ, Λ)
    end

    function FilterBank1dBlock(F::FilterBank1d)
        Λ = Array{Float64,2}(undef, 2^F.N, length(F.Λ))
        σs = Vector{Float64}(undef, length(F.Λ))
        ξs = Vector{Float64}(undef, length(F.Λ))
        js = Vector{Int64}(undef, length(F.Λ))
        for i=1:length(F.Λ)
            Λ[:,i].=F.Λ[i].ψ
            σs[i]=F.Λ[i].σ
            ξs[i]=F.Λ[i].ξ
            js[i]=F.Λ[i].j
        end
        return new(F.N, F.J, F.Q, σs, ξs, js, F.ω, F.ϕ.ϕ, Λ)
    end
end
