export cwt

"""
    cwt(signal, Λ)

Calculate the continuous wavelet transform of a signal given filter bank Λ.

# Arguments
 - `signal::Vector{Float64}`: Time series.
 - `Λ::Array{Float64, 2}`: Fourier transformed filters arranged by [frequency, scale].
"""
function cwt(signal::Vector{Float64}, Λ::Array{Float64,2})
    signal_ft = fft(signal, 1)
    W = ifft(signal_ft.*Λ, 1)
    return W
end

"""
    cwt(signal, Λ, ω0, δt, δj)

Convenience method for cwt that constructs the filter bank Λ given a signal.

# Arguments
 - `signal::Vector{Float64}`: Time series.
 - `F::FilterBank1dBlock`: Fourier transformed filters and metadata
"""
function cwt(signal::Vector{Float64}, F::FilterBank1dBlock)
    W = cwt(signal, F.Λ)
    return W
end
