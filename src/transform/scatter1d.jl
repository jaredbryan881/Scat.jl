export scatter1d, scatter1d_layer

# depth-first scattering tree traversal
"""
    scatter1d_layer(S, i, x_hat, Transform, FilterBanks)

Apply one layer of the wavelet scattering transform, calling future layers recursively.

# Arguments
 - `S::Scattered`: Scattered object for holding scattering coefficients.
 - `i::Int64`: Index of the current scattering layer.
 - `x_hat::Vector{Complex{Float64}}`: Fourier transformed input data vector.
 - `Transform::ScatteringTransform1d`: ScatteringTransform1d object containing information about the total transform.
 - `FilterBanks::Vector{FilterBank1d}`: Filter banks for each layer of the scattering transform.
"""
function scatter1d_layer(S::Scattered, i::Int64, nDone::Vector{Int64}, x_hat::Vector{Complex{Float64}}, Transform::ScatteringTransform1d, FilterBanks::Vector{FilterBank1d})
    for (j, ψ) in enumerate(FilterBanks[i].Λ)
        # convolve signal with morlet wavelet, transform back to the time domain, and take the modulus
        U_hat = fft(abs.(ifft(x_hat .* ψ.ψ)))

        # stop scattering when we reach the last layer
        if i < Transform.D
            scatter1d_layer(S, i+1, nDone, U_hat, Transform, FilterBanks)
        end

        # apply lowpass filter
        # result should be real-valued, so cast to real
        S.Coeff[i+1][:, nDone[i]] .= real.(ifft(U_hat .* FilterBanks[i].ϕ.ϕ))
        nDone[i]+=1
    end
end

"""
    scatter1d(x, Transform, FilterBanks)

Calculate the 1d wavelet scattering transform.

# Arguments
 - `x::Vector{Float64}`: Time series.
 - `Transform::ScatteringTransform1d`: Scattering transform object.
 - `FilterBanks::Vector{FilterBank1d}`: Filter banks for each layer of the scattering transform.
"""
function scatter1d(x::Vector{Float64}, Transform::ScatteringTransform1d, FilterBanks::Vector{FilterBank1d})
    @assert length(FilterBanks)==Transform.D "The number of filter banks should equal the scattering transform depth.
                                              Got $(length(FilterBanks)) filter banks and $(Transform.D) layers."

    # Transform the signal to the frequency domain
    x_hat = FFTW.fft(x)
    # create Scattered object to hold the scattering coefficients
    S = Scattered1d(FilterBanks)

    # zeroth order scattering coefficients are a simple averaging of the
    S.Coeff[1] .= real.(ifft(x_hat .* FilterBanks[1].ϕ.ϕ))
    nDone=ones(Int, Transform.D)
    # higher order scattering coefficients are computed recursively
    scatter1d_layer(S, 1, nDone, x_hat, Transform, FilterBanks)

    return S
end

# breadth-first scattering tree traversal
# TODO: allow any number of layers with this approach
"""
    scatter1d(x, Transform, FilterBanks)

Calculate the 1d wavelet scattering transform.

# Arguments
 - `x::Vector{Float64}`: Time series.
 - `Transform::ScatteringTransform1d`: Scattering transform object.
 - `FilterBanks::Vector{FilterBank1d}`: Filter banks for each layer of the scattering transform.
"""
function scatter1d(x::Vector{Float64}, Transform::ScatteringTransform1d, FilterBanks::Vector{FilterBank1dBlock})
    @assert length(FilterBanks)==Transform.D "The number of filter banks should equal the scattering transform depth.
                                              Got $(length(FilterBanks)) filter banks and $(Transform.D) layers."
    @assert Transform.D <=2 "The number of layers should be 2 or fewer. Got $(Transform.D)"

    # Transform the signal to the frequency domain
    x_hat = FFTW.fft(x)
    # create Scattered object to hold the scattering coefficients
    S = Scattered1d(FilterBanks)

    # zeroth order scattering coefficients are a simple averaging of the
    S.Coeff[1] .= real.(ifft(x_hat .* FilterBanks[1].ϕ))

    # first order scattering coefficients
    U1 = Array{Float64,2}(undef, size(FilterBanks[1].Λ))
    U1 .= abs.(cwt(x, FilterBanks[1]))
    S.Coeff[2] .= real.(ifft(fft(U1, 1) .* FilterBanks[1].ϕ, 1))

    # second order scattering coefficients
    U2i = Array{Float64,2}(undef, size(FilterBanks[2].Λ))
    nDone=0
    for i=1:size(U1,2)
        U2i .= abs.(cwt(U1[:,i], FilterBanks[2]))
        S.Coeff[3][:, nDone+1:nDone+size(U2i,2)] .= real.(ifft(fft(U2i,1) .* FilterBanks[2].ϕ, 1))
        nDone+=size(U2i,2)
    end

    return S
end
