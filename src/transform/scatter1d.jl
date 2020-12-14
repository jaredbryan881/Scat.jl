function scatter1d_layer(S::Scattered, i::Int64, x_hat::Vector{Complex{Float64}}, Transform::ScatteringTransform1d, FilterBanks::Vector{FilterBank1d})
    @assert size(S.Coeff[i+1], 2) == length(FilterBanks[i].Λ) "Array of scattering coefficients doesn't have the same dimensions as the number of filters.
                                                               Got $(size(S.Coeff[i+1], 2)) and $(length(FilterBanks[i].Λ)), respectively."
    for (j, ψ) in enumerate(FilterBanks[i].Λ)
        # convolve signal with morlet wavelet, transform back to the time domain, and take the modulus
        U_hat = fft(abs.(ifft(x_hat .* ψ.ψ)))

        # stop scattering when we reach the last layer
        if i < Transform.D
            scatter1d_layer(S, i+1, U_hat, Transform, FilterBanks)
        end

        # apply lowpass filter
        # result should be real-valued, so cast to real
        S.Coeff[i+1][:, j] .= real.(ifft(U_hat .* FilterBanks[i].ϕ.ϕ))
    end
end


function scatter1d(x::Vector{Float64}, Transform::ScatteringTransform1d, FilterBanks::Vector{FilterBank1d})
    @assert length(FilterBanks)==Transform.D "The number of filter banks should equal the scattering transform depth.
                                              Got $(length(FilterBanks)) filter banks and $(Transform.D) layers."

    # Transform the signal to the frequency domain
    x_hat = FFTW.fft(x)
    # create Scattered object to hold the scattering coefficients
    S = Scattered1d(length(x), FilterBanks)

    # zeroth order scattering coefficients are a simple averaging of the
    S.Coeff[1] .= real.(ifft(x_hat .* FilterBanks[1].ϕ.ϕ))
    # higher order scattering coefficients are computed recursively
    scatter1d_layer(S, 1, x_hat, Transform, FilterBanks)

    return S
end
