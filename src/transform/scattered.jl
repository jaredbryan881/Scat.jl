abstract type Scattered end

struct Scattered1d <: Scattered
    T::Int64 # number of samples in time
    nFilters::Vector{Int64} # number of filters in each layer
    Coeff::Vector{Array{Float64,2}} # scattering coefficients for each layer

    function Scattered1d(T::Int64, F::Vector{FilterBank1d})
        nFilters = Vector{Int64}(undef, length(F)+1)
        Coeff = Vector{Array{Float64,2}}(undef, length(F)+1)

        # 0th order scattering coefficients
        nFilters[1] = 1
        Coeff[1] = Array{Float64,2}(undef, T, 1)
        # higher order scattering coefficients
        for i=1:length(F)
            Coeff[i+1] = Array{Float64,2}(undef, T, length(F[i].Λ))
        end

        return new(T, nFilters, Coeff)
    end
end
