export Scattered1d

abstract type Scattered end

"""
    Scattered1d

A structure for scattering coefficients.

## Fields: Scattered1d
 | **Field** | **Description** |
 |:----------|:----------------|
 | :T        | signal length |
 | :nFilters | number of filters per layer |
 | :Coeff    | scattering coefficients |
"""
struct Scattered1d <: Scattered
    T::Int64 # number of samples in time
    nFilters::Vector{Int64} # number of filters in each layer
    Coeff::Vector{Array{Float64,2}} # scattering coefficients for each layer

    function Scattered1d(F::Vector{FilterBank1d})
        nFilters = Vector{Int64}(undef, length(F)+1)
        Coeff = Vector{Array{Float64,2}}(undef, length(F)+1)

        # 0th order scattering coefficients
        nFilters[1] = 1
        T = 2^F[1].N
        Coeff[1] = Array{Float64,2}(undef, T, 1)

        # higher order scattering coefficients
        for i=1:length(F)
            nFilters[i+1] = length(F[i].Λ)
            if i>=2
                # count the number of js in the second layer greater than js in the first layer
                nVect=sum([length([ψ2.j for ψ2 in F[i].Λ if ψ2.j>ψ1.j]) for ψ1 in F[i-1].Λ])
            else
                nVect=length(F[i].Λ)
            end
            Coeff[i+1] = Array{Float64,2}(undef, T, nVect)
        end

        return new(T, nFilters, Coeff)
    end

    # TODO: Reduce the number of convolutions by counting the js[i]<js[i-1]
    function Scattered1d(F::Vector{FilterBank1dBlock})
        nFilters = Vector{Int64}(undef, length(F)+1)
        Coeff = Vector{Array{Float64,2}}(undef, length(F)+1)

        # 0th order scattering coefficients
        nFilters[1] = 1
        T = 2^F[1].N
        Coeff[1] = Array{Float64,2}(undef, T, 1)

        # higher order scattering coefficients
        nRows=1
        for i=1:length(F)
            nFilters[i+1] = size(F[i].Λ, 2)
            nRows=nRows*size(F[i].Λ, 2)
            Coeff[i+1] = Array{Float64,2}(undef, T, nRows)
        end

        return new(T, nFilters, Coeff)
    end
end
