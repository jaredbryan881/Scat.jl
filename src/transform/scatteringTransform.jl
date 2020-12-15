export ScatteringTransform1d

abstract type ScatteringTransform end

"""
    ScatteringTransform1d

A structure for the wavelet scattering transform for 1d signals.

## Fields: ScatteringTransform1d
 | **Field** | **Description** |
 |:----------|:----------------|
 | :D        | number of layers |
 | :Q        | Number of wavelets per octave |
 | :J        | log2(max filter scale) |
"""
struct ScatteringTransform1d <: ScatteringTransform
    D::Int # depth/number of layers
    Q::Vector{Int} # number of filters per octave for each layer
    J::Vector{Int} # log2 maximum scale of filters for each layer
end
depth(S::ScatteringTransform1d) = S.D
