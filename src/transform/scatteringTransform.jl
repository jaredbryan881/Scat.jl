export ScatteringTransform1d, GraphScatteringTransform

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

"""
    GraphScatteringTransform

A structure for the graph scattering transform for a collection of 1d signals.

## Fields: GraphScatteringTransform
 | **Field** | **Description** |
 |:----------|:----------------|
 | :D        | number of layers |
 | :Q        | Number of wavelets per octave |
 | :J        | log2(max filter scale) |
 | :W        | graph adjacency matrix |
"""
struct GraphScatteringTransform <: ScatteringTransform
    D::Int # depth/number of layers
    Q::Vector{Int} # number of filters per octave for each layer
    J::Int # log2 maximum scale of filters for each layer
    W::Array{Float64,2} # adjacency matrix

    # make sure adjacency matrix is square
    GraphScatteringTransform(D,J,W) = size(W)[1]!=size(W)[2] ? error("Nonsquare adjacency matrix") : new(D,J,W)
end