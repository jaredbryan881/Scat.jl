abstract type ScatteringTransform end

struct ScatteringTransform1d <: ScatteringTransform
    D::Int # depth/number of layers
    J::Int # number of scales/number of wavelets in the filterbank
end
