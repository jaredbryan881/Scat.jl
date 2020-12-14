module Scat

using FFTW
using Wavelets

include("./filters/filter.jl")
include("./filters/filterBank.jl")
include("./transform/scatteringTransform.jl")
end
