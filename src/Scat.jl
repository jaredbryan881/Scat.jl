module Scat

using FFTW

include("./filters/filter.jl")
include("./filters/filterBank.jl")
include("./transform/scatteringTransform.jl")
include("./transform/cwt.jl")
include("./transform/scattered.jl")
include("./transform/scatter1d.jl")
end
