using Scat
using BenchmarkTools
using Plots

function time_scatteringTransform_N(Qs::Vector{Int64}, N::Int64)
    # define the scattering transform
    # first layer parameters and filter bank
    J=8
    σ0=0.1

    times = zeros(N)
    for n=1:N
        println(n)
        x=ones(2^n)
        FilterBanks = get_FilterBanks(n, Qs, J, σ0)
        Transform = ScatteringTransform1d(2, 1)
        times[n] = @elapsed scatter1d(x, Transform, FilterBanks)
    end

    return times
end
times = time_scatteringTransform_N([12, 1], 16)
scatter(collect(1:18), times, yscale=:log)

function time_scatteringTransform_Q(Q::Array{Int64,2}, N::Int64)
    # define the scattering transform
    # first layer parameters and filter bank
    J=8
    σ0=0.1

    times = zeros(size(Q,1))
    for i=1:size(Q,1)
        println(Q[i,:])
        x=ones(2^N)
        FilterBanks = get_FilterBanks(N, Q[i,:], J, σ0)
        Transform = ScatteringTransform1d(2, 1)
        times[i] = @elapsed scatter1d(x, Transform, FilterBanks)
    end

    return times
end

# vary the number of filters in the first layer
Q = ones(50, 2)
Q[:, 1] .= Q[:, 1].*collect(1:50)
Q = convert.(Int, Q)
N=12
times = time_scatteringTransform_Q(Q, N)
scatter(collect(2:2:50), times)

# vary the number of filters in the second layer
Q = ones(50, 2)
Q[:, 1] .= Q[:, 1].*collect(1:50)
Q = convert.(Int, Q)
N=12
times = time_scatteringTransform_Q(Q, N)
scatter(collect(2:2:50), times)

function time_scatteringTransform_J(J::Vector{Int64}, N::Int64)
    # define the scattering transform
    # first layer parameters and filter bank
    σ0=0.1
    Q=[12, 1]

    times = zeros(length(J))
    for i=1:length(J)
        println(J[i])
        x=ones(2^N)
        FilterBanks = get_FilterBanks(N, Q, J[i], σ0)
        Transform = ScatteringTransform1d(2, 1)
        times[i] = @elapsed scatter1d(x, Transform, FilterBanks)
    end

    return times
end
N=12
J = collect(4:2:32)
times = time_scatteringTransform_J(J, N)
scatter(J, times)
