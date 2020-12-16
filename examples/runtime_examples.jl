using Scat
using BenchmarkTools
using Profile
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
        Transform = ScatteringTransform1d(2, Q, [J,J])
        times[n] = @elapsed Scat.scatter1d(x, Transform, FilterBanks)
    end

    return times
end
times = time_scatteringTransform_N([12, 1], 18)
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
        Transform = ScatteringTransform1d(2, Q[i,:], [J,J])
        times[i] = @elapsed Scat.scatter1d(x, Transform, FilterBanks)
    end

    return times
end
# vary the number of filters in the first layer
Q = ones(50, 2)
Q[:, 1] .= Q[:, 1].*collect(1:50)
Q = convert.(Int, Q)
N=12
times = time_scatteringTransform_Q(Q, N)
scatter(collect(1:50), times)
# vary the number of filters in the second layer
Q = ones(50, 2)
Q[:, 1] .= Q[:, 1].*collect(1:50)
Q = convert.(Int, Q)
N=12
times = time_scatteringTransform_Q(Q, N)
scatter(collect(1:50), times)

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
        Transform = ScatteringTransform1d(2, Q, [J[i], J[i]])
        times[i] = @elapsed scatter1d(x, Transform, FilterBanks)
    end

    return times
end
N=12
J = collect(4:2:32)
times = time_scatteringTransform_J(J, N)
scatter(J, times)

function time_scatteringTransform_BreadthDepth(N)
    # define the scattering transform
    J=[8,8]
    σ0=0.1
    Q=[12,1]

    times_depth = zeros(N)
    times_breadth = zeros(N)
    for n=1:N
        println(n)
        x=ones(2^n)
        FilterBanks = get_FilterBanks(n, Q, J, σ0)
        FilterBankBlocks = get_FilterBanks(n, Q, J, σ0, typ=FilterBank1dBlock)
        Transform = ScatteringTransform1d(2, Q, J)
        times_depth[n] = @elapsed Scat.scatter1d(x, Transform, FilterBanks)
        times_breadth[n] = @elapsed Scat.scatter1d(x, Transform, FilterBankBlocks)
    end

    return times_depth, times_breadth
end
td, tb = time_scatteringTransform_BreadthDepth(18)
scatter(collect(1:18), td, yscale=:log, label="Depth-first")
scatter!(collect(1:18), tb, yscale=:log, label="Breadth-first")
