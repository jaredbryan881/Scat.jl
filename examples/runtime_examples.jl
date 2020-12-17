using Scat
using BenchmarkTools
using Plots

function time_subsamplingSize(Qs::Vector{Int64}, N::Int64; nTrials=1)
    # define the scattering transform
    # first layer parameters and filter bank
    J=8
    σ0=0.1

    times_subsamp = zeros(N)
    times_nosubsamp = zeros(N)
    for i=1:nTrials
        for n=1:N
            println(n)
            x=ones(2^n)
            Transform = ScatteringTransform1d(2, Qs, [J,J])
            FilterBanks = get_FilterBanks(n, Qs, J, σ0)
            times_subsamp[n] += @elapsed Scat.scatter1d(x, Transform, FilterBanks, subsample=true)
            times_nosubsamp[n] += @elapsed Scat.scatter1d(x, Transform, FilterBanks, subsample=false)
        end
    end
    times_subsamp./=nTrials
    times_nosubsamp./=nTrials

    return times_subsamp, times_nosubsamp
end
times_ss, times_nss = time_subsamplingSize([12, 1], 18, nTrials=5)
scatter(collect(1:18).+log2(64), times_nss, yscale=:log10, label="No subsampling", legend=:topleft)
scatter!(collect(1:18).+log2(64), times_ss, yscale=:log10, label="With subsampling")
xlabel!("log2(Input size [bits])")
ylabel!("Time [s]")
savefig("time_subsampling_size.png")
mean(times_nss./times_ss)

function time_scatteringTransform_Q(Q::Array{Int64,2}, N::Int64)
    # define the scattering transform
    # first layer parameters and filter bank
    J=8
    σ0=0.1
    x=ones(2^N)

    times = zeros(size(Q,1))
    nFilters = zeros(size(Q,1))
    for i=1:size(Q,1)
        FilterBanks = get_FilterBanks(N, Q[i,:], J, σ0)
        Transform = ScatteringTransform1d(2, Q[i,:], [J,J])
        println(length(FilterBanks[1].Λ) * length(FilterBanks[2].Λ))
        nFilters[i] = length(FilterBanks[1].Λ)*length(FilterBanks[2].Λ)
        times[i] = @elapsed Scat.scatter1d(x, Transform, FilterBanks)
    end

    return nFilters, times
end

# vary the number of filters in the first layer
nQ=100
Q = ones(nQ, 2)
#Q[:, 1] .= Q[:, 1].*collect(1:100)
Q[:,1].= rand(1:12,nQ)
Q[:,2].=rand(1:4,nQ)
Q = convert.(Int, Q)
N=12
nFilters, times = time_scatteringTransform_Q(Q, N)
scatter(Q[:,1].*Q[:,2], times, legend=:none)
xlabel!("Q1 ⋅ Q2")
ylabel!("Time [s]")

scatter(nFilters, times, legend=:none)
xlabel!("Total Number of Filters")
ylabel!("Time [s]")
savefig("nFilters_times.png")

function time_subsamplingQ(Q::Array{Int64,2}, N::Int64)
    # define the scattering transform
    # first layer parameters and filter bank
    J=8
    σ0=0.1
    x=ones(2^N)

    times_ss = zeros(size(Q,1))
    times_nss = zeros(size(Q,1))
    nFilters = zeros(size(Q,1))
    for i=1:size(Q,1)
        FilterBanks = get_FilterBanks(N, Q[i,:], J, σ0)
        Transform = ScatteringTransform1d(2, Q[i,:], [J,J])
        println(length(FilterBanks[1].Λ) * length(FilterBanks[2].Λ))
        nFilters[i] = length(FilterBanks[1].Λ)*length(FilterBanks[2].Λ)
        times_ss[i] = @elapsed Scat.scatter1d(x, Transform, FilterBanks)
        times_nss[i] = @elapsed Scat.scatter1d(x, Transform, FilterBanks, subsample=false)
    end

    return nFilters, times_ss, times_nss
end
# vary the number of filters in the first layer
nQ=100
Q = ones(nQ, 2)
#Q[:, 1] .= Q[:, 1].*collect(1:100)
Q[:,1].= rand(1:12,nQ)
Q[:,2].=rand(1:4,nQ)
Q = convert.(Int, Q)
N=12
nFilters, times_ss, times_nss = time_subsamplingQ(Q, N)
scatter(nFilters, times_nss, label="No subsampling")
scatter!(nFilters, times_ss, legend=:topleft, label="With subsampling")
xlabel!("Q1 ⋅ Q2")
ylabel!("Time [s]")
savefig("time_subsampling_nFilters.png")

function time_scatteringTransform_J(J::Vector{Int64}, N::Int64)
    # define the scattering transform
    # first layer parameters and filter bank
    σ0=0.1
    Q=[12, 2]

    times = zeros(length(J))
    for i=1:length(J)
        println(J[i])
        x=ones(2^N)
        FilterBanks = get_FilterBanks(N, Q, J[i], σ0)
        Transform = ScatteringTransform1d(2, Q, [J[i], J[i]])
        times[i] = @elapsed Scat.scatter1d(x, Transform, FilterBanks)
    end

    return times
end
N=12
J = collect(4:20)
times = time_scatteringTransform_J(J, N)
scatter(J, sqrt.(times), legend=:none)
xlabel!("log2(maximum scale)")
ylabel!("√Time")
savefig("maxScale_time.png")

function time_subsamplingJ(J::Vector{Int64}, N::Int64; nTrials=1)
    # define the scattering transform
    # first layer parameters and filter bank
    σ0=0.1
    Q=[12, 2]

    times = zeros(length(J))
    for n=1:nTrials
        for i=1:length(J)
            println(J[i])
            x=ones(2^N)
            FilterBanks = get_FilterBanks(N, Q, J[i], σ0)
            Transform = ScatteringTransform1d(2, Q, [J[i], J[i]])
            times[i] += @elapsed Scat.scatter1d(x, Transform, FilterBanks)
        end
    end

    return times./nTrials
end
N=12
J = collect(4:16)
times = time_subsamplingJ(J, N, nTrials=15)
scatter(J, sqrt.(times), legend=:none)
xlabel!("log2(maximum scale)")
ylabel!("√Time")
savefig("maxScale_time.png")

function time_scatteringTransform_BreadthDepth(N)
    # define the scattering transform
    J=[8,8]
    σ0=0.1
    Q=[12,6]

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
td, tb = time_scatteringTransform_BreadthDepth(14)
scatter(collect(1:14).+log2(64), td, yscale=:log10, label="Depth-first", legend=:topleft)
scatter!(collect(1:14).+log2(64), tb, yscale=:log10, label="Breadth-first")
xlabel!("log2(Input size) [bits]")
ylabel!("Time [s]")
savefig("depth_vs_breadth_time.png")
