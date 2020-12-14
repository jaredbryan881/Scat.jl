using FFTW
using Plots

# Plot filters for the first layer
N = 12 # log2 of the length of the signal
Q = 12 # number of wavelets per octave
J = 8 # log2 of the maximum scale
σ0 = 0.1

# get filter bank for the first layer
FilterBank1 = FilterBank1d(N, Q, J, σ0)

ω = FFTW.fftshift(FilterBank1.ω)
plot(ω, FFTW.fftshift(FilterBank1.Λ[1].ψ), xlims=(0.0, 0.5), color=:black, linewidth=2, legend=false)
for ψ in FilterBank1.Λ[2:end]
    plot!(ω, FFTW.fftshift(ψ.ψ), xlims=(0.0, 0.5), color=:black, linewidth=2, legend=false)
end
plot!(ω, FFTW.fftshift(FilterBank1.ϕ.ϕ), xlims=(0.0, 0.5), linewidth=2, color=:red, legend=false)
xlabel!("Frequency")


# Plot filters for the second layer
N = 12 # log2 of the length of the signal
Q = 1 # number of wavelets per octave
J = 8 # log2 of the maximum scale
σ0 = 0.1

# get filter bank for the second layer
FilterBank2 = FilterBank1d(N, Q, J, σ0)

ω = FFTW.fftshift(FilterBank2.ω)
plot(ω, FFTW.fftshift(FilterBank2.Λ[1].ψ), xlims=(0.0, 0.5), color=:black, linewidth=2, legend=false)
for ψ in FilterBank2.Λ[2:end]
    plot!(ω, FFTW.fftshift(ψ.ψ), xlims=(0.0, 0.5), color=:black, linewidth=2, legend=false)
end
plot!(ω, FFTW.fftshift(FilterBank2.ϕ.ϕ), xlims=(0.0, 0.5), linewidth=2, color=:red, legend=false)
xlabel!("Frequency")
