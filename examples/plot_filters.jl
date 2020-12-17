using FFTW
using Plots
using Scat

# plot filters in the frequency domain
# Plot filters for the first layer
N = 12 # log2 of the length of the signal
Q = [12, 1] # number of wavelets per octave
J = 8 # log2 of the maximum scale
σ0 = 0.1

FilterBank1, FilterBank2 = get_FilterBanks(N, Q, J, σ0)

ω = FFTW.fftshift(FilterBank1.ω)
plot(ω, FFTW.fftshift(FilterBank1.Λ[1].ψ), xlims=(0.0, 0.5), color=:steelblue, linewidth=2, legend=false)
for ψ in FilterBank1.Λ[2:end]
    plot!(ω, FFTW.fftshift(ψ.ψ), xlims=(0.0, 0.5), color=:steelblue, linewidth=2, legend=false)
end
plot!(ω, FFTW.fftshift(FilterBank1.ϕ.ϕ), xlims=(0.0, 0.5), linewidth=2, color=:crimson, legend=false)
xlabel!("Frequency")
savefig("FilterBank1.png")

# Plot filters for the second layer
ω = FFTW.fftshift(FilterBank2.ω)
plot(ω, FFTW.fftshift(FilterBank2.Λ[1].ψ), xlims=(0.0, 0.5), color=:steelblue, linewidth=2, legend=false)
for ψ in FilterBank2.Λ[2:end]
    plot!(ω, FFTW.fftshift(ψ.ψ), xlims=(0.0, 0.5), color=:steelblue, linewidth=2, legend=false)
end
plot!(ω, FFTW.fftshift(FilterBank2.ϕ.ϕ), xlims=(0.0, 0.5), linewidth=2, color=:crimson, legend=false)
xlabel!("Frequency")
savefig("FilterBank2.png")

# plot filter in the time domain
plot(real.(FFTW.fftshift(FFTW.ifft(FilterBank1.Λ[74].ψ))), linewidth=2, color=:steelblue, label="Real", axis=nothing)
plot!(imag.(FFTW.fftshift(FFTW.ifft(FilterBank1.Λ[74].ψ))), linewidth=2, linestyle=:dash, color=:black, label="Imaginary", axis=nothing)
xlabel!("Time")
savefig("MorletWavelet.png")
