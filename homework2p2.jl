using ElectrochemicalKinetics
using Plots

# Constants
R = 8.314  # Universal gas constant (J/mol.K)
F = 96485  # Faraday's constant (C/mol)
T = 298    # Temperature in Kelvin
alpha = 0.5  # Charge transfer coefficient
j0 = 1e-3   # Exchange current density (A/cm^2)

# Overpotential range
η = range(-0.1, 0.1, length=100)

# Compute the current densities
ja = j0 .* exp.(((1-alpha) * F .* η) ./ (R*T))
jc = -j0 .* exp.((-alpha * F .* η) ./ (R*T))
j = ja .+ jc  # Total current density

# First Plot: Current Density vs. Potential
plot(η, j .* 1000, label="j = ja + jc", color=:black)
plot!(η, ja .* 1000, label="jₐ", color=:red)
plot!(η, jc .* 1000, label="j𝚌", color=:blue)
xlabel!("Potential [V]")
ylabel!("Current Density [mA/cm²]")
title!("Butler-Volmer Equation")
savefig("butler_volmer_plot.png")

# Second Plot: Logarithmic Plot (Tafel Plot)
plot(η, log10.(abs.(j)), label="α = 0.5", linestyle=:solid)
plot!(η, log10.(abs.(j0 .* exp.((-0.3 * F .* η) ./ (R*T)))), label="α = 0.3", linestyle=:dot)
plot!(η, log10.(abs.(j0 .* exp.((-0.1 * F .* η) ./ (R*T)))), label="α = 0.1", linestyle=:dash)
xlabel!("Potential [V]")
ylabel!("Log (Current Density) [log(A/cm²)]")
title!("Tafel Plot")
savefig("tafel_plot.png")