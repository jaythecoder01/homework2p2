using Printf

# Constants
R = 8.314  # Universal gas constant (J/mol.K)
F = 96485  # Faraday constant (C/mol)

# Function to calculate the Nernst potential
function nernst_potential(E0, z, T, reactants, products)
    # Convert temperature to Kelvin if given in Celsius
    if T < 200
        T += 273.15
    end
    # Calculate reaction quotient Q
    Q = prod(products) / prod(reactants)
    # Calculate the Nernst potential
    E = E0 - (R * T / (z * F)) * log(Q)
    return E
end

# Example calculations
function calculate_ocv()
    println("Solid Oxide Fuel Cell (SOFC):")
    E_SOFC = nernst_potential(1.229, 2, 800, [0.21], [1.0])  # Example for oxygen reduction
    @printf("  OCV at 800K: %.3f V\n", E_SOFC)

    println("\nLi-Ion Battery:")
    E_Li = nernst_potential(3.8, 1, 298, [0.1], [1.0])  # Example for Li-ion chemistry
    @printf("  OCV at 25°C: %.3f V\n", E_Li)
end

calculate_ocv()
