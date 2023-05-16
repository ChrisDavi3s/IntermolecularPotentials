import numpy as np

# Physical and conversion constants
DEBYE_TO_COULOMB_METER = 3.33564e-30
ANGSTROM_TO_METER = 1e-10
EPSILON_0 = 8.8541878128e-12  # Permittivity of free space in F/m
K_E = 8.9875517923e9  # Coulomb's constant, N m^2 C^-2
ELEMENTARY_CHARGE = 1.602176634e-19

def polarizability_conversion_factor(angstrom_cubed):
    return (angstrom_cubed * ANGSTROM_TO_METER**3) * 4 * np.pi * EPSILON_0

def calculate_field(mus, positions, r, charge):
    num_molecules = len(mus)
    fields = np.zeros(num_molecules)

    for i in range(num_molecules):
        for j in range(num_molecules):
            if i != j:
                fields[i] += 2 * DEBYE_TO_COULOMB_METER * mus[j] / ((r * ANGSTROM_TO_METER * abs(i - j))**3)

        # Add field contribution from the charge
        dist_to_charge = np.linalg.norm(positions[i] - charge['position'])
        fields[i] += K_E * charge['magnitude'] * ELEMENTARY_CHARGE / (dist_to_charge**3)

    return fields

def calculate_induced_dipoles(alpha, fields):
    return alpha * fields

def calculate_induction_energy(induced_mus, fields):
    return -0.5 * np.sum(DEBYE_TO_COULOMB_METER * induced_mus * fields)

def iterate_induced_dipoles(alpha, initial_mus, positions, r, convergence_threshold, charge):
    total_mus = initial_mus.copy()
    change = np.inf
    step = 1

    while change > convergence_threshold:
        fields = calculate_field(total_mus, positions, r, charge)
        induced_mus = calculate_induced_dipoles(alpha, fields)
        total_mus_new = initial_mus + induced_mus
        change = np.max(np.abs(DEBYE_TO_COULOMB_METER * (total_mus_new - total_mus)))
        total_mus = total_mus_new

        induction_energy = calculate_induction_energy(induced_mus, fields)

        print(f"Induction step: {step}")
        print("Induced dipoles (C*m):", DEBYE_TO_COULOMB_METER * induced_mus)
        print("Total dipoles (C*m):", DEBYE_TO_COULOMB_METER * total_mus)
        print("Induction energy (Joules):", induction_energy)
        print()

        step += 1

    return total_mus

# Constants
convergence_threshold = 1e-4 * DEBYE_TO_COULOMB_METER
mu_initial = DEBYE_TO_COULOMB_METER * 1.8550
alpha_zz = polarizability_conversion_factor(1.4500)
num_molecules = 5
r = ANGSTROM_TO_METER * 2.5

# Initialize dipoles and positions
mus = np.ones(num_molecules) * mu_initial
positions = np.array([[0, 0, i * r] for i in range(num_molecules)])

# Define charge properties
charge = {
    'magnitude': 0,  # In units of elementary charge
    'position': np.array([0, 0, -r])  # In meters
}

# Calculate induced dipoles
total_mus = iterate_induced_dipoles(alpha_zz, mus, positions, r, convergence_threshold, charge)
induced_mus = total_mus - mu_initial

# Calculate the induction energy
fields = calculate_field(total_mus, positions, r, charge)
induction_energy = calculate_induction_energy(induced_mus, fields)

# Results
print("Final induced dipoles (C*m):", induced_mus * DEBYE_TO_COULOMB_METER)
print("Final total dipoles (C*m):", total_mus * DEBYE_TO_COULOMB_METER)
print("Final induction energy (Joules):", induction_energy)
