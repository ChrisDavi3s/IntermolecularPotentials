import numpy as np
import csv

def calculate_field(mus, r):
    '''Calculate the field at each molecule due to all other molecules'''
    num_molecules = len(mus)
    fields = np.zeros(num_molecules)
    for i in range(num_molecules):
        for j in range(num_molecules):
            if i != j:
                fields[i] += 2 * 1/(4*np.pi*epsilon_0_amstrong) *  mus[j] / ((r * abs(i - j))**3)
    #These field values have not been divived by 4*pi*epsilon_0
    return fields

def calculate_induced_dipoles(alpha, fields):
    '''Calculate the induced dipoles at each molecule due to the field at each molecule'''
    return alpha * fields

def induction_energy_per_dipole(fields):
    '''Calculate the induction energy at each molecule due to the induced dipoles at each molecule'''
    induction_energies = []
    for i in range(len(fields)):
        induction_energies.append(-0.5 * alpha_zz * fields[i]**2 * debye_to_eA**2 )
    return induction_energies

def electrostatic_energy_per_dipole(fields, mus, induction_energy):
    '''Calculate the electrostatic energy at each molecule due to the induced dipoles at each molecule'''
    electrostatic_energies = []
    for i in range(len(fields)):
        electrostatic_energies.append(induction_energy[i] - mus[i] * fields[i] * debye_to_eA)
    return electrostatic_energies

def iterate_induced_dipoles(alpha, initial_mus, r, convergence_threshold, filename):
    '''Iterate the induced dipoles until convergence'''
    previous_mus = np.ones(len(initial_mus)) * np.inf
    change = np.inf
    step = 1
    fields = calculate_field(initial_mus, r)
    
    with open(filename, 'w', newline='') as csvfile:

        writer = csv.writer(csvfile)
        writer.writerow(["Step"] + [f"Induced Dipole {i+1} (Debye)" for i in range(5)] + [f'Total Dipole {i+1} (Debye)' for i in range(5)] + ["Induction energy (eV)", "Electrostatic energy (eV)"])

        while change > convergence_threshold:
            induced_mus_new = calculate_induced_dipoles(alpha, fields)
            total_mus_new = induced_mus_new + initial_mus
            fields = calculate_field(total_mus_new, r)
            
            induction_energy = induction_energy_per_dipole(fields)

            electrostatic_energy = electrostatic_energy_per_dipole(fields, mus, induction_energy)
            
            writer.writerow([step] + list(np.round(induced_mus_new,4)) + list(total_mus_new) + [np.sum(induction_energy)] + [np.sum(electrostatic_energy)])
            
            change = np.max(np.abs(total_mus_new - previous_mus))
            previous_mus = total_mus_new.copy()

            step += 1

    return total_mus_new, fields, induction_energy, electrostatic_energy


if __name__ == '__main__':
    # Constants
    convergence_threshold = 1e-4
    mu_initial = 1.8550
    num_molecules = 5
    debye_to_eA = 0.20819434
    epsilon_0_amstrong = 0.0005526349 #is in eV*angstrom
    alpha_zz = 1.4500 * 4 * np.pi *epsilon_0_amstrong #is in Angstrom^3
    alpha_zz_to_meter = 1.4500e-30 #is in meter^3


    # Replace r with the actual distance in angstroms
    distances = [2.5,3.5,4.5] # Replace with the actual value

    for r in distances:
        # Initialize dipoles and positions
        mus = np.ones(num_molecules) * mu_initial

        # Calculate induced dipoles
        filename = f'results_q1_{r}.csv'
        total_mus, fields, induction_energy, electrostatic_energy = iterate_induced_dipoles(alpha_zz, mus, r, convergence_threshold, filename)
        induced_mus = total_mus - mu_initial

        # Print the results
        print("Results for r =", r)
        print("Induced dipoles (Debye):", induced_mus)
        print("Total dipoles (Debye):", total_mus)
        print("Induction energy (eV):", induction_energy)
        print("Electrostatic energy (eV):", electrostatic_energy)
        print("Total Electrostatic energy (eV):", np.sum(electrostatic_energy))
        print()
