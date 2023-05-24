import numpy as np
import csv

def calculate_field(mus, r):
    '''Calculate the field at each molecule due to all other molecules'''
    num_molecules = len(mus)
    fields = np.zeros(num_molecules)
    for i in range(num_molecules):
        for j in range(num_molecules):
            if i != j:
                #Field due to other dipoles
                fields[i] += 2 * 1/(4*np.pi*epsilon_0_amstrong) *  mus[j] / ((r * abs(i - j))**3)
        if q != 0:
            #Field due to the charge
            # The charge is at a distance of -r from the first molecule, so it's at a distance of -r + i*r from the i-th molecule
            fields[i] -= (q / ((r + i*r) ** 2)) / (eA_per_Debye * 4 * np.pi * epsilon_0_amstrong)
    return fields

def calculate_induced_dipoles(alpha, fields):
    '''Calculate the induced dipoles at each molecule due to the field at each molecule'''
    return alpha * fields

def induction_energy_per_dipole(fields):
    '''Calculate the induction energy at each molecule due to the induced dipoles at each molecule'''
    induction_energies = []
    for i in range(len(fields)):
        induction_energies.append(-0.5 * alpha_zz * fields[i]**2 * eA_per_Debye**2 )
    return induction_energies

def electrostatic_energy_per_dipole(fields, mus):
    '''Calculate the electrostatic energy at each molecule due to the induced dipoles at each molecule'''
    electrostatic_energies = []
    for i in range(len(fields)):
        electrostatic_energies.append(- mus[i]* fields[i] * eA_per_Debye**2)
    return electrostatic_energies

def iterate_induced_dipoles(alpha, initial_mus, r, convergence_threshold, filename):
    '''Iterate the induced dipoles until convergence'''
    previous_mus = np.ones(len(initial_mus)) * np.inf
    change = np.inf
    step = 1
    fields = calculate_field(initial_mus, r)
    
    with open(filename, 'w', newline='') as csvfile:

        # Write the header
        writer = csv.writer(csvfile)
        writer.writerow(["Step"] + [f"Induced Dipole {i+1} (Debye)" for i in range(5)] + [f'Total Dipole {i+1} (Debye)' for i in range(5)] + ["Induction energy (eV)", "Electrostatic energy (eV)"])
        
        # Iterate until convergence
        while change > convergence_threshold:
            
            induced_mus_new = calculate_induced_dipoles(alpha, fields)
            total_mus_new = induced_mus_new + initial_mus
            fields = calculate_field(total_mus_new, r)
            induction_energy = induction_energy_per_dipole(fields)
            electrostatic_energy = electrostatic_energy_per_dipole(fields, total_mus_new)
            
            # Write the results to the csv file
            writer.writerow([step] + list(np.round(induced_mus_new,4)) + list(np.round(total_mus_new,4)) + [np.round(np.sum(induction_energy),4)] + [np.round(np.sum(electrostatic_energy),4)])
            
            change = np.max(np.abs(total_mus_new - previous_mus))
            previous_mus = total_mus_new.copy()
            step += 1

    return total_mus_new, fields, induction_energy, electrostatic_energy

if __name__ == '__main__':
    # Constants
    convergence_threshold = 1e-4
    mu_initial = 1.8550
    num_molecules = 5
    eA_per_Debye = 0.20819434
    epsilon_0_amstrong = 0.005526349 
    alpha_zz = 1.4500 * 4 * np.pi *epsilon_0_amstrong

    # Parameters
    distances = [2.5,3.5,4.5] #distances in Angstrom
    qs = [0,1] #charges

    for q in qs:
            for r in distances:
                # Initialize dipoles and positions for each molecule
                mus = np.ones(num_molecules) * mu_initial

                # Calculate induced dipoles
                filename = f'results_q({q})_r({r}).csv' #filename for the results
                total_mus, fields, induction_energy, electrostatic_energy = iterate_induced_dipoles(alpha_zz, mus, r, convergence_threshold, filename)
                induced_mus = total_mus - mu_initial

                # Print the results - same information as in the csv file
                print("Results for r =", r)
                print("q =", q)
                print("Induced dipoles (Debye):", np.round(induced_mus,4))
                print("Total dipoles (Debye):", np.round(total_mus,4))
                print("Induction energy (eV):", np.round(induction_energy,4))
                print("Electrostatic energy (eV):", np.round(electrostatic_energy,4))
                print("Total Induction energy (eV):", np.round(np.sum(induction_energy),4))
                print("Total Electrostatic energy (eV):", np.round(np.sum(electrostatic_energy),4))
                print()