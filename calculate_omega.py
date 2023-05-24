import numpy as np

#constants
distances = [2.5,3.5,4.5] 
debye_to_eA = 0.20819434
epsilon_0_amstrong = 0.005526349 
alpha = 2 * np.pi *epsilon_0_amstrong

#inputs
E = [6.9482,1.258,0.4921]
#convert to units used in main code
E = [E[i] / debye_to_eA**2 for i in range(len(E))]

for r in distances:
    omega = 0
    for i in range(5):
        for j in range(5):
            if i != j:
                # Sum over all molecules where i != j
                omega += 1 / ((r * abs(i - j))**3)
    index = distances.index(r)
    static_dipole = np.sqrt(alpha*E[index] / (omega))

    #print results
    print(f"r = {r} Å, Omega = {omega}, Static dipole = {static_dipole} Debye")