import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV data into pandas DataFrames
df_2_5 = pd.read_csv('results_q(1)_r(2.5).csv')

# For each of the dipoles
for i in range(1, 6):
    # Plot induced dipole change over steps
    plt.plot(df_2_5['Step'], df_2_5[f'Total Dipole {i} (Debye)'], label=f'Dipole at {i}')

plt.xlabel('Step')
plt.ylabel('Dipole (Debye)')
plt.legend()
plt.show()
