import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV data into pandas DataFrames
df_2_5 = pd.read_csv('results_q(0)_r(2.5).csv')
df_3_5 = pd.read_csv('results_q(0)_r(3.5).csv')
df_4_5 = pd.read_csv('results_q(0)_r(4.5).csv')
df_2_5_1 = pd.read_csv('results_q(1)_r(2.5).csv')
df_3_5_1 = pd.read_csv('results_q(1)_r(3.5).csv')
df_4_5_1 = pd.read_csv('results_q(1)_r(4.5).csv')

# Plot induction energy against step for each DataFrame
plt.plot(df_2_5['Step'], df_2_5['Induction energy (eV)'], label='q=0 r=2.5')
plt.plot(df_3_5['Step'], df_3_5['Induction energy (eV)'], label='q=0 r=3.5')
plt.plot(df_4_5['Step'], df_4_5['Induction energy (eV)'], label='q=0 r=4.5')

# Plot induction energy against step for each DataFrame
plt.plot(df_2_5_1['Step'], df_2_5_1['Induction energy (eV)'], label='q=1 r=2.5')
plt.plot(df_3_5_1['Step'], df_3_5_1['Induction energy (eV)'], label='q=1 r=3.5')
plt.plot(df_4_5_1['Step'], df_4_5_1['Induction energy (eV)'], label='q=1 r=4.5')

# Add labels and legend
plt.xlabel('Step')
plt.ylabel('Induction Energy (eV)')
plt.legend()

# Display the plot
plt.show()
