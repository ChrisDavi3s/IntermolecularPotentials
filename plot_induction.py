import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV data into pandas DataFrames
df_2_5 = pd.read_csv('results_q3_2.5.csv')
df_3_5 = pd.read_csv('results_q3_3.5.csv')
df_4_5 = pd.read_csv('results_q3_4.5.csv')

# Plot induction energy against step for each DataFrame
plt.plot(df_2_5['Step'], df_2_5['Induction energy (Debye*angstrom)'], label='r=2.5')
plt.plot(df_3_5['Step'], df_3_5['Induction energy (Debye*angstrom)'], label='r=3.5')
plt.plot(df_4_5['Step'], df_4_5['Induction energy (Debye*angstrom)'], label='r=4.5')

# Add labels and legend
plt.xlabel('Step')
plt.ylabel('Induction Energy (Debye*angstrom)')
plt.title('Induction Energy vs Step for different r')
plt.legend()

# Display the plot
plt.show()
