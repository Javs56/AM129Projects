import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('data/solution_n-3_alpha-0.0_C1.0.dat')

time = data[:, 0]
positions = data[:, 1:] 
nMasses = positions.shape[1] - 2 
 
interior_positions = positions[:, 1:-1]

plt.figure(figsize=(12, 6))

# Plotting all masses (not walls)
for i in range(nMasses):
    plt.plot(time, interior_positions[:, i], 
             label=f'Mass {i}', alpha=0.7, linewidth=2, marker='o', markersize=3)

plt.xlabel('Time', fontsize=12)
plt.ylabel('Position', fontsize=12)
plt.title(f'Position vs Time - {nMasses}', fontsize=14)
plt.legend(loc='best', fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('multiple_masses.png', dpi=150)
plt.show()