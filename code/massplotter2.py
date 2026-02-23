import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('data/solution_n-3_alpha-0.0_C1.0.dat')

time = data[:, 0]
all_masses = data[:, 1:] #all masses including walls
N_masses = all_masses.shape[1] - 2 #just the main masses (from simul.init)
interior_positions = all_masses[:, 1:-1]

half_mass = N_masses // 2  #mass at halfway point

T_f = time[-1] # Final time 

plt.figure(figsize=(14, 7))

#main plot
plt.plot(time, interior_positions[:, half_mass], 
         linewidth=2, label=f'Mass {half_mass} (Middle)', color='blue', alpha=0.7)

#labels for time increments (red lines)
plt.text(T_f/4,-0.185, 'Tf/4',color='r')
plt.text(T_f/2,-0.185, 'Tf/2',color='r')
plt.text((3*T_f)/4,-0.185, '3Tf/4',color='r')
plt.text(T_f,-0.185, 'Tf',color='r')
#plotting red lines at time increments
plt.axvline(x=T_f/4,color='r')
plt.axvline(x=T_f/2,color='r')
plt.axvline(x=(3*T_f)/4,color='r')
plt.axvline(x=T_f,color='r')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Position', fontsize=12)
plt.title(f'Position vs Time - Middle Mass (Mass {half_mass} of {N_masses}), C = 0.9515, alpha = -n/10', fontsize=14)
plt.legend(fontsize=11, loc='upper right')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('halfmass.png', dpi=150)
plt.show()