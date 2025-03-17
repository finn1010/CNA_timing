import numpy as np
import matplotlib.pyplot as plt

mu_values = np.linspace(0, 0.1, 1000)  
T_values = 1 / (2 * mu_values)

mu_point = 0.025
T_point = 1 / (2 * mu_point)

plt.figure(figsize=(7, 5))

plt.fill_between(mu_values, 0, T_values, color='blue', alpha=0.5, label='Identifiable')
plt.fill_between(mu_values, T_values, 100, color='red', alpha=0.5, label='Non-identifiable')

plt.plot(mu_values, T_values, color='black', label=r'$T = \frac{1}{\mu+\gamma}$', linewidth=1)

plt.plot([mu_point, mu_point], [0, T_point], linestyle='--', color='black', linewidth=1, label=r'$\mu = 0.025$, $T=20$')
plt.plot([0, mu_point], [T_point, T_point], linestyle='--', color='black', linewidth=1)

plt.scatter(mu_point, T_point, color='black', zorder=5)
current_yticks = plt.yticks()[0]
new_yticks = np.sort(np.append(current_yticks, 25))
plt.xlabel(r'Switching rate per year ($\mu = \gamma$)', fontsize=13)
plt.ylabel(r'$T_{after}$ (years)', fontsize=13)
plt.ylim(0,90) 
plt.xlim(0, 0.1)

plt.tick_params(axis='both', labelsize=11)
plt.legend(loc='upper right', fontsize=11)


plt.grid(True, linestyle='--', alpha=0.5)
plt.savefig('/Users/finnkane/Desktop/Plots/mu_regime.png', dpi=300, bbox_inches="tight")

plt.show()
