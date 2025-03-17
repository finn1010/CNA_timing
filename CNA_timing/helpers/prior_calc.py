import numpy as np
import matplotlib.pyplot as plt

P = 1000000

mu_mean = np.random.normal(-5.1, 0.35, P)
mu_sdd = np.random.normal(0, 0.2, P)
mu_raw = np.random.normal(0, 0.1, P)

mu_p = np.exp(mu_mean + mu_sdd * mu_raw)

plt.figure(figsize=(8, 6))
plt.hist(mu_p, bins=100, density=True, alpha=0.75, edgecolor='black')
plt.yscale('log')
plt.xlabel('$\mu[p]$', fontsize=13)
plt.ylabel('Density', fontsize=13)
plt.title('Distribution of $\mu[p]$', fontsize=13)
plt.show()