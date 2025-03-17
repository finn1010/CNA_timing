import arviz as az
import numpy as np
import matplotlib.pyplot as plt
import glob

# def get_nc_files(folder):
#     return sorted(glob.glob(f"{folder}/*.nc"))
# folder_path_loh = '/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)'
# folder_path_tri = '/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)'
# folder_path_tet = '/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)'
# loh_az_files = get_nc_files(folder_path)
# mu_values = [float(file.split('/')[-1].replace('.nc', '')) for file in loh_az_files]

# az_fits = [az.from_netcdf(file) for file in loh_az_files]

# def ci_interval(az_fit, hdi_prob=0.95):
#     event_time_samples = az_fit.posterior['t']
#     credible_intervals = az.hdi(event_time_samples, hdi_prob=hdi_prob)
#     lower_bounds = credible_intervals['t'].sel(hdi="lower").values
#     upper_bounds = credible_intervals['t'].sel(hdi="higher").values
#     return upper_bounds - lower_bounds

# ci_widths = [ci_interval(az_fit) for az_fit in az_fits]
# mean_ci_widths = [np.mean(ci) for ci in ci_widths]

# mu_positions = np.arange(len(mu_values)) + 1 

# plt.figure(figsize=(8, 6))
# plt.boxplot(ci_widths, positions=mu_positions, widths=0.5) 

# plt.scatter(mu_positions, mean_ci_widths, color='red', label='Mean CI Width', zorder=3)

# plt.xticks(mu_positions, labels=[f'{mu:.5f}' for mu in mu_values], rotation=45, ha="right")
# plt.xlabel('Switching Rate $\mu$')
# plt.ylabel('Credible Interval Width')
# plt.title('Box Plot of CI Widths vs. $\mu$')
# plt.legend()
# plt.grid(True, linestyle="--", alpha=0.5)
# plt.show()

def get_nc_files(folder):
    return sorted(glob.glob(f"{folder}/*.nc"))

def load_az_files(folder_path):
    az_files = get_nc_files(folder_path)
    mu_values = [float(file.split('/')[-1].replace('.nc', '')) for file in az_files]
    az_fits = [az.from_netcdf(file) for file in az_files]
    return mu_values, az_fits



def get_nc_files(folder):
    return sorted(glob.glob(f"{folder}/*.nc"))

def load_az_files(folder_path):
    az_files = get_nc_files(folder_path)
    mu_values = [float(file.split('/')[-1].replace('.nc', '')) for file in az_files]
    az_fits = [az.from_netcdf(file) for file in az_files]
    return mu_values, az_fits

def ci_interval(az_fit, hdi_prob=0.95):
    event_time_samples = az_fit.posterior['t']
    credible_intervals = az.hdi(event_time_samples, hdi_prob=hdi_prob)
    lower_bounds = credible_intervals['t'].sel(hdi="lower").values
    upper_bounds = credible_intervals['t'].sel(hdi="higher").values
    return upper_bounds - lower_bounds

def compute_ci_widths(az_fits):
    ci_widths = [ci_interval(az_fit) for az_fit in az_fits]
    return ci_widths

folder_path_loh = '/Volumes/files/az/pool/loh/mu_regime/norm(0,0.01)'
folder_path_tri = '/Volumes/files/az/pool/tri/mu_regime/norm(0,0.01)'
folder_path_tet = '/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)'

mu_values_loh, az_fits_loh = load_az_files(folder_path_loh)
mu_values_tri, az_fits_tri = load_az_files(folder_path_tri)
mu_values_tet, az_fits_tet = load_az_files(folder_path_tet)

ci_widths_loh = compute_ci_widths(az_fits_loh)
ci_widths_tri = compute_ci_widths(az_fits_tri)
ci_widths_tet = compute_ci_widths(az_fits_tet)

all_mu_values = sorted(set(mu_values_loh + mu_values_tri + mu_values_tet))

def align_data(all_mu, mu_values, ci_widths):
    ci_aligned = []
    ci_means = []
    for mu in all_mu:
        if mu in mu_values:
            ci_values = ci_widths[mu_values.index(mu)]
            ci_aligned.append(ci_values)
            ci_means.append(np.mean(ci_values))
        else:
            ci_aligned.append(np.array([]))
            ci_means.append(None)
    return ci_aligned, ci_means

loh_data, loh_means = align_data(all_mu_values, mu_values_loh, ci_widths_loh)
tri_data, tri_means = align_data(all_mu_values, mu_values_tri, ci_widths_tri)
tet_data, tet_means = align_data(all_mu_values, mu_values_tet, ci_widths_tet)

plt.figure(figsize=(15, 7))

width = 0.2
spacing = 0.25
positions_loh = np.arange(len(all_mu_values)) - spacing
positions_tri = np.arange(len(all_mu_values))
positions_tet = np.arange(len(all_mu_values)) + spacing

plt.boxplot(
    loh_data, 
    positions=positions_loh, 
    widths=width, 
    patch_artist=True, 
    boxprops=dict(facecolor='lightblue', alpha=1), 
    medianprops=dict(color='black'), 
    zorder=2, 
    manage_ticks=False
)
plt.boxplot(
    tri_data, 
    positions=positions_tri, 
    widths=width, 
    patch_artist=True, 
    boxprops=dict(facecolor='lightgreen', alpha=1), 
    medianprops=dict(color='black'), 
    zorder=2, 
    manage_ticks=False
)
plt.boxplot(
    tet_data, 
    positions=positions_tet, 
    widths=width, 
    patch_artist=True, 
    boxprops=dict(facecolor='lightcoral', alpha=1), 
    medianprops=dict(color='black'), 
    zorder=2, 
    manage_ticks=False
)
plt.plot([], [], linestyle='-', linewidth=8, color='lightblue', alpha=1, label='cn-LOH')
plt.plot([], [], linestyle='-', linewidth=8, color='lightgreen', alpha=1, label='Trisomy')
plt.plot([], [], linestyle='-', linewidth=8, color='lightcoral', alpha=1, label='Tetraploidy')

plt.scatter(positions_loh, loh_means, color='black', s=25, label='LOH Mean', zorder=3)
plt.scatter(positions_tri, tri_means, color='black', s=25, label='TRI Mean', marker='s', zorder=3)
plt.scatter(positions_tet, tet_means, color='black', s=25, label='TET Mean', marker='^', zorder=3)



plt.xticks(
    np.arange(len(all_mu_values)),
    labels=[f'{mu:.5f}' for mu in all_mu_values],
)

plt.xlabel('Switching Rate (per year)',fontsize=13)
plt.ylabel('95% Credible Interval (years)',fontsize=13)
plt.legend(loc='upper center', bbox_to_anchor=(0.65, 1.0),fontsize=12)
plt.tick_params(axis='both',labelsize=12)
plt.tight_layout()
plt.savefig('/Users/finnkane/Desktop/Plots/mu_regime_vals.png', dpi=300, bbox_inches="tight")

plt.show()

