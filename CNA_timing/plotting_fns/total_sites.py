import arviz as az
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.optimize as opt

def get_nc_files(folder):
    return sorted(glob.glob(f"{folder}/*.nc"))

folder_path_loh_sites = '/Volumes/files/az/pool/loh/num_sites'
folder_path_tri_sites = '/Volumes/files/az/pool/tri/num_sites'
folder_path_tet_sites = '/Volumes/files/az/pool/tet/num_sites'

loh_az_files_sites = get_nc_files(folder_path_loh_sites)
tri_az_files_sites = get_nc_files(folder_path_tri_sites)
tet_az_files_sites = get_nc_files(folder_path_tet_sites)

loh_num_sites = [float(file.split('/')[-1].replace('.nc', '')) for file in loh_az_files_sites]
tri_num_sites = [float(file.split('/')[-1].replace('.nc', '')) for file in tri_az_files_sites]
tet_num_sites = [float(file.split('/')[-1].replace('.nc', '')) for file in tet_az_files_sites]

loh_az_fits_sites = [az.from_netcdf(file) for file in loh_az_files_sites]
tri_az_fits_sites = [az.from_netcdf(file) for file in tri_az_files_sites]
tet_az_fits_sites = [az.from_netcdf(file) for file in tet_az_files_sites]

folder_path_loh_patients = '/Volumes/files/az/pool/loh/number_patients'
folder_path_tri_patients = '/Volumes/files/az/pool/tri/number_patients'
folder_path_tet_patients = '/Volumes/files/az/pool/tet/number_patients'

loh_az_files_patients = get_nc_files(folder_path_loh_patients)
tri_az_files_patients = get_nc_files(folder_path_tri_patients)
tet_az_files_patients = get_nc_files(folder_path_tet_patients)

loh_num_patients = [float(file.split('/')[-1].replace('.nc', '')) for file in loh_az_files_patients]
tri_num_patients = [float(file.split('/')[-1].replace('.nc', '')) for file in tri_az_files_patients]
tet_num_patients = [float(file.split('/')[-1].replace('.nc', '')) for file in tet_az_files_patients]

loh_az_fits_patients = [az.from_netcdf(file) for file in loh_az_files_patients]
tri_az_fits_patients = [az.from_netcdf(file) for file in tri_az_files_patients]
tet_az_fits_patients = [az.from_netcdf(file) for file in tet_az_files_patients]

def ci_interval(az_fit, hdi_prob=0.95):
    event_time_samples = az_fit.posterior['t']
    credible_intervals = az.hdi(event_time_samples, hdi_prob=hdi_prob)
    lower_bounds = credible_intervals['t'].sel(hdi="lower").values
    upper_bounds = credible_intervals['t'].sel(hdi="higher").values
    return upper_bounds - lower_bounds

loh_ci_sites = [ci_interval(az_fit) for az_fit in loh_az_fits_sites]
tri_ci_sites = [ci_interval(az_fit) for az_fit in tri_az_fits_sites]
tet_ci_sites = [ci_interval(az_fit) for az_fit in tet_az_fits_sites]

loh_ci_patients = [ci_interval(az_fit) for az_fit in loh_az_fits_patients]
tri_ci_patients = [ci_interval(az_fit) for az_fit in tri_az_fits_patients]
tet_ci_patients = [ci_interval(az_fit) for az_fit in tet_az_fits_patients]

loh_total_sites = [30 * num for num in loh_num_sites]
tri_total_sites = [30 * num for num in tri_num_sites]
tet_total_sites = [30 * num for num in tet_num_sites]

loh_total_patients = [500 * num for num in loh_num_patients]
tri_total_patients = [500 * num for num in tri_num_patients]
tet_total_patients = [500 * num for num in tet_num_patients]

loh_mean_ci_sites = [np.mean(ci) for ci in loh_ci_sites]
tri_mean_ci_sites = [np.mean(ci) for ci in tri_ci_sites]
tet_mean_ci_sites = [np.mean(ci) for ci in tet_ci_sites]

loh_mean_ci_patients = [np.mean(ci) for ci in loh_ci_patients]
tri_mean_ci_patients = [np.mean(ci) for ci in tri_ci_patients]
tet_mean_ci_patients = [np.mean(ci) for ci in tet_ci_patients]


def power_law(x, a, b):
    return a * np.power(x, b)

def plot_nonlinear_fit(x, y, color, linestyle, label):
    if len(x) > 1: 
        popt, _ = opt.curve_fit(power_law, x, y, maxfev=10000) 
        x_fit = np.linspace(min(x), max(x), 100)
        plt.plot(x_fit, power_law(x_fit, *popt), linestyle=linestyle, color=color, linewidth=2, label=label)

plt.figure(figsize=(10, 8))

plt.scatter(loh_total_sites, loh_mean_ci_sites, color='red', label='cn-LOH (sites)', marker='o')
plt.scatter(tri_total_sites, tri_mean_ci_sites, color='red', label='Trisomy (sites)', marker='^')
plt.scatter(tet_total_sites, tet_mean_ci_sites, color='red', label='Tetraploidy (sites)', marker='s')

plt.scatter(loh_total_patients, loh_mean_ci_patients, color='blue', label='cn-LOH (patients)', marker='o')
plt.scatter(tri_total_patients, tri_mean_ci_patients, color='blue', label='Trisomy (patients)', marker='^')
plt.scatter(tet_total_patients, tet_mean_ci_patients, color='blue', label='Tetraploidy (patients)', marker='s')

plot_nonlinear_fit(loh_total_sites, loh_mean_ci_sites, 'red', '-', 'cn-LOH (sites) Fit')
plot_nonlinear_fit(tri_total_sites, tri_mean_ci_sites, 'red', '--', 'Trisomy (sites) Fit')
plot_nonlinear_fit(tet_total_sites, tet_mean_ci_sites, 'red', ':', 'Tetraploidy (sites) Fit')

plot_nonlinear_fit(loh_total_patients, loh_mean_ci_patients, 'blue', '-', 'cn-LOH (patients) Fit')
plot_nonlinear_fit(tri_total_patients, tri_mean_ci_patients, 'blue', '--', 'Trisomy (patients) Fit')
plot_nonlinear_fit(tet_total_patients, tet_mean_ci_patients, 'blue', ':', 'Tetraploidy (patients) Fit')

plt.xlabel('Total number of sites', fontsize=14.5)
plt.ylabel('Mean 95% credible interval (years)', fontsize=14.5)
plt.legend(fontsize=12.5, loc='upper right')
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.savefig('/Users/finnkane/Desktop/Plots/sim/mean_ci_vs_total_sites.png', dpi=300, bbox_inches="tight")
plt.show()
