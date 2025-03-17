import numpy as np
import matplotlib.pyplot as plt
import arviz as az
import glob

def get_nc_files(folder):
    return sorted(glob.glob(f"{folder}/*.nc"))

folder_path_loh = '/Volumes/T7/ICR1/k4'
folder_path_tri = '/Volumes/T7/ICR1/k5'
folder_path_tet = '/Volumes/T7/ICR1/k6'

loh_az_files = get_nc_files(folder_path_loh)
tri_az_files = get_nc_files(folder_path_tri)
tet_az_files = get_nc_files(folder_path_tet)

loh_num_cna = [float(file.split('/')[-1].replace('.nc', '')) for file in loh_az_files]
tri_num_cna = [float(file.split('/')[-1].replace('.nc', '')) for file in tri_az_files]
tet_num_cna = [float(file.split('/')[-1].replace('.nc', '')) for file in tet_az_files]

loh_az_fits = [az.from_netcdf(file) for file in loh_az_files]
tri_az_fits = [az.from_netcdf(file) for file in tri_az_files]
tet_az_fits = [az.from_netcdf(file) for file in tet_az_files]

def compute_posterior_mean(az_fit):
    return az_fit.posterior['t'].mean().values

loh_posterior_means = np.array([compute_posterior_mean(az_fit) for az_fit in loh_az_fits])
tri_posterior_means = np.array([compute_posterior_mean(az_fit) for az_fit in tri_az_fits])
tet_posterior_means = np.array([compute_posterior_mean(az_fit) for az_fit in tet_az_fits])

loh_sorted_idx = np.argsort(loh_num_cna)
tri_sorted_idx = np.argsort(tri_num_cna)
tet_sorted_idx = np.argsort(tet_num_cna)

loh_num_cna_sorted = np.array(loh_num_cna)[loh_sorted_idx]
tri_num_cna_sorted = np.array(tri_num_cna)[tri_sorted_idx]
tet_num_cna_sorted = np.array(tet_num_cna)[tet_sorted_idx]

loh_posterior_means_sorted = loh_posterior_means[loh_sorted_idx]
tri_posterior_means_sorted = tri_posterior_means[tri_sorted_idx]
tet_posterior_means_sorted = tet_posterior_means[tet_sorted_idx]

def compute_rate_of_change(x, y):
    return np.gradient(y, x)

loh_rate_of_change = compute_rate_of_change(loh_num_cna_sorted, loh_posterior_means_sorted)
tri_rate_of_change = compute_rate_of_change(tri_num_cna_sorted, tri_posterior_means_sorted)
tet_rate_of_change = compute_rate_of_change(tet_num_cna_sorted, tet_posterior_means_sorted)

plt.figure(figsize=(10, 8))
plt.plot(loh_num_cna_sorted, loh_rate_of_change, '-o', color='blue', label='cn-LOH')
plt.plot(tri_num_cna_sorted, tri_rate_of_change, '-^', color='green', label='Trisomy')
plt.plot(tet_num_cna_sorted, tet_rate_of_change, '-s', color='red', label='Tetraploidy')

plt.xlabel('Kappa ($\kappa$)', fontsize=14.5)
plt.ylabel('Rate of change of mean credible interval width (years/$\kappa$)', fontsize=14.5)
plt.legend(fontsize=12.5, loc='upper right')
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig('f.png', dpi=300, bbox_inches="tight")
plt.show()
