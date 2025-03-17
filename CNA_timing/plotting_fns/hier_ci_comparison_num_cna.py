import arviz as az
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit

def get_nc_files(folder):
    return sorted(glob.glob(f"{folder}/*.nc"))

folder_path_loh = '/Volumes/files/az/hier/loh'
folder_path_tri = '/Volumes/files/az/hier/tri'
folder_path_tet = '/Volumes/files/az/hier/tet'

loh_az_files = get_nc_files(folder_path_loh)
tri_az_files = get_nc_files(folder_path_tri)
tet_az_files = get_nc_files(folder_path_tet)

loh_num_sites = [float(file.split('/')[-1].replace('.nc', '')) for file in loh_az_files]
tri_num_sites = [float(file.split('/')[-1].replace('.nc', '')) for file in tri_az_files]
tet_num_sites = [float(file.split('/')[-1].replace('.nc', '')) for file in tet_az_files]

loh_az_fits = [az.from_netcdf(file) for file in loh_az_files]
tri_az_fits = [az.from_netcdf(file) for file in tri_az_files]
tet_az_fits = [az.from_netcdf(file) for file in tet_az_files]

def ci_interval(az_fit, hdi_prob=0.95):
    event_time_samples = az_fit.posterior['t']
    credible_intervals = az.hdi(event_time_samples, hdi_prob=hdi_prob)
    lower_bounds = credible_intervals['t'].sel(hdi="lower").values
    upper_bounds = credible_intervals['t'].sel(hdi="higher").values
    return upper_bounds - lower_bounds

loh_ci_widths = [ci_interval(az_fit) for az_fit in loh_az_fits]
tri_ci_widths = [ci_interval(az_fit) for az_fit in tri_az_fits]
tet_ci_widths = [ci_interval(az_fit) for az_fit in tet_az_fits]

loh_sorted_idx = np.argsort(loh_num_sites)
tri_sorted_idx = np.argsort(tri_num_sites)
tet_sorted_idx = np.argsort(tet_num_sites)

loh_num_sites_sorted = np.array(loh_num_sites)[loh_sorted_idx]
tri_num_sites_sorted = np.array(tri_num_sites)[tri_sorted_idx]
tet_num_sites_sorted = np.array(tet_num_sites)[tet_sorted_idx]

loh_ci_widths_sorted = [loh_ci_widths[i] for i in loh_sorted_idx]
tri_ci_widths_sorted = [tri_ci_widths[i] for i in tri_sorted_idx]
tet_ci_widths_sorted = [tet_ci_widths[i] for i in tet_sorted_idx]

offset = 0.25
loh_positions = np.arange(len(loh_num_sites_sorted)) + 1 - offset
tri_positions = np.arange(len(tri_num_sites_sorted)) + 1
tet_positions = np.arange(len(tet_num_sites_sorted)) + 1 + offset

plt.figure(figsize=(10, 8))


loh_color = 'lightblue'
tri_color = 'lightgreen'
tet_color = 'lightcoral'

plt.boxplot(loh_ci_widths_sorted, 
            positions=loh_positions, 
            widths=0.2, 
            patch_artist=True, 
            boxprops=dict(facecolor=loh_color, alpha=1, color='blue'),
            medianprops=dict(color='blue'),
            whiskerprops=dict(color='blue'),
            capprops=dict(color='blue'))

plt.boxplot(tri_ci_widths_sorted, 
            positions=tri_positions, 
            widths=0.2, 
            patch_artist=True, 
            boxprops=dict(facecolor=tri_color, alpha=1, color='green'),
            medianprops=dict(color='green'),
            whiskerprops=dict(color='green'),
            capprops=dict(color='green'))

plt.boxplot(tet_ci_widths_sorted, 
            positions=tet_positions, 
            widths=0.2, 
            patch_artist=True, 
            boxprops=dict(facecolor=tet_color, alpha=1, color='red'),
            medianprops=dict(color='red'),
            whiskerprops=dict(color='red'),
            capprops=dict(color='red'))

plt.plot([], [], linestyle='-', linewidth=8, color=loh_color, alpha=1, label='cn-LOH')
plt.plot([], [], linestyle='-', linewidth=8, color=tri_color, alpha=1, label='Trisomy')
plt.plot([], [], linestyle='-', linewidth=8, color=tet_color, alpha=1, label='Tetraploidy')

plt.scatter(loh_positions, [np.mean(ci) for ci in loh_ci_widths_sorted], 
            color='black', label='LOH Mean', zorder=3, marker='o')
plt.scatter(tri_positions, [np.mean(ci) for ci in tri_ci_widths_sorted], 
            color='black', label='Trisomy Mean', zorder=3, marker='s')
plt.scatter(tet_positions, [np.mean(ci) for ci in tet_ci_widths_sorted], 
            color='black', label='Tetraploidy Mean', zorder=3, marker='^')



tri_total_sites = [30 * num for num in tri_num_sites_sorted]

plt.xticks(
    ticks=tri_positions.astype(int),
    labels=[str(i + 1) for i in range(len(tri_positions))], 
    fontsize=13
)

plt.xlabel('Number CNA events per patient', fontsize=14.5)

plt.ylim(0,70)
# plt.xlabel('Number of fCpG Sites', fontsize=14.5)
plt.ylabel('95% Credible Interval (years)', fontsize=14.5)
plt.legend(fontsize=12.5, loc='upper right')
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()


plt.savefig('hier_num_cna', dpi=300, bbox_inches="tight")
plt.show()
