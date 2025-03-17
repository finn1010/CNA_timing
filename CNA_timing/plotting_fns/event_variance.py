import arviz as az
import numpy as np
import matplotlib.pyplot as plt
import glob

def get_nc_files(folder):
    return sorted(glob.glob(f"{folder}/*.nc"))

folder_path_loh = '/Volumes/files/az/pool/variance/loh'
folder_path_tri = '/Volumes/files/az/pool/variance/tri'
folder_path_tet = '/Volumes/files/az/pool/variance/tet'

loh_az_files = get_nc_files(folder_path_loh)
tri_az_files = get_nc_files(folder_path_tri)
tet_az_files = get_nc_files(folder_path_tet)

loh_var = [float(file.split('/')[-1].replace('.nc', '')) for file in loh_az_files]
tri_var = [float(file.split('/')[-1].replace('.nc', '')) for file in tri_az_files]
tet_var = [float(file.split('/')[-1].replace('.nc', '')) for file in tet_az_files]

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

loh_sorted_idx = np.argsort(loh_var)
tri_sorted_idx = np.argsort(tri_var)
tet_sorted_idx = np.argsort(tet_var)

loh_var_sorted = np.array(loh_var)[loh_sorted_idx]
tri_var_sorted = np.array(tri_var)[tri_sorted_idx]
tet_var_sorted = np.array(tet_var)[tet_sorted_idx]

loh_ci_widths_sorted = [loh_ci_widths[i] for i in loh_sorted_idx]
tri_ci_widths_sorted = [tri_ci_widths[i] for i in tri_sorted_idx]
tet_ci_widths_sorted = [tet_ci_widths[i] for i in tet_sorted_idx]

offset = 0.25
loh_positions = np.arange(len(loh_var_sorted)) + 1 - offset
tri_positions = np.arange(len(tri_var_sorted)) + 1
    
tet_positions = np.arange(len(tet_var_sorted)) + 1 + offset

loh_color = 'lightblue'
tri_color = 'lightgreen'
tet_color = 'lightcoral'

custom_positions = np.array([1, 2, 3, 4])  
custom_labels = [12, 24, 29, 31]

offset = 0.19
loh_positions = custom_positions - offset  
tri_positions = custom_positions 
tet_positions = custom_positions + offset 

plt.figure(figsize=(10, 8))

plt.boxplot(loh_ci_widths_sorted, 
            positions=loh_positions, 
            widths=0.15, 
            patch_artist=True, 
            boxprops=dict(facecolor=loh_color, alpha=1, color='blue'),
            medianprops=dict(color='blue'),
            whiskerprops=dict(color='blue'),
            capprops=dict(color='blue'))

plt.boxplot(tri_ci_widths_sorted, 
            positions=tri_positions, 
            widths=0.15, 
            patch_artist=True, 
            boxprops=dict(facecolor=tri_color, alpha=1, color='green'),
            medianprops=dict(color='green'),
            whiskerprops=dict(color='green'),
            capprops=dict(color='green'))

plt.boxplot(tet_ci_widths_sorted, 
            positions=tet_positions, 
            widths=0.15, 
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
            color='black', label='Trisomy Mean', zorder=3, marker='^')
plt.scatter(tet_positions, [np.mean(ci) for ci in tet_ci_widths_sorted], 
            color='black', label='Tetraploidy Mean', zorder=3, marker='s')

plt.xticks(ticks=custom_positions, labels=custom_labels, fontsize=13)

plt.xlabel('IQR of event times (years)', fontsize=14.5)
plt.ylabel('95% Credible Interval (years)', fontsize=14.5)
plt.legend(fontsize=11, loc='upper right', bbox_to_anchor=(0.8357,1))
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig('variance.png', dpi=300, bbox_inches="tight")
plt.show()



folder_path_loh = '/Volumes/files/az/pool/variance/loh'

event_times_dict = {
    1: [22,30,27,33,16,27,28,18,16,25,20,24,29,24,15,31,14,13,20,10,29,13,17,23,17,29,28,12,15,29],
    4: [15,58,25,13,23,11,27,5,7,26,29,27,27,25,9,21,12,19,26,8,29,11,11,35,27,43,18,10,18,44],
    3: [11,24,11,5,29,54,23,31,23,8,29,48,28,57,11,60,17,9,24,7,20,18,31,13,5,14,60,15,10,45],
    4: [24,58,5,15,7,10,35,24,17,15,5,9,14,8,14,53,27,16,30,11,5,13,19,51,13,42,30,8,28,57],
    5: [39,8,48,5,29,54,9,8,20,5,5,29,5,57,5,64,12,9,41,4,37,23,15,37,27,56,10,5,28,47]
}

loh_az_files = get_nc_files(folder_path_loh)

loh_var = [float(file.split('/')[-1].replace('.nc', '')) for file in loh_az_files]

loh_az_fits = [az.from_netcdf(file) for file in loh_az_files]

loh_differences = []
for loh_val, az_fit in zip(loh_var, loh_az_fits):
    posterior_mean_t = az_fit.posterior['t'].mean(dim=["chain", "draw"]).values.flatten()
    
    true_event_times = np.array(event_times_dict.get(int(loh_val), []))

    loh_differences.append(posterior_mean_t - true_event_times)


loh_sorted_idx = np.argsort(loh_var)

loh_num_cna_sorted = np.array(loh_var)[loh_sorted_idx]

loh_differences_sorted = [loh_differences[i] for i in loh_sorted_idx]

fig, ax = plt.subplots(figsize=(10, 6))

ax.boxplot(loh_differences_sorted, vert=True, patch_artist=True)

ax.set_xticks(range(1, len(loh_num_cna_sorted) + 1))
ax.set_xticklabels(loh_num_cna_sorted, rotation=45)

ax.set_xlabel("LOH Variance")
ax.set_ylabel("Difference (Posterior Mean - True Event Time)")
ax.set_title("Box Plot of Differences Between Posterior Mean and True Event Times")

plt.show()
