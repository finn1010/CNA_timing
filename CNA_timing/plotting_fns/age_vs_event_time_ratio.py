import arviz as az
import numpy as np
import matplotlib.pyplot as plt

loh_az_files = [
    '/Volumes/files/az/pool/loh/0.01/0.3/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.4/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.5/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.6/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.7/600.nc'
]
tri_az_files = [
    '/Volumes/files/az/pool/loh/0.01/0.3/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.4/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.5/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.6/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.7/600.nc'
]
tet_az_files = [
    '/Volumes/files/az/pool/loh/0.01/0.3/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.4/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.5/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.6/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.7/600.nc'
]

loh_az_files = [
    '/Volumes/files/az/pool/loh/0.01/0.3/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.4/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.5/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.6/600.nc',
    '/Volumes/files/az/pool/loh/0.01/0.7/600.nc'
]
tri_az_files = [
    '/Volumes/files/az/pool/tri/0.01/0.3.nc',
    '/Volumes/files/az/pool/tri/0.01/0.4.nc',
    '/Volumes/files/az/pool/tri/0.01/0.5.nc',
    '/Volumes/files/az/pool/tri/0.01/0.6.nc',
    '/Volumes/files/az/pool/tri/0.01/0.7.nc'
]
tet_az_files = [
    '/Volumes/files/az/pool/tet/0.01/0.3.nc',
    '/Volumes/files/az/pool/tet/0.01/0.4.nc',
    '/Volumes/files/az/pool/tet/0.01/0.5.nc',
    '/Volumes/files/az/pool/tet/0.01/0.6.nc',
    '/Volumes/files/az/pool/tet/0.01/0.7.nc'
]

def load_az_from_files(file_list):
    return [az.from_netcdf(file) for file in file_list]

def ci_interval(az_fit, hdi_prob=0.95):
    event_time_samples = az_fit.posterior['t']
    credible_intervals = az.hdi(event_time_samples, hdi_prob=hdi_prob)
    lower_bounds = credible_intervals['t'].sel(hdi="lower").values
    upper_bounds = credible_intervals['t'].sel(hdi="higher").values
    return upper_bounds - lower_bounds

def compute_ci_widths(az_fits):
    return [ci_interval(az_fit) for az_fit in az_fits]

loh_widths = compute_ci_widths(load_az_from_files(loh_az_files))
tri_widths = compute_ci_widths(load_az_from_files(tri_az_files))
tet_widths = compute_ci_widths(load_az_from_files(tet_az_files))
labels = ['0.3', '0.4', '0.5', '0.6', '0.7']
positions = np.arange(len(labels))
offset = 0.2
width = 0.15

plt.figure(figsize=(10, 8))

loh_color = 'lightblue'
tri_color = 'lightgreen'
tet_color = 'lightcoral'

plt.boxplot(loh_widths, 
            positions=positions - offset, 
            widths=width, 
            patch_artist=True, 
            boxprops=dict(facecolor=loh_color, alpha=1, color='blue'),
            medianprops=dict(color='blue'),
            whiskerprops=dict(color='blue'),
            capprops=dict(color='blue'))

plt.boxplot(tri_widths, 
            positions=positions, 
            widths=width, 
            patch_artist=True, 
            boxprops=dict(facecolor=tri_color, alpha=1, color='green'),
            medianprops=dict(color='green'),
            whiskerprops=dict(color='green'),
            capprops=dict(color='green'))

plt.boxplot(tet_widths, 
            positions=positions + offset, 
            widths=width, 
            patch_artist=True, 
            boxprops=dict(facecolor=tet_color, alpha=1, color='red'),
            medianprops=dict(color='red'),
            whiskerprops=dict(color='red'),
            capprops=dict(color='red'))

plt.scatter(positions - offset, [np.mean(ci) for ci in loh_widths], 
            color='black', label='LOH Mean', zorder=3, marker='o')
plt.scatter(positions, [np.mean(ci) for ci in tri_widths], 
            color='black', label='Trisomy Mean', zorder=3, marker='s')
plt.scatter(positions + offset, [np.mean(ci) for ci in tet_widths], 
            color='black', label='Tetraploidy Mean', zorder=3, marker='^')

plt.plot([], [], linestyle='-', linewidth=8, color=loh_color, alpha=1, label='LOH')
plt.plot([], [], linestyle='-', linewidth=8, color=tri_color, alpha=1, label='Trisomy')
plt.plot([], [], linestyle='-', linewidth=8, color=tet_color, alpha=1, label='Tetraploidy')

plt.xticks(ticks=positions, labels=labels, fontsize=13)
plt.xlabel('Event Time to Patient Age Ratio', fontsize=14.5)
plt.ylabel('95% Credible Interval Width', fontsize=14.5)

plt.legend(fontsize=12.5, loc='upper right')

plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig('/Users/finnkane/Desktop/Plots/sim/ratio.png', dpi=300, bbox_inches="tight")

plt.show()
