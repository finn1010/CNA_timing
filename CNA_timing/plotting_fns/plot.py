from scipy import linalg
import numpy as np
import matplotlib.pyplot as plt
from helpers.colours import pallet_dip,pallet_tri
import arviz as az
from scipy.stats import linregress
import seaborn as sns
import pandas as pd
import scipy.stats as stats
from scipy.stats import ttest_ind


def stat_creds(az_fit, event_times, patient_ages, patient_idx, event_types):
    event_time_samples = az_fit.posterior['t']
    cred_int = az.hdi(event_time_samples, hdi_prob=0.95)

    lower_bounds = cred_int['t'].sel(hdi="lower").values
    upper_bounds = cred_int['t'].sel(hdi="higher").values
    event_times = np.array(event_times)
    patient_ages = np.array(patient_ages)
    patient_idx = np.array(patient_idx)
    event_types = np.array(event_types)

    credible_interval_sizes = upper_bounds - lower_bounds

    trisomy_events = {"chr8", "chr12", "chr19", "chr15", "chr17", "chr2", "chr1"}
    cn_loh_events = {"chr11", "chr13", "chr20", "chr6"}

    trisomy_sizes = [size for size, event in zip(credible_interval_sizes, event_types) if event in trisomy_events]
    cn_loh_sizes = [size for size, event in zip(credible_interval_sizes, event_types) if event in cn_loh_events]

    trisomy_mean = np.median(trisomy_sizes)
    cn_loh_mean = np.median(cn_loh_sizes)

    t_stat, p_value = stats.ttest_ind(trisomy_sizes, cn_loh_sizes, equal_var=False)

    return {
        "tri_mean_int": trisomy_mean,
        "loh_mean_int": cn_loh_mean,
        "t_stat": t_stat,
        "p_val": p_value
    }



def plot_event_time_uncertainty_hier_switched(fig_name, az_fit, event_times, patient_ages, patient_idx, event_types, extracted_data, df_final):
    event_time_samples = az_fit.posterior['t']
    cred_int = az.hdi(event_time_samples, hdi_prob=0.95)

    lower_bounds = cred_int['t'].sel(hdi="lower").values
    upper_bounds = cred_int['t'].sel(hdi="higher").values
    
    event_times = np.array(event_times)
    patient_ages = np.array(patient_ages)
    patient_idx = np.array(patient_idx)
    
    patient_idx_map = {pid: i for i, pid in enumerate(df_final['index'])}
    patient_idx_mapped = np.array([patient_idx_map.get(pid, -1) for pid in extracted_data['patient_idx']])
    event_types = np.array(extracted_data['cna_types'])
    chromosomes = np.array(df_final['chr'])  

    ci_sizes = upper_bounds - lower_bounds
    
    ci_sizes_type_4 = ci_sizes[event_types == 4]
    ci_sizes_type_5 = ci_sizes[event_types == 5]
    print(len(ci_sizes_type_4))
    print(len(ci_sizes_type_5))

    mean_ci_4 = np.mean(ci_sizes_type_4) if len(ci_sizes_type_4) > 0 else np.nan
    mean_ci_5 = np.mean(ci_sizes_type_5) if len(ci_sizes_type_5) > 0 else np.nan
    median_ci_4 = np.median(ci_sizes_type_4) if len(ci_sizes_type_4) > 0 else np.nan
    median_ci_5 = np.median(ci_sizes_type_5) if len(ci_sizes_type_5) > 0 else np.nan
    
    if len(ci_sizes_type_4) > 1 and len(ci_sizes_type_5) > 1:
        t_stat, p_value = ttest_ind(ci_sizes_type_4, ci_sizes_type_5, equal_var=False)
    else:
        t_stat, p_value = np.nan, np.nan 

    print(f"Mean Credible Interval Size for CNA Type 4: {mean_ci_4:.4f}")
    print(f"Mean Credible Interval Size for CNA Type 5: {mean_ci_5:.4f}")
    print(f"Mean Credible Interval Size for CNA Type 4: {median_ci_4:.4f}")
    print(f"Mean Credible Interval Size for CNA Type 5: {median_ci_5:.4f}")
    print(f"T-Statistic: {t_stat:.4f}, P-Value: {p_value:.4e}")

    marker_options = ["o", "s", "D", "^", "v", "p", "*", "<", "d", "H", "d", ">"]
    unique_chromosomes = np.unique(chromosomes)

    chromosome_marker_map = {chrom: marker_options[i % len(marker_options)] for i, chrom in enumerate(unique_chromosomes)}

    color_options = plt.cm.get_cmap("tab10", len(np.unique(event_types)))
    event_color_map = {event: color_options(i) for i, event in enumerate(np.unique(event_types))}
    event_color_map[4] = "red"
    event_color_map[5] = "blue"

    trisomy_events = {"chr8", "chr12", "chr19", "chr15", "chr17", "chr2", "chr1"}
    cn_loh_events = {"chr11", "chr13", "chr20", "chr6", "chr1", "chr17"}  

    unique_patients = np.unique(patient_idx_mapped[patient_idx_mapped != -1])
    sorted_patient_indices = unique_patients[np.argsort([patient_ages[patient_idx_mapped == p][0] for p in unique_patients])]

    patient_idx_sorted_map = {patient: i + 1 for i, patient in enumerate(sorted_patient_indices)}
    patient_idx_sorted = np.array([patient_idx_sorted_map.get(p, -1) for p in patient_idx_mapped if p != -1])

    plt.figure(figsize=(11, 7))

    unique_patients_sorted, patient_counts = np.unique(patient_idx_sorted, return_counts=True)
    sorted_offsets = {patient: list(np.linspace(-0.15, 0.15, count)) if count > 1 else [0]
                      for patient, count in zip(unique_patients_sorted, patient_counts)}

    offsets = [sorted_offsets[patient].pop(0) for patient in patient_idx_sorted]

    trisomy_legend_handles = []
    cn_loh_legend_handles = []

    for chrom in sorted(unique_chromosomes):
        marker = chromosome_marker_map[chrom]
        handle = plt.Line2D([0], [0], marker=marker, color='black', markersize=5, linestyle='None', label=chrom)
        
        if chrom in trisomy_events:
            trisomy_legend_handles.append(handle)
        if chrom in cn_loh_events:
            cn_loh_legend_handles.append(handle)

    for lower, upper, event_time, index, offset, event_type, chrom in zip(
        lower_bounds, upper_bounds, event_times, patient_idx_sorted, offsets, event_types, chromosomes
    ):
        color = event_color_map.get(event_type, "black")
        marker = chromosome_marker_map.get(chrom, "o")  
        plt.plot([index + offset, index + offset], [lower, upper], color=color, alpha=0.7, linewidth=2)
        plt.plot(index + offset, event_time, marker=marker, markersize=5, color="black", alpha=0.8)

    patient_age_handles = []
    for patient in unique_patients_sorted:
        var = patient_idx_sorted == patient 
        patient_age = patient_ages[var][0] 
        min_x = patient - 0.15
        max_x = patient + 0.15 * (np.sum(var) - 1) 
        plt.plot([min_x, max_x], [patient_age, patient_age], 'g-', linewidth=2.5, alpha=0.7)

    patient_age_handles.append(plt.Line2D([0], [0], color='blue', linewidth=2.5, label="95% Credible Interval (Trisomy)"))
    patient_age_handles.append(plt.Line2D([0], [0], color='red', linewidth=2.5, label="95% Credible Interval (cn-LOH)"))
    patient_age_handles.append(plt.Line2D([0], [0], color='green', linewidth=2.5, label="Patient age"))

    trisomy_legend = plt.legend(handles=trisomy_legend_handles, fontsize=12, loc='upper left', 
                                title="Trisomy:", title_fontsize=13)

    cn_loh_legend = plt.legend(handles=cn_loh_legend_handles, fontsize=12, loc='upper left', bbox_to_anchor=(0.11, 1), 
                               title="cn-LOH:", title_fontsize=13)

    patient_age_legend = plt.legend(handles=patient_age_handles, fontsize=12, loc='upper left', bbox_to_anchor=(0.22, 1))

    plt.gca().add_artist(trisomy_legend)
    plt.gca().add_artist(cn_loh_legend)

    plt.yticks(fontsize=14)
    plt.xticks([])
    plt.ylabel("Time (years)", fontsize=18)

    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_name, dpi=300, bbox_inches="tight")
    plt.show()



def plot_prob_dist(time_points, probabilities, plt_title):
    plt.figure(figsize=(10, 6))
    for state in range(probabilities.shape[1]):
        plt.plot(time_points, probabilities[:, state], label=f'State {state + 1}')
    plt.xlabel('Time (years)',fontsize=13)  
    plt.ylabel('Probability',fontsize=13)
    plt.legend()
    plt.show()

def hist_plot(beta_vals_before, beta_vals_after, event_type, pre_time,post_time,fig_name):
    plt.hist(beta_vals_before, bins=30, edgecolor='black', alpha=0.5, label=f'Time before event ({pre_time} years)', density=True, align='mid')
    plt.hist(beta_vals_after, bins=30, edgecolor='black', alpha=0.5, label=f'Time after event ({post_time} years)', density=True, align='mid')
    
    plt.xlabel('Fraction methylated',fontsize=14)
    plt.ylabel('Probability Density',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='upper right', fontsize=12, bbox_to_anchor=(0.9,1))
    plt.savefig(f'{fig_name}.pdf', format='pdf', dpi=300)
    plt.show()



def plot_ppc(az_fit, file_path):
    ax = az.plot_ppc(
        data=az_fit,
        data_pairs={"y": "y_rep"},
        kind="kde",
    )

    lines = ax.get_lines() 
    for line in lines:
        if "Posterior predictive mean" in line.get_label():
            x_data, y_data = line.get_data()
            var = x_data > 0.001  
            line.set_data(x_data[var], y_data[var])

    ax.set_xlabel("Fraction of alleles methylated", fontsize=13)
    ax.set_ylabel("Probability Density", fontsize=13)
    
    ax.tick_params(axis='y', which='both', direction='in', length=6, width=1, labelsize=12)

    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.tight_layout()
    plt.show()


def plot_pair(az_fit,patient_num):
    t_num1 = patient_num * 3
    t_num2 = patient_num * 3 + 1
    t_num3 = patient_num * 3 + 2

    az.plot_pair(
        az_fit,
        var_names=["gamma", "mu", "t"], 
        coords={"gamma_dim_0": [patient_num], "mu_dim_0": [patient_num], "t_dim_0": [t_num1, t_num2, t_num3]}, 
    )
    plt.show()


def plot_pair_pool(fig_name, az_fit, patient_num):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    pairs = [("gamma_raw", "mu"), ("gamma_raw", "t"), ("mu", "t")]
    t = [
        'x',
        'x',
        'x'
    ]

    for ax, (x_var, y_var), title in zip(axes, pairs, t):
        az.plot_pair(
            az_fit,
            var_names=[x_var, y_var],
            coords={"t_dim_0": [patient_num]},
            ax=ax,
        )
        ax.set_xlabel(
            "gamma_raw" if x_var == "gamma_raw" else 
            "Mu (years)" if x_var == "mu" else 
            "Time (years)",
            fontsize=14
        )
        ax.set_ylabel(
            "gamma_raw" if y_var == "gamma_raw" else 
            "Mu (years)" if y_var == "mu" else 
            "Time (years)",
            fontsize=14
        )
        ax.tick_params(axis='both', labelsize=12)

    plt.tight_layout()
    plt.savefig(fig_name, dpi=300, bbox_inches="tight")

    plt.show()




def plot_mu_distribution(fig_name, az_fit, true_mu=0.01, J=4000):
    prior_mu = np.random.normal(0, 0.01, J)
    posterior_mu = az_fit.posterior['mu'].values.flatten()
    plt.figure(figsize=(7, 5))
    plt.hist(prior_mu, bins=500, histtype="step", color='b', label='Prior', density=False, linewidth=1)
    plt.hist(posterior_mu, bins=200, histtype="step", color='r', label='Posterior', density=False, linewidth=1)

    plt.axvline(true_mu, color='red', linestyle='--', linewidth=1, label=r"True $\mu$")
    plt.xlim(0,0.025)
    plt.xlabel(r"Rate of methylation ($\mu$)", fontsize=14)
    plt.ylabel("Density", fontsize=14)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(fig_name, dpi=300, bbox_inches="tight")

    plt.show()

def plot_distributions(az_fit, mu_values, gamma_values, J=4000):
    mu = np.abs(np.random.normal(0, 0.1, J))
    gamma = np.abs(np.random.normal(1, 0.5, J) * mu)
    
    samples = {
        "mu": mu,
        "gamma": gamma,
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    posterior_mu = az_fit.posterior["mu"].values.flatten()
    posterior_gamma = az_fit.posterior["gamma"].values.flatten()

    axes[0].hist(samples["mu"], bins=200, histtype="step", color='b', label='Prior', density=True, linewidth=1)
    axes[0].hist(posterior_mu, bins=200, histtype="step", color='r', label='Posterior', density=True, linewidth=1)
    
    for mu in mu_values:
        axes[0].axvline(mu, color='red', linestyle='--', linewidth=0.5)
    
    axes[0].set_xlim(0, 0.09)
    axes[0].set_xlabel(r"Rate of methylation ($\mu$)",fontsize=13)
    axes[0].set_ylabel("Density",fontsize=13)
    axes[0].legend(labels=["Prior", "Posterior", r"True $\mu$ values"], fontsize=11)

    axes[1].hist(samples["gamma"], bins=400, histtype="step", color='b', label='Prior', density=True, linewidth=1)
    axes[1].hist(posterior_gamma, bins=400, histtype="step", color='r', label='Posterior', density=True, linewidth=1)
    
    for gamma in gamma_values:
        axes[1].axvline(gamma, color='red', linestyle='--', linewidth=0.5)
    
    axes[1].set_xlim(0, 0.09)
    axes[1].set_xlabel(r"Rate of demethylation ($\gamma$)",fontsize=13)
    axes[1].set_ylabel("Density",fontsize=13)

    axes[1].legend(labels=["Prior", "Posterior", r"True $\gamma$ values"], fontsize=11)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_event_time_uncertainty_hier(fig_name, az_fit, event_times, patient_ages, group_size=1):
    event_time_samples = az_fit.posterior['t']
    cred_int = az.hdi(event_time_samples, hdi_prob=0.95)

    lower_bounds = cred_int['t'].sel(hdi="lower").values
    upper_bounds = cred_int['t'].sel(hdi="higher").values
    event_times = np.array(event_times)
    patient_ages = np.array(patient_ages)

    plt.figure(figsize=(10, 7))

    num_patients = len(lower_bounds) // group_size 
    patient_indices = np.repeat(np.arange(1, num_patients + 1), group_size)[:len(lower_bounds)]

    offsets = np.linspace(-0.15, 0.15, group_size) if group_size > 1 else [0]

    for i, (lower, upper, event_time, patient_age, index, offset) in enumerate(
        zip(lower_bounds, upper_bounds, event_times, patient_ages, patient_indices, 
            np.tile(offsets, len(patient_indices) // group_size + 1)[:len(lower_bounds)])
    ):
        plt.plot([lower, upper], [index + offset, index + offset], 'b-', alpha=0.7, linewidth=2)

        plt.plot(event_time, index + offset, 'rx', markersize=9.5, label="True event time" if i == 0 else None)

    unique_patient_indices = np.unique(patient_indices)
    for patient_idx in unique_patient_indices:
        var = patient_indices == patient_idx
        patient_age = patient_ages[var][0]  
        min_y = patient_idx - 0.15
        max_y = patient_idx + 0.15 * (group_size - 1)
        
        plt.plot([patient_age, patient_age], [min_y, max_y], 'g-', linewidth=2.5, alpha=0.7, label="Patient Age" if patient_idx == 1 else None)

    plt.xticks(fontsize=14)

    plt.yticks(unique_patient_indices, fontsize=14)

    plt.xlabel("Time (years)", fontsize=18)
    plt.ylabel("Patient Index", fontsize=18)
    plt.legend(fontsize=14, loc='lower right', bbox_to_anchor=(1, 0.22))
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_name, dpi=300, bbox_inches="tight")

    plt.show()

def plot_event_time_uncertainty(fig_name, az_fit, event_times, patient_ages):
    event_time_samples = az_fit.posterior['t']
    event_time_array = event_time_samples.values
    cred_int = az.hdi(event_time_samples, hdi_prob=0.95)

    lower_bounds = cred_int['t'].sel(hdi="lower").values
    upper_bounds = cred_int['t'].sel(hdi="higher").values
    event_times = np.array(event_times)
    patient_ages = np.array(patient_ages)

    credible_interval_sizes = upper_bounds - lower_bounds
    avg_credible_interval = np.median(credible_interval_sizes)
    print(f"Average credible interval size: {avg_credible_interval:.2f} years")

    plt.figure(figsize=(10, 7))

    num_patients = len(lower_bounds)

    for i, (lower, upper, event_time, patient_age) in enumerate(zip(lower_bounds, upper_bounds, event_times, patient_ages)):
        patient_index = i + 1  

        plt.plot([lower, upper], [patient_index, patient_index], 'b-', alpha=0.7, linewidth=2)  
        plt.plot(event_time, patient_index, 'rx', markersize=9.5, label="True event time" if i == 0 else None) 
        plt.plot(patient_age, patient_index, 'go', markersize=6.5, label="Patient age" if i == 0 else None)  

    plt.xticks(fontsize=14)

    yticks = np.arange(0, num_patients + 5, 5)
    plt.yticks(yticks, fontsize=14)

    plt.xlabel("Time (years)", fontsize=18)
    plt.ylabel("Patient Index", fontsize=18)
    plt.legend(fontsize=14, loc='lower right', bbox_to_anchor=(1, 0.22))
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_name, dpi=300, bbox_inches="tight")

    plt.show()



def plot_mu_uncertainty(fig_name, az_fit, mu_values, patient_ages):
    mu_samples = az_fit.posterior['mu']
    cred_int = az.hdi(mu_samples, hdi_prob=0.95)

    lower_bounds = cred_int['mu'].sel(hdi="lower").values
    upper_bounds = cred_int['mu'].sel(hdi="higher").values
    posterior_means = mu_samples.mean(dim=("chain", "draw")).values 

    mu_values = np.array(mu_values)


    errors = posterior_means - mu_values
    ME = np.mean(errors)
    RMSE = np.sqrt(np.mean(errors**2))

    num_patients = len(mu_values)

    plt.figure(figsize=(10, 7))

    for i, (lower, upper, mu_value, posterior_mean) in enumerate(zip(lower_bounds, upper_bounds, mu_values, posterior_means), start=1):
        plt.plot([lower, upper], [i, i], 'b-', alpha=0.7, linewidth=2, label="95% CI" if i == 1 else None) 
        plt.plot(mu_value, i, 'rx', markersize=9.5, label="True methylation rate" if i == 1 else None)  
        plt.plot(posterior_mean, i, 'ko', markersize=6.5, label="Posterior Mean" if i == 1 else None)  

    plt.xticks(fontsize=12)
    plt.yticks(range(1, num_patients + 1), fontsize=12) 
    plt.xlabel("Methylation rate (per year)", fontsize=14.5)
    plt.ylabel("Patient Index", fontsize=14.5)

    legend_text = [
        "95% CI",
        "True methylation rate",
        "Posterior Mean"
    ]

    plt.legend(legend_text, fontsize=13, loc='lower right', bbox_to_anchor=(1, 0.2))
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_name, dpi=300, bbox_inches="tight")

    plt.show()




def mu_inf_vs_evo(df_epirate_mu):

    df_unique = df_epirate_mu.drop_duplicates(subset=['Sample_id']).reset_index(drop=True)

    df_unique = df_unique.sort_values(by='t_mean', ascending=True).reset_index(drop=True)

    indices = np.arange(len(df_unique))
    epi_rates_unique = df_unique['epiRate']
    mu_means = df_unique['gamma']
    mu_lower = df_unique['gamma_lower']
    mu_upper = df_unique['gamma_upper']
    num_kappa_sorted = df_unique['t_mean']

    mu_error_lower = mu_means - mu_lower
    mu_error_upper = mu_upper - mu_means

    plt.figure(figsize=(8, 6))

    plt.errorbar(
        indices, mu_means, 
        yerr=[mu_error_lower, mu_error_upper], 
        fmt='o', color='blue', alpha=0.5, capsize=3, label='Inferred mu'
    )

    plt.scatter(indices, epi_rates_unique, color='red', marker='x', s=80, label='True mu')

    plt.xlabel('epiRate',fontsize=13)
    plt.ylabel('mu (Inferred and True)',fontsize=13)
    plt.title('Comparison of Inferred and True mu Values with Error Bars (Sorted by num_kappa)')

    plt.xticks(indices, num_kappa_sorted) 
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid(True)
    plt.legend(fontsize=11)

    plt.show()


def mu_diff_regplot(fig_name, df_epirate_mu,num_sites):

    mu_diff = df_epirate_mu['mu'] - df_epirate_mu['epiRate']
    
    mu_error_lower = df_epirate_mu['mu'] - df_epirate_mu['mu_lower']
    mu_error_upper = df_epirate_mu['mu_upper'] - df_epirate_mu['mu']

    plt.figure(figsize=(12, 9))
    
    ax = sns.regplot(
        x=num_sites, 
        y=abs(mu_diff), 
        scatter_kws={'s': 50, 'color': 'blue', 'label': 'Data points'}, 
        line_kws={'color': 'black', 'label': 'Linear fit'}
    )

    plt.errorbar(
        num_sites, mu_diff, 
        yerr=[mu_error_lower, mu_error_upper], 
        fmt='o', color='blue', alpha=0.5, capsize=3, label='95% CI'
    )

    slope, intercept, r_value, p_value, std_err = linregress(num_sites, mu_diff)
    
    print(f"Slope: {slope}")
    print(f"Intercept: {intercept}")
    print(f"R-squared: {r_value**2}")
    print(f"P-value: {p_value}")
    print(f"Standard Error: {std_err}")

    plt.axhline(0, color='red', linestyle='--', linewidth=2, label='No difference')
    plt.xlabel('Number of Sites', fontsize=14.5)
    plt.ylabel('Difference (mu - epiRate)', fontsize=14.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.grid(True)
    plt.legend(fontsize=13)
    plt.tight_layout()
    plt.savefig(f'{fig_name}', dpi=300, bbox_inches="tight")
    plt.show()



def mu_regplot(fig_name, df_epirate_mu):

    epi_rates_unique = df_epirate_mu['epiRate']
    mu_means = df_epirate_mu['mu']
    mu_lower = df_epirate_mu['mu_lower']
    mu_upper = df_epirate_mu['mu_upper']

    mu_error_lower = mu_means - mu_lower
    mu_error_upper = mu_upper - mu_means

    plt.figure(figsize=(12, 9))
    ax = sns.regplot(
        x=epi_rates_unique, 
        y=mu_means, 
        scatter_kws={'s': 50, 'color': 'blue', 'label': 'Data points'}, 
        line_kws={'color': 'black', 'label': 'Linear fit'}
    )

    plt.errorbar(
        epi_rates_unique, mu_means, 
        yerr=[mu_error_lower, mu_error_upper], 
        fmt='o', color='blue', alpha=0.5, capsize=3, label='95% CI'
    )

    slope, intercept, r_value, p_value, std_err = linregress(epi_rates_unique, mu_means)
    
    print(f"Slope: {slope}")
    print(f"Intercept: {intercept}")
    print(f"R-squared: {r_value**2}")
    print(f"P-value: {p_value}")
    print(f"Standard Error: {std_err}")

    x_vals = np.linspace(0, 0.02, 100)
    plt.plot(x_vals, x_vals, 'r--', linewidth=2, label='y = x')

    plt.xlim(0.006, 0.0205)
    plt.xlabel('True switching rate (per year)', fontsize=14.5)
    plt.ylabel('Inferred methylation rate (per year)', fontsize=14.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.grid(True)
    plt.legend(fontsize=13)
    plt.tight_layout()
    plt.savefig(f'{fig_name}', dpi=300, bbox_inches="tight")
    plt.show()


def chr_eval(fig_name, df, x_col='chr', y_col='t_mean'):
    plt.figure(figsize=(8, 6))
    sns.boxplot(x=x_col, y=y_col, data=df)
    sns.stripplot(x=x_col, y=y_col, data=df, color='black', alpha=0.5)
    plt.xlabel('Chromosome Number',fontsize=14.5)
    plt.ylabel('Inferred event time',fontsize=14.5)
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12) 
    plt.tight_layout()
    plt.savefig(f'{fig_name}', format='pdf', dpi=300)
    plt.show()

