import cmdstanpy
import numpy as np
import pandas as pd
from scipy.special import softmax
import pickle
import arviz as az
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from plot import plot_event_time_uncertainty_hier_switched, chr_eval, mu_regplot, mu_inf_vs_evo, plot_event_time_uncertainty, plot_ppc, plot_pair

tri_az = '/Volumes/T7/ICR1/real/tri_real(4.3)kap3.4,0.2_filt.nc'
tri_sum = '/Volumes/T7/ICR1/real/tri_real(4.3)kap3.4,0.2_filt.pkl'
type = 2


if type == 2:
    az_fit = az.from_netcdf(tri_az)
   
    with open(tri_sum, 'rb') as f:
        summary_df = pickle.load(f)

    with open('dict_tri_cll.pkl', 'rb') as f:
        dict_tri = pickle.load(f)
    selected_dict = dict_tri 

all_beta_vals = [] 

counter=0

for patient_id, patient_data in selected_dict.items():
    chr_data = patient_data.get('chr_data', {})
    for chr_id, chr_info in chr_data.items():
        beta_vals = chr_info.get('beta', [])
        all_beta_vals.extend(beta_vals)
        print('sss',len(beta_vals))
        counter+=1

# density = gaussian_kde(all_beta_vals, bw_method=0.12)  
# x_vals = np.linspace(0, 1, 1000)  
# y_vals = density(x_vals)

# plt.hist(all_beta_vals, bins=50, density=True, alpha=0.5, label="Histogram")
# plt.plot(x_vals, y_vals, label="KDE Fit", linewidth=2, color='red')

# plt.xlabel("Fraction of alleles methylated", fontsize=14.5)
# plt.ylabel("Probability Density", fontsize=14.5)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# # plt.legend(loc='upper right', bbox_to_anchor=(0.8,1))
# plt.tight_layout()
# plt.savefig('cll_hist', dpi=300, bbox_inches="tight")
# plt.show()
# eps = 1e-6
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import gaussian_kde





# posterior_pred = az_fit.posterior_predictive["y_rep"]

# observed_data = {"y": np.array(all_beta_vals)}  
# posterior_predictive = {"y_rep": az_fit.posterior_predictive["y_rep"].values}
# new_az_data = az.from_dict(observed_data=observed_data, posterior_predictive=posterior_predictive)
# az.plot_ppc(new_az_data, data_pairs={"y": "y_rep"}, kind="kde", figsize=(8, 6))
# plt.show()



# num_sites = []
# vals = []
# patient_ages = []
# patient_idx = []
# cna_per_patient = []

# def count_nested_keys(data):
#     return {sample: len(sample_data) for sample, sample_data in data.items()}

# nested_key_counts = count_nested_keys(selected_dict)

# for count, value in enumerate(nested_key_counts.values(), 1):
#     cna_per_patient.append(value)
#     for i in range(value):
#         patient_idx.append(count)

# for patient, chr_data in selected_dict.items():
#     for chr_values in chr_data.values():
#         beta_values = chr_values['beta']
#         num_sites.append(len(beta_values))
#         vals.extend(beta_values)
#         patient_ages.append(chr_values['age'])

# with open('/Users/finnkane/Desktop/ICR/dict_tri.pkl', 'rb') as f:
#         d = pickle.load(f)
# df_csv = pd.read_csv('CLL-fCpGs-by-CNA-v2.csv')
# df_csv.columns = df_csv.columns.str.strip()
# print(df_csv.columns)
ages_list = []

for patient_id, patient_data in dict_tri.items():
    for chr_id, chr_data in patient_data['chr_data'].items():
        ages_list.append(chr_data['age'])
print(ages_list)



unique_chr_counts = {patient_id: data['num_unique_chr'] for patient_id, data in dict_tri.items()}

for patient_id, num_unique_chr in unique_chr_counts.items():
    print(f"Patient {patient_id}: {num_unique_chr} unique chromosomes")



kappa_values = summary_df[summary_df.index.str.contains(r'^kappa\[\d+\]$')]['Mean'].values

sample_chr_counts = list(unique_chr_counts.items())   # e.g. [(sample1, nchr1), (sample2, nchr2), ...]
split_kappa_values = np.split(
    kappa_values,
    np.cumsum([count for _, count in sample_chr_counts])[:-1]
)


epi_rates = [
    chr_data['epiRate']
    for sample_id, sample_data in dict_tri.items()
    if isinstance(sample_data, dict) and 'chr_data' in sample_data
    for chr_data in sample_data['chr_data'].values()
    if isinstance(chr_data, dict) and 'epiRate' in chr_data
]

epi_rates_unique = list(dict.fromkeys(epi_rates))

sample_ids = list(selected_dict.keys())
epi_rates_dict = {
    sample_id: epi_rates_unique[i]
    for i, sample_id in enumerate(sample_ids)
}


mu_summary = summary_df[summary_df.index.str.contains(r'^mu\[\d+\]$')].copy()
mu_summary = mu_summary[['Mean', '5%', '95%']]

filtered_mu_summary = mu_summary[mu_summary['Mean'] < 0.075]

print(filtered_mu_summary)


mu_summary['Sample_id'] = sample_ids[:len(mu_summary)]
mu_dict = mu_summary.set_index('Sample_id')[['Mean', '5%', '95%']].to_dict(orient='index')


gamma_summary = summary_df[summary_df.index.str.contains(r'^gamma\[\d+\]$')].copy()
gamma_summary = gamma_summary[['Mean', '5%', '95%']]  

gamma_summary['Sample_id'] = sample_ids[:len(gamma_summary)]
gamma_dict = gamma_summary.set_index('Sample_id')[['Mean', '5%', '95%']].to_dict(orient='index')


final_dict = {}

for (sample_id, nchr), kappa_vals in zip(sample_chr_counts, split_kappa_values):
    sample_data = dict_tri[sample_id]
    chr_data_dict = sample_data['chr_data'] 


    chr_entries = list(chr_data_dict.values())

    per_chr_info = []
    for kappa_val, chr_info in zip(kappa_vals, chr_entries):
        per_chr_info.append({
            'chr':       chr_info['chr'],
            'kappa':     kappa_val,
            'epiRate':   chr_info['epiRate'],
            'age':       chr_info['age'],
            'num_sites': chr_info['num_sites'],
            'cancerAge': chr_info['cancerAge'],
            'restrictAge': chr_info['restrictAge']
        })

    final_dict[sample_id] = {
        'gamma':    gamma_dict.get(sample_id, {}).get('Mean', None),
        'gamma_lower': mu_dict.get(sample_id, {}).get('5%', None),
        'gamma_upper': mu_dict.get(sample_id, {}).get('95%', None),
        'mu':       mu_dict.get(sample_id, {}).get('Mean', None),
        'mu_lower': mu_dict.get(sample_id, {}).get('5%', None),
        'mu_upper': mu_dict.get(sample_id, {}).get('95%', None),
        'epiRate':  epi_rates_dict.get(sample_id, None),  # sample-level epiRate
        'chr_data': per_chr_info

    }

rows = []
sample_idx = 0

for sample_id, sample_info in final_dict.items():
    sample_idx += 1  
    for chr_info in sample_info['chr_data']:
        rows.append({
            'Sample_id':  sample_id,
            'index':      sample_idx,
            'chr':        chr_info['chr'],
            'num_kappa': chr_info['kappa'],
            'epiRate': chr_info['epiRate'],   
            'gamma':      sample_info['gamma'],
            'gamma_lower': mu_dict.get(sample_id, {}).get('5%', None),
            'gamma_upper': mu_dict.get(sample_id, {}).get('95%', None),
            'mu':         sample_info['mu'],
            'mu_lower':   sample_info['mu_lower'],
            'mu_upper':   sample_info['mu_upper'],
            'diff':     chr_info['epiRate'] - sample_info['mu']
        })

df_final = pd.DataFrame(rows)
df_final = df_final.sort_values(by='diff', ascending=False)
print(df_final.head(10))

# print(df_final[df_final['mu'] < 0.0075])

patient_errors = {}

unique_true_mu = []
unique_inferred_mu = []


for sample_id, group in df_final.groupby('Sample_id'):
    inferred_mu = group['mu'].unique()[0] 
    
    true_mu = group['mu'].unique()[0] 
    
    unique_true_mu.append(true_mu)
    unique_inferred_mu.append(inferred_mu)

absolute_errors = [abs(inferred - true) for inferred, true in zip(unique_inferred_mu, unique_true_mu)]
mae = sum(absolute_errors) / len(absolute_errors) if absolute_errors else 0

print(f"Mean Absolute Error (MAE) Across Patients: {mae:.6f}")


t_summary = summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]
eta = summary_df[summary_df.index.str.contains(r'^eta\[\d+\]$')]['Mean'].values
delta  = summary_df[summary_df.index.str.contains(r'^delta\[\d+\]$')]['Mean'].values
t = t_summary['Mean']
print(t_summary)
df_final['eta'] = eta
df_final['delta'] = delta
df_final['age'] = ages_list
df_final['t_mean'] = list(t)
df_final['bounds'] = list(t_summary['StdDev'])

# chr12_df = df_final[df_final['chr'] == 'chr12']
# print(chr12_df)
# print(len(list(t)))
# chr_eval(df_final)

az.plot_pair(
    az_fit,
    var_names=["gamma", "mu", "t"],  
    coords={"gamma_dim_0": [3], "mu_dim_0": [3], "t_dim_0": [4]}, 
)
t = t_summary['Mean']

df_final['eta'] = eta
df_final['delta'] = delta
df_final['age'] = ages_list
df_final['t_mean'] = list(t)
event_times = df_final['t_mean']
patient_ages = df_final['age']
patient_idx = df_final['index']
chr = df_final['chr']

plot_event_time_uncertainty_hier_switched('inf', az_fit, event_times, patient_ages, patient_idx, chr)

plt.show()
# mu_linear_regression(df_final)
mu_regplot(df_final)
mu_inf_vs_evo(df_final)
                                                                                                                                                                                                                                                                                                  
t_summary = summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]
mu_summary = summary_df[summary_df.index.str.contains(r'^mu\[\d+\]$')]
gamma_summary = summary_df[summary_df.index.str.contains(r'^gamma\[\d+\]$')]
kappa_summary = summary_df[summary_df.index.str.contains(r'^kappa\[\d+\]$')]
mu_sdd = summary_df[summary_df.index.str.contains(r'mu_sdd')]
print(mu_summary)